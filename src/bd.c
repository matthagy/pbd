/* -*- Mode: c -*-
 * bd.c - Brownian Dynamics Integrator
 *--------------------------------------------------------------------------
 * Copyright (C) 2009, Matthew Hagy (hagy@gatech.edu)
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* Parallel Simulation Algorithm
 * -------------------------------------------------------------------
 * The parallel simulation algorithm can best be described as follows.
 * 
 *   o Initially allocate particles into cells
 *
 *   o Build Neighbor Lists
 *       Here we figure out which pairs of particles must be 
 *       considered in force evaluations.  This involves communication
 *       between threads to figure out which particles in one cell are
 *       neighbors in junctioned cells.
 *
 *   o Simulation Cycles
 *       For each time step we send the current coordinates of these
 *       neighboring particles between threads.  The neighbor lists are
 *       then used to evaluate the forces on each particle, and from this
 *       we perform an integration with respect to time.  After each 
 *       integration we check if particles have moved far enough for the
 *       neighbor lists to be invalid, and once they are we rebuild them.
 */

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>

#include "opt.h"
#include "debug.h"
#include "vector.h"
#include "random.h"
#include "mem.h"
#include "array.h"
#include "comm.h"
#include "periodic.h"
#include "cells.h"
#include "constants.h"
#include "bd.h"
#include "init.h"

/* Data from bd.h
 *---------------*/
double CEX_fric_gamma=-1;
double CEX_T=-1;
double CEX_kT=-1;
double CEX_dt=-1;
double CEX_dU_max=-1;
double CEX_r_pair_cutoff=-1;

linterp_table_t CEX_pair_potential = {-1,-1, NULL};
linterp_table_t CEX_pair_force = {-1,-1, NULL};
int CEX_force_update_rate=-1;

double CEX_r_neighbor=-1;
double CEX_r_neighbor_sqr=-1;

array_t *CEX_positions=NULL;
array_t *CEX_new_positions=NULL;
int CEX_N_internal_particles=-1;
array_t *CEX_tags=NULL;
array_t *CEX_forces=NULL;
array_t *CEX_random_vectors=NULL;
array_t *CEX_nl_displace=NULL;

array_t *CEX_internal_neighbors=NULL;
array_t *CEX_external_neighbors=NULL;


/* Building Neighbor Lists
 * -----------------------------------------------------------
 * The building of neighbor lists is a multistep process that 
 * requires communication between threads.  This process can be 
 * described as follows:
 *
 *  o Update Particle Membership
 *     Particles may have moved outside of this threads's cell 
 *     since the last time neighbor lists where built. We find 
 *     these particles, and send them to the proper thread 
 *     based on their positions and the extent of our junctioned 
 *     cells.
 *
 *  o Determine Particles Close to Junctioned Cells
 *     For each junctioning cell mode, determine which particles in
 *     this cell have to be considered by the junctioned cell. 
 *     This information is then sent to the corresponding thread.  
 * 
 *  o Build Neighbor Lists
 *     Use the location of internal particles and particles of junctioned cells
 *     to build neighbor lists for the particles internal to this cell.
 *
 *  o Remove Un-needed External Particles
 *     Find out which external particles are not used in any force
 *     evaluations, as they are not included in the neighbor lists.
 *     We requests these particles be removed from the send indices
 *     table to minimize communication overhead.
 */

static void update_particle_membership(void);
static void determine_possible_neighbors(void);
static void rebuild_neighborlists(void);
static void remove_unneeded_external_particles(void);
static void allocate_external_exchange_buffers(void);
static void sort_neighbor_list(array_t *);
static void sort_send_indices(void);
static void setup_force_aux(void);

void
CEX_thread_update_neighbors(void)
{
        REQ_INIT();
        if (HAVE_JUNCTIONS()) {
                update_particle_membership();
                determine_possible_neighbors();
        }
        rebuild_neighborlists();
        sort_neighbor_list(CEX_internal_neighbors);
        if (HAVE_JUNCTIONS()) {
                remove_unneeded_external_particles();
                allocate_external_exchange_buffers();
                sort_neighbor_list(CEX_external_neighbors);
                sort_send_indices();
        }
        setup_force_aux();
        /* don't clear CEX_send_indices, as we'll uses these 
         * indices durring simulation to communicate new positions
         * of external particles durring simulation */
}


/* For both redistributing particles and for exchaining neighboring positions, 
 * we have to associate with each communicator an array of particle indices 
 * internal to this thread to send to the corresponding thread. */
array_t *CEX_send_indices=NULL;
#define GET_SEND_INDICES(comm) GET_COMM_DATA(array_t *,  CEX_send_indices,  comm)

/* On the recprical end, we need to know the length of array's recieved */
array_t *CEX_recv_lengths=NULL;
#define GET_RECV_LENGTH(comm) GET_COMM_DATA(int,  CEX_recv_lengths,  comm)
#define SET_RECV_LENGTH(comm,len) (GET_RECV_LENGTH(comm) = (len))

static void clear_send_indices(void);
static void exchange_send_lengths(void);


/* Updating Particle Membership
 *-------------------------------------------------------------
 * Redistribute exited particles to adjacent cells, s.t. that 
 * each particle is within the extent of its new cell.
 */

/* add indices of particle to send to each comm's send_indices */
static int annote_exited_particles(void);
/* temporary arrays to recieve positions and tags */
static int setup_exchange_exit_temps(void);
/* send/recieve particles from pair communicators */
static void exchange_exited_particles(void);
/* update internal data structures */
static void remove_exited_particles(int n_sent);
static void insert_entered_particles(void);

static void
update_particle_membership(void)
{
        /* lose all concept of external particles */
        truncate_array(CEX_positions, CEX_N_internal_particles);
        clear_send_indices();
        int n_sent = annote_exited_particles();
        exchange_send_lengths();
        int n_recv = setup_exchange_exit_temps();
        exchange_exited_particles();
        remove_exited_particles(n_sent);
        insert_entered_particles();
        clear_send_indices();
        /* update internal data structures */
        CEX_N_internal_particles += n_recv - n_sent;
        assert(ARR_LENGTH(CEX_positions) == CEX_N_internal_particles);
        CEX_prealloc_array(CEX_forces, CEX_N_internal_particles);
        ARR_LENGTH(CEX_forces) = CEX_N_internal_particles;
        CEX_prealloc_array(CEX_nl_displace, CEX_N_internal_particles);
        ARR_LENGTH(CEX_nl_displace) = CEX_N_internal_particles;
        CEX_prealloc_array(CEX_new_positions, CEX_N_internal_particles);
        ARR_LENGTH(CEX_new_positions) = CEX_N_internal_particles;
        CEX_prealloc_array(CEX_random_vectors, CEX_N_internal_particles);
        ARR_LENGTH(CEX_random_vectors) = CEX_N_internal_particles;
}

static int 
annote_exited_particles(void)
{
        int inx, n_sent=0;
        for (inx=0; inx<CEX_N_internal_particles; inx++) {
                vec_t position = ARR_INDEX_AS(vec_t, CEX_positions,inx);
                if (unlikely(!CELL_CONTAINS(CEX_this_cell, position))) {
                        cell_t *dest_cell = find_jcell_containing(position);
                        if (unlikely(dest_cell==NULL)) {
                                Fatal("lost particle " Vec3E_FRMT
                                      "in cell " Vec3E_FRMT " to " Vec3E_FRMT,
                                      Vec3_ARGS(position),
                                      Vec3_ARGS(CEX_this_cell->min_extent),
                                      Vec3_ARGS(CEX_this_cell->max_extent));
                        }
                        array_t *send_indices = GET_SEND_INDICES(dest_cell->comm);
                        IARR_APPEND(send_indices, inx);
                        n_sent ++;
                }
        }
        return n_sent;
}

/* temporary arrays for each communicator to recieve particles */
array_t *CEX_tmp_recv_positions=NULL;
array_t *CEX_tmp_recv_tags=NULL;

#define GET_TMP_RECV_POSITIONS(comm) \
        GET_COMM_DATA(array_t *,  CEX_tmp_recv_positions,  comm)
#define GET_TMP_RECV_TAGS(comm) \
        GET_COMM_DATA(array_t *,  CEX_tmp_recv_tags,  comm)

static inline void
setup_temp_array(int size, array_t *arr)
{
        clear_array(arr);
        CEX_prealloc_array(arr, size);
}

static int
setup_exchange_exit_temps(void)
{
        comm_t *comm;
        int counter, n_recv=0;
        
        COMM_FOREACH(comm, counter) {
                int size = GET_RECV_LENGTH(comm);
                n_recv += size;
                setup_temp_array(size, GET_TMP_RECV_POSITIONS(comm));
                setup_temp_array(size, GET_TMP_RECV_TAGS(comm));
        }
        return n_recv;
}

static void 
exchange_exited_particles(void)
{
        DO_COMM(comm, /* positions */
        /* send */ comm_send_vecs_by_index(comm, CEX_positions, GET_SEND_INDICES(comm)),
        /* recv */ comm_recv_vecs(comm, GET_TMP_RECV_POSITIONS(comm)));
        DO_COMM(comm, /* tags */
        /* send */ comm_send_ints_by_index(comm, CEX_tags, GET_SEND_INDICES(comm)),
        /* recv */ comm_recv_ints(comm, GET_TMP_RECV_TAGS(comm)));
}

static void 
remove_exited_particles(int n_sent)
{
        /* this could be more efficient, don't really care */
        array_t *all_sent = CEX_make_int_array(n_sent);
        comm_t *comm;
        int counter;
        COMM_FOREACH(comm, counter) {
                CEX_extend_array(all_sent, GET_SEND_INDICES(comm));
        }
        CEX_remove_array_of_indices(CEX_positions, all_sent);
        CEX_remove_array_of_indices(CEX_tags, all_sent);
        CEX_free_array(all_sent);
}

static void 
insert_entered_particles(void)
{
        comm_t *comm;
        int counter;
        
        COMM_FOREACH(comm, counter) {
                CEX_extend_array(CEX_positions, GET_TMP_RECV_POSITIONS(comm));
                CEX_extend_array(CEX_tags, GET_TMP_RECV_TAGS(comm));
        }
}

/* Send Indices */
static void
clear_send_indices(void)
{
        comm_t *comm;
        int counter;
        COMM_FOREACH(comm, counter) {
                clear_array(GET_SEND_INDICES(comm));
        }
}

static void
exchange_send_lengths(void)
{
        DO_COMM(comm, 
        /* send */ comm_send_int(comm, (int)ARR_LENGTH(GET_SEND_INDICES(comm))),
        /* recv */ SET_RECV_LENGTH(comm, comm_recv_int(comm)));
}

static void
sort_send_indices(void)
{
        comm_t *comm;
        int counter;
        COMM_FOREACH(comm, counter) {
                CEX_sort_int_array(GET_SEND_INDICES(comm));
        }
}


/* Determine Particles Close to Junctioned Cells
 *-------------------------------------------------------------*/

/* For every particle, check which junctions the particle is close enough
 * to, to possibly be a neighbor to a particle inside of that adjacent cell.
 * We add the index of these particles to the send indices of the 
 * communicators for the threads corresponding to these adjacent cells. */
static void record_particles_near_junctions(void);

/* Once we know how many positions we'll be recieving from each of
 * our adjacnet threads, allocte space in CEX_positions to hold 
 * these positions  */
static void allocate_for_external_neighbors(void);

/* Following allocation, we exchange these positions */
static void exchange_external_positions(void);

static void 
determine_possible_neighbors(void)
{
        clear_send_indices();
        record_particles_near_junctions();
        exchange_send_lengths();
        allocate_for_external_neighbors();
        exchange_external_positions();
}

/* record all cells that this position is close enough to, to be a neighbor */
static inline void check_position_junctions(vec_t position, array_t *cells_near)
        GCC_ATTRIBUTE((always_inline));

static void
record_particles_near_junctions(void)
{
        /* all particles should be internal at this point */
        assert(ARR_LENGTH(CEX_positions)==CEX_N_internal_particles);
        vec_t *posp;
        int pos_counter;
        array_t *cells_near = CEX_make_array(sizeof(void *), ARR_LENGTH(CEX_jcells));
        VARR_FOREACH(CEX_positions, posp, pos_counter) {
                clear_array(cells_near);
                check_position_junctions(*posp, cells_near);
                if (unlikely(ARR_LENGTH(cells_near))) {
                        int pos_index = (int)(posp - ARR_DATA_AS(vec_t, CEX_positions));
                        cell_t **cellp; int cell_counter;
                        XARR_FOREACH(cells_near, cellp, cell_counter) {
                                #if 0
                                xprintf(Vec3_FRMT("%.2f") " [%d] near %d", 
                                        Vec3_ARGS_SCALED(1e9, *posp),
                                        pos_index, (*cellp)->comm->comm_rank);
                                #endif
                                array_t *send_indices = GET_SEND_INDICES((*cellp)->comm);
                                IARR_APPEND(send_indices, pos_index);
                        }
                }
        }
        CEX_free_array(cells_near);
}

static void 
allocate_for_external_neighbors(void)
{
        comm_t *comm;
        int counter, n_recv = 0;
        COMM_FOREACH(comm, counter) {
                n_recv += GET_RECV_LENGTH(comm);
        }
        CEX_prealloc_array(CEX_positions, ARR_LENGTH(CEX_positions) + n_recv);
}

/* For every sequence of external positions recieved from a neighboring cell
 * we need to know where these positions are stored in CEX_positions for both
 * re-recieving position for evaluating forces and for removing unneeded 
 * external indices. */
array_t *CEX_ext_positions_offset=NULL;

#define GET_EXT_POSITIONS_OFFSET(comm) \
        GET_COMM_DATA(int,  CEX_ext_positions_offset,  comm)
#define SET_EXT_POSITIONS_OFFSET(comm,index) \
        (GET_EXT_POSITIONS_OFFSET(comm) = (index))

static void 
exchange_external_positions(void)
{
        DO_COMM(comm, 
        /* send */ comm_send_vecs_by_index(comm, CEX_positions, 
                                           GET_SEND_INDICES(comm)),
        /* recv */ ({SET_EXT_POSITIONS_OFFSET(comm, ARR_LENGTH(CEX_positions));
                     comm_recv_extend_vecs(comm,  CEX_positions);}));
}


/* for 1 position, check every junction, and record nearby ones */
static inline void add_uniq(array_t *cells_near, cell_t* cell) 
        GCC_ATTRIBUTE((always_inline));
static inline int check_surface_junction(surface_junction_t *, vec_t)
        GCC_ATTRIBUTE((always_inline));
static inline int check_line_junction(line_junction_t *, vec_t)
        GCC_ATTRIBUTE((always_inline));
static inline int check_point_junction(point_junction_t *, vec_t)
        GCC_ATTRIBUTE((always_inline));

static inline void
check_position_junctions(vec_t position, array_t *cells_near)
{
        int counter;
        surface_junction_t *sjp;
        XARR_FOREACH(CEX_surface_junctions, sjp, counter) {
                if (unlikely(check_surface_junction(sjp, position))) {
                        add_uniq(cells_near, sjp->cell);
                }
        }
        line_junction_t *ljp;
        XARR_FOREACH(CEX_line_junctions, ljp, counter) {
                if (unlikely(check_line_junction(ljp, position))) {
                        add_uniq(cells_near, ljp->cell);
                }
        }
        point_junction_t *pjp;
        XARR_FOREACH(CEX_point_junctions, pjp, counter) {
                if (unlikely(check_point_junction(pjp, position))) {
                        add_uniq(cells_near, pjp->cell);
                }
        }
}

static inline void
add_uniq(array_t *cells_near, cell_t* cell)
{
        cell_t **cellp; int cell_counter;
        XARR_FOREACH(cells_near, cellp, cell_counter) {
                if (unlikely(*cellp == cell)) {
                        return;
                }
        }
        ARR_APPEND(cell_t *, cells_near, cell);
}

static inline int 
check_surface_junction(surface_junction_t *sjp, vec_t pos)
{
        vec_t extent, r;
        extent = sjp->dir==-1 ? CEX_this_cell->min_extent : CEX_this_cell->max_extent;
        PERIODIC_SEPARATION_VECTOR(r, extent, pos);
        return fabs(INDEX_AXIS(&r, sjp->axis)) <= CEX_r_neighbor;
}

static inline int 
check_line_junction(line_junction_t *lnp, vec_t pos)
{
        vec_t point, r;

        INDEX_AXIS(&pos, lnp->axis) = 0;
        switch (lnp->axis) {
        case AXIS_X: 
                Vec3_SET(point, 0, lnp->offset1, lnp->offset2);
                break;
        case AXIS_Y:
                Vec3_SET(point, lnp->offset1, 0, lnp->offset2);
                break;
        case AXIS_Z:
                Vec3_SET(point, lnp->offset1, lnp->offset2, 0);
                break;
        default:
                Fatal("bad axis %c", lnp->axis);
        }
        PERIODIC_SEPARATION_VECTOR(r, point, pos);
        return Vec3_SQR(r) <= CEX_r_neighbor_sqr;
}

static inline int 
check_point_junction(point_junction_t *pnp, vec_t position)
{
        vec_t r;

        PERIODIC_SEPARATION_VECTOR(r, pnp->offset, position);
        return Vec3_SQR(r) <= CEX_r_neighbor_sqr;
}


/* Build Neighbor Lists
 *-------------------------------------------------------------
 * Since each cell should only have one the order of 100s of 
 * particles, we naively iterate through every pair of particles
 * and record those that are within the neighbor distance.
 */

/* square maximum absolute distance any single particle can move
 * before neighbor lists are invalid.  This is tracked by
 * CEX_nl_displace */
static double r_delta_2_sqr;

static void rebuild_internal_neighborlists(void);
/* this function is called twice on each neighbor rebuilding cycle.
 * first with all possible external positions, and then
 * again after we've removed the external positions that aren't need  */
static void rebuild_external_neighborlists(void);

static void 
rebuild_neighborlists(void)
{
        double r_delta_2 = (CEX_r_neighbor - CEX_r_pair_cutoff) / 2;
        r_delta_2_sqr = r_delta_2 * r_delta_2;
        assert(ARR_LENGTH(CEX_nl_displace) == CEX_N_internal_particles);
        CEX_zero_array_elements(CEX_nl_displace);
        rebuild_internal_neighborlists();
        rebuild_external_neighborlists();
}

static void
rebuild_internal_neighborlists(void)
{
        clear_array(CEX_internal_neighbors);
        vec_t * CEX_RESTRICT positions = ARR_DATA_AS(vec_t, CEX_positions);
        for (int i=CEX_N_internal_particles - 1; i>=0; i--) {
                vec_t pos_i = positions[i];
                for (int j=i-1; j>=0; j--) {
                        vec_t r, pos_j = positions[j];
                        PERIODIC_SEPARATION_VECTOR(r, pos_i, pos_j);
                        if (Vec3_SQR(r) <= CEX_r_neighbor_sqr) {
                                IARR_APPEND(CEX_internal_neighbors, i);
                                IARR_APPEND(CEX_internal_neighbors, j);
                        }
                }
        }
}

static void
rebuild_external_neighborlists(void)
{
        clear_array(CEX_external_neighbors);
        vec_t * CEX_RESTRICT positions = ARR_DATA_AS(vec_t, CEX_positions);
        /* i: loop over external particles
         * j: loop over interanl particles */
        for (int i=ARR_LENGTH(CEX_positions) - 1; i>=CEX_N_internal_particles; i--) {
                vec_t pos_i = positions[i];
                for (int j=CEX_N_internal_particles-1; j>=0; j--) {
                        vec_t r, pos_j = positions[j];
                        PERIODIC_SEPARATION_VECTOR(r, pos_i, pos_j);
                        if (Vec3_SQR(r) <= CEX_r_neighbor_sqr) {
                                /* add internal first */
                                IARR_APPEND(CEX_external_neighbors, j);
                                IARR_APPEND(CEX_external_neighbors, i);
                        }
                }
        }
}

static int
cmp_neighbors(int *a, int *b)
{
        int x = a[0] - b[0];
        if (unlikely(x==0)) {
                x = a[1] - b[1];
        }
        return x;
}

/* sort neighbor lists to imporve memory locality to
 * imporve  L2 cache performance */
static void 
sort_neighbor_list(array_t *nl)
{
        assert(ARR_LENGTH(nl) % 2 == 0);
        qsort(ARR_DATA(nl), ARR_LENGTH(nl)>>1, 2*sizeof(int),
              (array_el_comparer)&cmp_neighbors);
}


/* Remove Un-needed External Particles
 *-------------------------------------------------------------
 * to remove external particles from those that we will recieve; 
 * we record the index in the CEX_send_indices that we no longer
 * wish to recieve */

static inline void clear_remove_indices(void);
static void find_remove_indices(void);
static void exchange_removes_indices(void);

static void 
remove_unneeded_external_particles(void)
{
        clear_remove_indices();
        find_remove_indices();
        exchange_removes_indices();
        /* remove all external particles and repeat end of 
         * determine_possible_neighbors after the indices
         * of external particles we don't need are removed
         * from CEX_send_indices */
        ARR_LENGTH(CEX_positions) = CEX_N_internal_particles;
        exchange_send_lengths();
        /* allocate shouldn't be necessary */
        exchange_external_positions();
        /* rebuild external neighbor lists using new indices
         * of external particles */
        rebuild_external_neighborlists();
        clear_remove_indices();
}

array_t *CEX_remove_indices=NULL;
#define GET_REMOVE_INDICES(comm) GET_COMM_DATA(array_t *,  CEX_remove_indices,  comm)

static inline void
clear_remove_indices(void)
{
        comm_t *comm; int counter;
        
        COMM_FOREACH(comm, counter) {
                clear_array(GET_REMOVE_INDICES(comm));
        }
}

static void 
find_remove_indices(void)
{
        array_t *uses;
        int n_external;
        int *enp, counter;

        /* record number of uses of particles by their
         * index in CEX_positions */
        n_external = ARR_LENGTH(CEX_positions) - CEX_N_internal_particles;
        uses = CEX_make_int_array(n_external);
        ARR_LENGTH(uses) = n_external;
        CEX_zero_array_elements(uses);
         /* internal index in CEX_external_neighbors is first, skip it using slice */
        XARR_FOREACH_SLICE(CEX_external_neighbors, 
                           1, ARR_LENGTH(CEX_external_neighbors), 2,
                           enp, counter) {
                int external_index = *enp - CEX_N_internal_particles;
                ARR_INDEX_AS(int, uses, external_index) ++;
        }
        /* for every external particle we recieve, check its uses.
         * if its not used, record its send indices in remove_indices */
        comm_t *comm;
        COMM_FOREACH(comm, counter) {
                int external_offset = GET_EXT_POSITIONS_OFFSET(comm) - CEX_N_internal_particles;
                int n_send_indices = GET_RECV_LENGTH(comm);
                array_t *remove_indices = GET_REMOVE_INDICES(comm);
                for (int send_index=0; send_index<n_send_indices; send_index++) {
                        int external_index = external_offset + send_index;
                        assert(external_index < n_external);
                        int n_uses = ARR_INDEX_AS(int, uses, external_index);
                        if (unlikely(n_uses==0)) {
                                IARR_APPEND(remove_indices, send_index);
                        }
                }
        }
        CEX_free_array(uses);
}

static void 
exchange_removes_indices(void)
{
        array_t *remove_indices = CEX_make_int_array(10);
        DO_COMM(comm,
        /* send */  comm_send_ints(comm, GET_REMOVE_INDICES(comm)),
        /* recv */({clear_array(remove_indices);
                    comm_recv_ints(comm, remove_indices);
                    CEX_remove_array_of_indices(GET_SEND_INDICES(comm), remove_indices);}));
        CEX_free_array(remove_indices);
}


/* Updating External Positions
 *----------------------------------------------------------------
 * Update external positions for every force evaluation */

array_t *CEX_send_positions_buffers;

#define GET_SEND_POSITION_BUFFER(comm)                                  \
        GET_COMM_DATA(array_t *,  CEX_send_positions_buffers,  comm)

static void 
allocate_external_exchange_buffers(void)
{
        int counter;
        comm_t *comm;
        COMM_FOREACH(comm, counter) {
                array_t * send_buffer = GET_SEND_POSITION_BUFFER(comm);
                int length = ARR_LENGTH(GET_SEND_INDICES(comm));
                CEX_prealloc_array(send_buffer, length);
                ARR_LENGTH(send_buffer) = length;
        }
}

/* this function is the main bottleneck in parallel applications
 * we therefore use non-blocking IO and synchronize everything 
 * at the end 
 */
static void
update_external_positions(void)
{
        /* copy to external buffers */ {
        int counter;
        comm_t *comm;
        vec_t * CEX_RESTRICT _positions = ARR_DATA_AS(vec_t, CEX_positions);
        COMM_FOREACH(comm, counter) {
                array_t * send_buffer = GET_SEND_POSITION_BUFFER(comm);
                array_t * send_indices = GET_SEND_INDICES(comm);
                for (int i=ARR_LENGTH(send_indices)-1; i>=0; i--) {
                        ARR_INDEX_AS(vec_t, send_buffer, i) =
                                _positions[ARR_INDEX_AS(int, send_indices, i)];
                }
        }}
        /* exchange */
        DO_COMM(comm,
                /* send */
                comm_send_vecs(comm, GET_SEND_POSITION_BUFFER(comm)),
                /* recv */
                comm_recv(comm, ARR_ADDRESS_ELEMENT(CEX_positions,
                                                    GET_EXT_POSITIONS_OFFSET(comm)),
                          3*GET_RECV_LENGTH(comm), MPI_DOUBLE));
}

/* Force Evaluation
 *----------------------------------------------------------------
 * Using neighbor lists, evaluate the forces on each internal 
 * particle.
 */

/* helper macros */

#define _PER_SEP(r, pos_i, pos_j)                                       \
        XPERIODIC_SEPARATION_VECTOR(r, pos_i, pos_j, _box_size, _box_half)

#ifdef SINGLE_POINT_INTERPOLATION
#  define _SETUP_INTERPOLATE_FORCE_LOCALS         \
          double _linterp_x_min = CEX_pair_force.x_min;                   \
          double _inv_linterp_x_prec = 1.0 / CEX_pair_force.x_prec;       \
          const double * CEX_RESTRICT _linterp_table =                    \
                  ARR_DATA_AS(double, CEX_pair_force.table);

#  define _INTERPOLATE_FORCE(r)                                           \
          _linterp_table[(int)((r - _linterp_x_min) * _inv_linterp_x_prec)]
#else

#define _SETUP_INTERPOLATE_FORCE_LOCALS         \
          double _linterp_x_min = CEX_pair_force.x_min;                   \
          double _inv_linterp_x_prec = 1.0 / CEX_pair_force.x_prec;       \
          const double * CEX_RESTRICT _linterp_table =                    \
                  ARR_DATA_AS(double, CEX_pair_force.table);

#  define _INTERPOLATE_FORCE(r)  ({                                       \
           double k = (r - _linterp_x_min) * _inv_linterp_x_prec;         \
           double tblinx = floor(k);                                    \
           int inx = (int)tblinx;                                       \
           double weight = k - tblinx;                                  \
           _linterp_table[inx] * (1 - weight) + _linterp_table[inx+1] * weight; })

#endif

#define _SETUP_FORCE_LOCALS                                             \
        double r_pair_cutoff_sqr = CEX_r_pair_cutoff * CEX_r_pair_cutoff; \
        vec_t _box_size=CEX_box_size, _box_half=CEX_box_half;           \
        vec_t * CEX_RESTRICT _positions=ARR_DATA_AS(vec_t, CEX_positions); \
        vec_t * CEX_RESTRICT _forces=ARR_DATA_AS(vec_t, CEX_forces);    \
        _SETUP_INTERPOLATE_FORCE_LOCALS

static void evaluate_forces(void) GCC_ATTRIBUTE((noinline));


#ifdef OMP_PARALLELIZE_FORCES
# ifdef OMP_CONCURRENT_FORCE_EVALUATION
#   error "can't perform concurrent force evaluation and position exchange"
# endif
# include "eval-forces-openmp.c"
#else
# include "eval-forces-simple.c"
#endif

static inline void update_random(void);

#ifndef OMP_CONCURRENT_FORCE_EVALUATION
static inline void
update_forces(void)
{
        if (HAVE_JUNCTIONS()) {
                update_external_positions();
        }
        evaluate_forces();
}
#else
/* when we have multiple threads and we also have junctions,
 * thread id 0 will communicate external positions to junctioned
 * cells while thread 1 evaluates internal forces.  After this
 * external forces are evaluated.  Otherwise we thread 0 uses 
 * the logic in the above thread*/
static inline void 
update_forces(void)
{

        if (!HAVE_JUNCTIONS()) {
                evaluate_forces();
        } else {
#pragma omp parallel num_threads(2)
                {
                        switch (omp_get_thread_num()) {
                        case 0:
                                update_external_positions();
                                break;
                        case 1:
                                evaluate_internal_forces();
                                update_random();
                                break;
                        default:
                                break;
                        }
                } /* end parallel */
                evaluate_external_forces();
        }
}
#endif /*OMP_CONCURRENT_FORCE_EVALUATION*/

void 
CEX_thread_update_forces(void)
{
        REQ_INIT();
        update_forces();
}



static inline vec_t
eval_one_force(int target_index, vec_t position)
{
        double r_pair_cutoff_sqr = CEX_r_pair_cutoff * CEX_r_pair_cutoff;
        vec_t _box_size=CEX_box_size, _box_half=CEX_box_half;
        vec_t * CEX_RESTRICT _positions=ARR_DATA_AS(vec_t, CEX_positions);
        double _linterp_x_min = CEX_pair_force.x_min;
        double _inv_linterp_x_prec = 1.0 / CEX_pair_force.x_prec;
        const double * CEX_RESTRICT _linterp_table = 
                ARR_DATA_AS(double, CEX_pair_force.table);
        
        vec_t sum_force;
        Vec3_CLEAR(sum_force);
        /* loop over internal neighbors */
        for (int n_counter=ARR_LENGTH(CEX_internal_neighbors) >> 1,
             *n_ptr=ARR_DATA_AS(int, CEX_internal_neighbors);
             n_counter -- > 0;) {
                int part_i = *(n_ptr++);
                int part_j = *(n_ptr++);
                int other_index;
                if (unlikely(part_i == target_index)) {
                        other_index = part_j; 
                } else if (unlikely(part_j == target_index)) {
                        other_index = part_i;
                } else {
                        continue;
                }
                vec_t r;
                _PER_SEP(r, position, _positions[other_index]);
                double rsqr = Vec3_SQR(r);
                if (rsqr<r_pair_cutoff_sqr) {
                        /* XXX Assumes CEX_pair_force has already
                         * been divied by vector length, s.t. this multiplication
                         * also normalizes the force vector to direction unit vector */
                        double force_div_rlen = _INTERPOLATE_FORCE(sqrt(rsqr));
                        vec_t force;
                        Vec3_MUL(force, r, force_div_rlen);
                        Vec3_SUBTO(sum_force, force);
                }
        }
        /* loop over external neighbors */        
        int N_positions GCC_ATTRIBUTE((unused)) = ARR_LENGTH(CEX_positions);
        for (int n_counter=ARR_LENGTH(CEX_external_neighbors) >> 1,
             *n_ptr=ARR_DATA_AS(int, CEX_external_neighbors);
             n_counter -- > 0;) {
                int part_i = *(n_ptr++); /* inner particle */
                int part_j = *(n_ptr++); /* extern particle */
                if (likely(part_i != target_index)) {
                        continue;
                }
                assert(part_i < CEX_N_internal_particles);
                assert(part_j >= CEX_N_internal_particles);
                assert(part_j < N_positions);
                vec_t r;
                _PER_SEP(r, position, _positions[part_j]);
                double rsqr = Vec3_SQR(r);
                if (rsqr<r_pair_cutoff_sqr) {
                        double force_div_rlen = _INTERPOLATE_FORCE(sqrt(rsqr));
                        vec_t force;
                        Vec3_MUL(force, r, force_div_rlen);
                        Vec3_SUBTO(sum_force, force);
                }
        }
        return sum_force;
}

#undef _PER_SEP
#undef _INTERPOLATE_FORCE

/* * * * * * * *
 * Integration *
 *---------------------------------------------
 * performs one CEX_dt time step integration.  use subcyles where necessary (see below).
 * returns 1 if neighbor lists are now invalid as a result of this integration
 * cycle.  0 otherwise 
 */

/* paramters describing an integration of time dt */
typedef struct {
        double dt_inv_gamma;
        double B2;
} subcycle_parameters;

static inline subcycle_parameters gen_subcycle_parameters(int divisions)
        GCC_ATTRIBUTE((always_inline));

static vec_t integrate_brownian_subcycle(int i_particle, double dU, vec_t rnd);

static int random_numbers_fresh = 0;

static inline void
update_random(void)
{
        if (!random_numbers_fresh) {
                subcycle_parameters sp0 = gen_subcycle_parameters(1);
                CEX_generate_gauss(CEX_random_vectors, sp0.B2);
                random_numbers_fresh = 1;
        }
}

/* performs one CEX_dt time step integration.  use subcyles where necessary (see below).
 * returns 1 if neighbor lists are now invalid as a result of this integration
 * cycle.  0 otherwise 
 */

/* rewrite this to be vectorized, including generation of random numbers 
 * in one large, specifically use intel math kernal library (MKL)
 */
static int
integrate_cycle(void)
{
        assert(ARR_LENGTH(CEX_positions) >= CEX_N_internal_particles);
        assert(ARR_LENGTH(CEX_forces) == CEX_N_internal_particles);
        assert(ARR_LENGTH(CEX_nl_displace) == CEX_N_internal_particles);
        assert(ARR_LENGTH(CEX_new_positions) == CEX_N_internal_particles);
        assert(ARR_LENGTH(CEX_random_vectors) == CEX_N_internal_particles);

        vec_t * CEX_RESTRICT _positions = ARR_DATA_AS(vec_t, CEX_positions);
        vec_t * CEX_RESTRICT _new_positions = ARR_DATA_AS(vec_t, CEX_new_positions);
        vec_t * CEX_RESTRICT _forces = ARR_DATA_AS(vec_t, CEX_forces);
        vec_t * CEX_RESTRICT _nl_displace = ARR_DATA_AS(vec_t, CEX_nl_displace);
        vec_t * CEX_RESTRICT _rnd_force = ARR_DATA_AS(vec_t, CEX_random_vectors);
        vec_t _box_size = CEX_box_size;
        subcycle_parameters sp0 = gen_subcycle_parameters(1);
        int displace_beyond_nl = 0;
        double _dU_max = CEX_dU_max;

        update_random();

#ifdef USE_OMP_INTEGRATE
#  pragma omp parallel for schedule(static) firstprivate(_positions, _new_positions, _forces, _nl_displace, _rnd_force, \
                                                        _box_size, sp0, displace_beyond_nl, _dU_max)
#endif
        for (int i_particle=CEX_N_internal_particles-1; i_particle>=0; i_particle--) {
                // dx = dt/gamma * F(x,t) + sqrt(2kT*dt/gamma)*R_gauss
                // dx = s.dt_inv_gamma * F(x,t) + sp.B2*R_gauss
                vec_t force = _forces[i_particle];

                vec_t rforce = _rnd_force[i_particle];
                vec_t delta;
                Vec3_MUL(delta, force, sp0.dt_inv_gamma);
                Vec3_ADDTO(delta, rforce);
                double dU = fabs(Vec3_DOT(delta, force));
                if (unlikely(dU > _dU_max)) {
                        vec_t r;
                        Vec3_DIV(r, rforce, sp0.B2);
                        delta = integrate_brownian_subcycle(i_particle, dU, r);
                } else {
                        vec_t position = _positions[i_particle];
                        Vec3_ADDTO(position, delta);
                        XPERIODIZE_LOCATION(position, _box_size);
                        if (unlikely(position.x > _box_size.x ||
                                     position.y > _box_size.y ||
                                     position.z > _box_size.z)) {
                                printf("Bad position " Vec3_FRMT("%.3f") " (nm)",
                                       Vec3_ARGS_SCALED(1e9, position));
                                       
                                abort();
                        }
                        _new_positions[i_particle] = position;
                }
                vec_t nl_displace = _nl_displace[i_particle];
                Vec3_ADDTO(nl_displace, delta);
                _nl_displace[i_particle] = nl_displace;
                displace_beyond_nl |= Vec3_SQR(nl_displace) > r_delta_2_sqr;
        }
        XMEMCPY(vec_t, _positions, _new_positions, CEX_N_internal_particles);
        random_numbers_fresh = 0;
        return displace_beyond_nl;
}

static inline subcycle_parameters
gen_subcycle_parameters(int divisions)
{
        subcycle_parameters sp;
        assert(divisions>=1);
        double sdt = CEX_dt / divisions;
        sp.dt_inv_gamma = sdt / CEX_fric_gamma;
        sp.B2 = sqrt(2*CEX_kT*sdt / CEX_fric_gamma);
        return sp;
}

/* when a particle is in an unusually high gradient, we perform
 * the integration as a series of sub-integration.  this allows
 * us to use a moderatley large integration time step, without
 * incurring a leaky hamiltonian.  this can introduce problems 
 * with a stochastic force field, in that sub-integration becomes
 * correlated with large fluctuations.  to properly handle this 
 * correlation, we never toss out the stochastic forces that caused
 * an individual cycle to go over our tolerance.  instead we record
 * the thermal force trajectory and simply apply them over smaller 
 * and smaller time steps, until each individual integration is within
 * our tolerance.  real fun bookkeeping!
 */

typedef struct rl_cell rl_cell;
struct rl_cell {
        vec_t rnd;
        rl_cell * next;
};
static rl_cell * rl_free_cells;
static rl_cell * rl_head=NULL;

#ifdef USE_OMP_INTEGRATE
# pragma opm threadprivate(rl_free_cells, rl_head)
#endif

#define RL_BLOCK_SIZE 16
static void
rl_grow(void)
{
        int i;
        rl_cell *block, *head, *ptr;

        assert(rl_free_cells==NULL);
        block = XNEW(rl_cell, RL_BLOCK_SIZE);
        XBZERO(rl_cell, block, RL_BLOCK_SIZE);
        head = NULL;
        for (i=RL_BLOCK_SIZE, ptr=block; i--;) {
                ptr->next = head;
                head = ptr++;
        }
        rl_free_cells = head;
}

static void
rl_clear(void)
{
        rl_cell *free, *head, *cnext;
        
        free = rl_free_cells;
        head = rl_head;
        rl_head = NULL;
        while (head) {
                cnext = head->next;
                head->next = free;
                free = head;
                head = cnext;
        }
        rl_free_cells = free;
}

static inline void
rl_push(vec_t rnd)
{
        rl_cell *cell;

        if (unlikely(rl_free_cells==NULL)) {
                rl_grow();
        }
        cell = rl_free_cells;
        rl_free_cells = cell->next;
        cell->rnd = rnd;
        cell->next = rl_head;
        rl_head = cell;
}

static inline vec_t
rl_pop_or_create(void)
{
        if (unlikely(rl_head)) {
                rl_cell *top;
                top = rl_head;
                rl_head = top->next;
                top->next = rl_free_cells;
                rl_free_cells = top;
                return top->rnd;
        } else {
                vec_t v;
                CEX_generate_gauss_vector(&v, 1.0);
                return v;
        }
}

typedef struct {
        vec_t delta;
        double dU_max;
} subcycle_result;

static subcycle_result 
do_integrate_subcycle(subcycle_parameters sp, int i_particle, 
                      vec_t position, int n_subcycles);

static vec_t
integrate_brownian_subcycle(int i_particle, double dU, vec_t rnd)
{
        int n_subcycles;
        subcycle_parameters sp;
        subcycle_result results;
        vec_t position;

        position = ARR_INDEX_AS(vec_t, CEX_positions, i_particle);
        n_subcycles = 1 + (int)ceil(CEX_dU_max/dU);
        rl_push(rnd);
        for (;;) {
                sp = gen_subcycle_parameters(n_subcycles);
                results = do_integrate_subcycle(sp, i_particle, position, n_subcycles);
                if (results.dU_max < CEX_dU_max)
                        break;
                n_subcycles <<= 1;
                assert(n_subcycles > 1);
                if (unlikely(n_subcycles > 4)) {
                        //xprintf("subcycle integration with n_subcycles=%d", n_subcycles);
                }
        }
        rl_clear();
        Vec3_ADDTO(position, results.delta);
        PERIODIZE_LOCATION(position);
        ARR_INDEX_AS(vec_t, CEX_new_positions, i_particle) = position;
        return results.delta;
}

static subcycle_result
do_integrate_subcycle(subcycle_parameters sp, int i_particle, 
                     vec_t position, int n_subcycles)
{
        subcycle_result res;
        if (n_subcycles==0) {
                Vec3_CLEAR(res.delta);
                res.dU_max = 0.0;
        } else {
                vec_t force = eval_one_force(i_particle, position);
                vec_t r = rl_pop_or_create();
                vec_t delta, rforce;
                Vec3_MUL(delta, force, sp.dt_inv_gamma);
                Vec3_MUL(rforce, r, sp.B2);
                Vec3_ADDTO(delta, rforce);
                double dU = fabs(Vec3_DOT(delta, rforce));
                Vec3_ADDTO(position, delta);
                PERIODIZE_LOCATION(position);
                res = do_integrate_subcycle(sp, i_particle, position, n_subcycles-1);
                if (dU > res.dU_max) {
                        res.dU_max = dU;
                }
                Vec3_ADDTO(res.delta, delta);
                rl_push(r);
        }
        return res;
}

/* * * * * * * * * *
 * simulation loop *
 * * * * * * * * * */
#define CMD_UPDATE_NEIGHBORS 1
#define CMD_UPDATE_FORCES 2
#define CMD_INTEGRATE_ONE 3
#define CMD_EXIT_LOOP 255

void
CEX_slave_simulation_loop(void)
{
        int exit_loop, cmd, ret, recv;

        REQ_SLAVE();
        REQ_INIT();
        for (exit_loop=0; !exit_loop;) {
                MPI_Bcast(&cmd, 1, MPI_INT, 0, MPI_COMM_WORLD);
                switch (cmd) {
                case CMD_UPDATE_NEIGHBORS:
                        CEX_thread_update_neighbors();
                        ret = 0;
                        break;
                case CMD_UPDATE_FORCES:
                        update_forces();
                        ret = 0;
                        break;
                case CMD_INTEGRATE_ONE:
                        ret = integrate_cycle();
                        break;
                case CMD_EXIT_LOOP:
                        exit_loop = 1;
                        ret = 0;
                        break;
                default:
                        Fatal("unkown command %d", cmd);
                }
                MPI_Reduce(&ret, &recv, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        }
}

static inline void
tell_slaves(int cmd)
{
        if (HAVE_JUNCTIONS()) {
                MPI_Bcast(&cmd, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
}

static inline int
poll_slaves(void)
{
        int send, res;

        res = send = 0;
        if (HAVE_JUNCTIONS()) {
                MPI_Reduce(&send, &res, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        return res;
}

static inline void
update_neighbors_everywhere(void)
{
        tell_slaves(CMD_UPDATE_NEIGHBORS);
        CEX_thread_update_neighbors();
        poll_slaves();
}

static inline void 
update_forces_everywhere(void)
{
        tell_slaves(CMD_UPDATE_FORCES);
        update_forces();
        poll_slaves();
}

static inline int
integrate_every_where(void)
{
        int displace_beyond_nl;

        tell_slaves(CMD_INTEGRATE_ONE);
        displace_beyond_nl = integrate_cycle();
        return displace_beyond_nl | poll_slaves();
}

static inline void
exit_loop_everywhere(void)
{
        tell_slaves(CMD_EXIT_LOOP);
        poll_slaves();
}

void
CEX_master_simulate_cycles(int cycles)
{
        int integrate_cycles;
        int displace_beyond_nl;

        REQ_MASTER();
        REQ_INIT();
        if (cycles<0) {
                Fatal("bad number of cycles %d", cycles);
        }
        update_neighbors_everywhere();
        while (likely(cycles)) {
                update_forces_everywhere();
                integrate_cycles = CEX_force_update_rate;
                while (likely(integrate_cycles--)) {
                        if ((cycles%50)==0) {
                                //xprintf("%d cyles", cycles);
                        }
                        displace_beyond_nl = integrate_every_where();
                        cycles --;
                        if (displace_beyond_nl) {
                                integrate_cycles = 0;
                                update_neighbors_everywhere();
                        }
                }
        }
        exit_loop_everywhere();
}

