/* -*- Mode: c -*-
 * init.c - Routines to initialize simulation from control process
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

#include <math.h>
#include <float.h>
#include <omp.h>

#include "constants.h"
#include "debug.h"
#include "mem.h"
#include "random.h"
#include "periodic.h"
#include "comm.h"
#include "cells.h"
#include "bd.h"
#include "msg.h"
#include "init.h"

static const char *init_state = "uninitialized";

#define REQ_STATE(state) do {                   \
        if (strcmp(state, init_state)!=0) {     \
                Fatal("require state " state \
                      "but currently in state %.50s", init_state);      \
        }                                                               \
} while(0)


/* Helper Functions */
static void
check_name(msg_t *msg, const char *name)
{
        //xprintf("reading %.50s", name);
        array_t *char_arr = CEX_msg_read_char_array(msg);
        char buffer[ARR_LENGTH(char_arr)+1];
        CEX_char_array_as_string(char_arr, buffer);
        if (strcmp(buffer, name)!=0) {
                Fatal("alignment issue; expected field named %.50s; read %.50s",
                      name, buffer);
        }
}

#define VALUE_CHECK(frmt) ({ do{                                        \
        if (unlikely( value < mn || value > mx)) {                      \
                if (index!=-1) {                                        \
                   Fatal("bad value " frmt " for parameter %.50s[%d]; " \
                         "must be in the range [" frmt ":" frmt "]",    \
                         value, name, index, mn, mx);                   \
           } else {                                                     \
                   Fatal("bad value " frmt "for parameter %.50s;"       \
                         "must be in range [" frmt ":" frmt "]",        \
                         value, name, mn, mx);                          \
           }                                                            \
        }                                                               \
} while (0); value;})
  
static inline int
check_int(const char *name, int index, int value, int mn, int mx)
{
        return VALUE_CHECK("%d");
}

static inline double
check_double(const char *name, int index, double value, double mn, double mx)
{
        return VALUE_CHECK("%.6g");
}

static inline vec_t
check_vec(const char *name, int index, vec_t value, vec_t mn, vec_t mx)
{
        int len=strlen(name);
        char buffer[len+3];
        memcpy(buffer, name, len);
        buffer[len+2] = '\0';
        buffer[len] = '.';
        buffer[len+1] = 'x';
        check_double(buffer, index, value.x, mn.x, mx.x);
        buffer[len+1] = 'y';
        check_double(buffer, index, value.y, mn.y, mx.y);                
        buffer[len+1] = 'z';
        check_double(buffer, index, value.z, mn.z, mx.z);
        return value;
}

static inline int
read_int(msg_t *msg, const char *name, int mn, int mx)
{
        check_name(msg, name);
        return check_int(name, -1, CEX_msg_read_int(msg), mn, mx);
}

static inline double
read_double(msg_t *msg, const char * name, double mn, double mx)
{
        check_name(msg, name);
        return check_double(name, -1, CEX_msg_read_double(msg), mn, mx);
}

static inline array_t *
read_int_array(msg_t *msg, const char *name, int mn, int mx)
{
        check_name(msg, name);
        array_t *arr = CEX_msg_read_int_array(msg);
        for (int i=0; i<ARR_LENGTH(arr); i++) {
                check_int(name, i, ARR_INDEX_AS(int, arr, i), mn, mx);
        }
        return arr;
}

static inline array_t *
read_double_array(msg_t *msg, const char *name, double mn, double mx)
{
        check_name(msg, name);
        array_t *arr = CEX_msg_read_double_array(msg);
        for (int i=0; i<ARR_LENGTH(arr); i++) {
                check_double(name, i, ARR_INDEX_AS(double, arr, i), mn, mx);
        }
        return arr;
}

static inline vec_t
read_vec(msg_t *msg, const char *name, vec_t mn, vec_t mx)
{
        check_name(msg, name);
        return check_vec(name, -1, CEX_msg_read_vec(msg), mn, mx);
}

static inline vec_t
read_vec_sclims(msg_t *msg, const char *name, double smn, double smx)
{
        vec_t vmn, vmx;
        Vec3_SETALL(vmn, smn);
        Vec3_SETALL(vmx, smx);
        return read_vec(msg, name, vmn, vmx);
}

static inline array_t *
read_vec_array(msg_t *msg, const char *name, vec_t mn, vec_t mx)
{
        check_name(msg, name);
        array_t *arr = CEX_msg_read_vec_array(msg);
        for (int i=0; i<ARR_LENGTH(arr); i++) {
                check_vec(name, i, ARR_INDEX_AS(vec_t, arr, i), mn, mx);
        }
        CEX_align_array(arr, sizeof(double));
        return arr;
}

static void
read_linterp(msg_t *msg, const char *name, linterp_table_t *tbl)
{
        check_name(msg, name);
        tbl->x_min = read_double(msg, "x_min", 0, DBL_MAX);
        tbl->x_prec = read_double(msg, "x_prec", 0, DBL_MAX);
        tbl->table = read_double_array(msg, "table", -DBL_MAX, DBL_MAX);
        CEX_align_array(tbl->table, sizeof(double));
}

void
CEX_initialize_system(msg_t *msg)
{
        REQ_STATE("uninitialized");
        CEX_box_size = read_vec_sclims(msg, "box_size", 2*CEX_R_particle, 1e5*CEX_R_particle);
        Vec3_MUL(CEX_box_half, CEX_box_size, 0.5);
        CEX_T = read_double(msg, "T", 1, 2000);
        CEX_kT = CEX_kB * CEX_T;
        CEX_dt = read_double(msg, "dt", CEX_ps, 1000*CEX_ns);
        CEX_dU_max = read_double(msg, "dU_max", 1e-3*CEX_kT, 200*CEX_kT);
        CEX_fric_gamma = read_double(msg, "fric_gamma", 1e-14, 1e-8);
        CEX_force_update_rate = read_int(msg, "force_update", 1, 1000);
        CEX_r_pair_cutoff = read_double(msg, "r_pair_cutoff", 2*CEX_R_particle, 5*CEX_R_particle);
        read_linterp(msg, "pair_potential", &CEX_pair_potential);
        read_linterp(msg, "pair_force", &CEX_pair_force);
        CEX_r_neighbor = read_double(msg, "r_neighbor", 2*CEX_R_particle, 10*CEX_R_particle);
        CEX_r_neighbor_sqr = CEX_r_neighbor * CEX_r_neighbor;
        xprintf("system init: box_size " Vec3_FRMT("%.1f")  "(nm) "
                "T=%.2F (K) "
                "dt=%.1f (ps) "
                "dU_max=%.1f (kT) "
                "fric_gamma=%.2f (pN*ns/nm) ",
                Vec3_ARGS_SCALED((1.0/CEX_nm), CEX_box_size),
                CEX_T, 
                CEX_dt/CEX_ps,
                CEX_dU_max / CEX_kT,
                CEX_fric_gamma / CEX_pN / CEX_ns * CEX_nm);
        xprintf("r_pair_cutoff=%.2f (nm) "
                "r_pair_neighor=%.2f (nm) ",
                CEX_r_pair_cutoff / CEX_nm,
                CEX_r_neighbor / CEX_nm);              
        init_state = "system";
}

void 
CEX_initialize_random(msg_t *msg)
{
        REQ_STATE("system");
        unsigned int seed = CEX_msg_read_uint(msg);
        CEX_seed_random(seed);
        vec_t pull;
        CEX_generate_gauss_vector(&pull, 1.0);
        xprintf("random seeded with 0x%X; first gaussian vector " Vec3_FRMT("%.3f"),
                seed, Vec3_ARGS(pull));
        init_state = "random";
}

static vec_t
read_extent(msg_t *msg, const char *name)
{
       vec_t vzero;
       Vec3_SETALL(vzero, 0);
       return read_vec(msg, name, vzero, CEX_box_size);
}

void
CEX_initialize_cell_state(msg_t *msg)
{
        REQ_STATE("random");
        CEX_this_cell = XNEW(cell_t, 1);
        CEX_this_cell->comm = NULL;
        CEX_this_cell->min_extent = read_extent(msg, "min_extent");
        CEX_this_cell->max_extent = read_extent(msg, "max_extent");
        array_t *positions = read_vec_array(msg, "positions",  
                                            CEX_this_cell->min_extent, 
                                            CEX_this_cell->max_extent);
        array_t *tags = read_int_array(msg, "tags", 0, 10000000);
        if (ARR_LENGTH(positions) != ARR_LENGTH(tags)) {
                Fatal("inconsistent positions and tags length: %lu and %lu respectively",
                      ARR_LENGTH(positions), ARR_LENGTH(tags));
        }
        CEX_N_internal_particles = ARR_LENGTH(positions);
        CEX_positions = positions;
        CEX_new_positions = CEX_make_vec_array(CEX_N_internal_particles);
        CEX_tags = tags;
        CEX_align_array(CEX_tags, sizeof(int));
        CEX_forces = CEX_make_vec_array(CEX_N_internal_particles);
        CEX_align_array(CEX_forces, sizeof(double));
        CEX_random_vectors = CEX_make_vec_array(CEX_N_internal_particles);
        CEX_align_array(CEX_random_vectors, sizeof(double));
        CEX_nl_displace = CEX_make_vec_array(CEX_N_internal_particles);
        CEX_align_array(CEX_nl_displace, sizeof(double));
        ARR_LENGTH(CEX_forces) = CEX_N_internal_particles;
        ARR_LENGTH(CEX_random_vectors) = CEX_N_internal_particles;
        ARR_LENGTH(CEX_nl_displace) = CEX_N_internal_particles;
        CEX_internal_neighbors = CEX_make_int_array(2*CEX_N_internal_particles);
        CEX_align_array(CEX_internal_neighbors, sizeof(int));
        CEX_external_neighbors = CEX_make_int_array(0);
        CEX_align_array(CEX_external_neighbors, sizeof(int));
        xprintf("initialized cell with extent " 
                Vec3_FRMT("%.2f") " -> " Vec3_FRMT("%.2f") "(nm) "
                "%d internal particles",
                Vec3_ARGS_SCALED((1.0/CEX_nm), CEX_this_cell->min_extent),
                Vec3_ARGS_SCALED((1.0/CEX_nm), CEX_this_cell->max_extent),
                CEX_N_internal_particles);
        init_state = "cell-state";
}

static void
read_comm_el(msg_t *msg, array_t *arr)
{
        comm_t comm;
        comm.comm_rank = read_int(msg, "comm_rank", 0, CEX_size-1);
        comm.arr_inx = ARR_LENGTH(arr);
        comm.current_rule = NULL;
        ARR_APPEND(comm_t, arr, comm);
}

static comm_t*
read_comm_ptr(msg_t *msg)
{
       int comm_index = read_int(msg, "comm_index", 0, ARR_LENGTH(CEX_comms)-1);
       return (comm_t *)ARR_ADDRESS_ELEMENT(CEX_comms, comm_index);
       //return &ARR_INDEX_AS(comm_t, CEX_comms, comm_index);
}

static void
read_comm_rule_el(msg_t *msg, array_t *arr)
{
        comm_rule_t rule;
        rule.inst = read_int(msg, "inst", 1, 2);
        rule.comm = read_comm_ptr(msg);
        rule.tag = read_int(msg, "tag", 0, 1000000);
        ARR_APPEND(comm_rule_t, arr, rule);
}

static array_t *
read_array(msg_t *msg, const char *name, size_t el_size, msg_element_reader reader)
{
        check_name(msg, name);
        return CEX_msg_read_array(msg, el_size, reader);
}

static array_t *
make_array_of_arrayps(size_t el_size, int length)
{
        array_t *arr = CEX_make_array(sizeof(void*), length);
        for (int i=length; i-->0;) {
                array_t *el = CEX_make_array(el_size, 0);
                CEX_align_array(el, el_size);
                ARR_APPEND(array_t *, arr, el);
        }
        return arr;
}

void
CEX_initialize_cell_comm(msg_t *msg)
{
        REQ_STATE("cell-state");
        CEX_comms = read_array(msg, "comms", sizeof(comm_t), &read_comm_el);
        CEX_comm_rules = read_array(msg, "comm_rules", sizeof(comm_rule_t), 
                                    &read_comm_rule_el);
        xprintf("initialized %lu communicators and %lu communication rules",
                ARR_LENGTH(CEX_comms), ARR_LENGTH(CEX_comm_rules));
        /* setup auxillary data structures */
        int N_comms = ARR_LENGTH(CEX_comms);
        CEX_send_indices = make_array_of_arrayps(sizeof(int), N_comms);
        CEX_recv_lengths = CEX_make_int_array(N_comms);
        CEX_tmp_recv_positions = make_array_of_arrayps(sizeof(vec_t), N_comms);
        CEX_tmp_recv_tags = make_array_of_arrayps(sizeof(int), N_comms);
        CEX_ext_positions_offset = make_array_of_arrayps(sizeof(int), N_comms);
        CEX_remove_indices = make_array_of_arrayps(sizeof(int), N_comms);
        CEX_send_positions_buffers = make_array_of_arrayps(sizeof(vec_t), N_comms);
        init_state = "cell-comm";
}

static void
read_jcell_el(msg_t *msg, array_t *arr)
{
        cell_t jcell;
        jcell.comm = read_comm_ptr(msg);
        jcell.min_extent = read_extent(msg, "min_extent");
        jcell.max_extent = read_extent(msg, "max_extent");
        ARR_APPEND(cell_t, arr, jcell);
}

static cell_t *
read_cell(msg_t *msg)
{
       int cell_index = read_int(msg, "cell_index", 0, ARR_LENGTH(CEX_jcells)-1);
       //xprintf("cell_index %d %p", cell_index, ARR_ADDRESS_ELEMENT(CEX_jcells, cell_index));
       return (cell_t *)ARR_ADDRESS_ELEMENT(CEX_jcells, cell_index);
}

static int
read_axis(msg_t *msg)
{
        return read_int(msg, "axis", 0, 2);
}

static void
read_surface_junction_el(msg_t *msg, array_t *arr)
{
        surface_junction_t sj;
        sj.cell = read_cell(msg);
        sj.axis = read_axis(msg);
        sj.dir = read_int(msg, "dir", -1, 1);
        if (sj.dir==0) {
                Fatal("bad direction 0");
        }
        ARR_APPEND(surface_junction_t, arr, sj);
}

static void
read_line_junction_el(msg_t *msg, array_t *arr)
{
        line_junction_t lj;
        lj.cell = read_cell(msg);
        lj.axis = read_axis(msg);
        int axis1, axis2;
        switch(lj.axis) {
        case AXIS_X:
                axis1 = AXIS_Y;
                axis2 = AXIS_Z;
                break;
        case AXIS_Y:
                axis1 = AXIS_X;
                axis2 = AXIS_Z;
                break;
        case AXIS_Z:
                axis1 = AXIS_X;
                axis2 = AXIS_Y;
                break;
        default:
                Fatal("bad axis %d", lj.axis);
        }
        lj.offset1 = read_double(msg, "offset1", 
                                 INDEX_AXIS(&CEX_this_cell->min_extent, axis1),
                                 INDEX_AXIS(&CEX_this_cell->max_extent, axis1));
        lj.offset2 = read_double(msg, "offset2", 
                                 INDEX_AXIS(&CEX_this_cell->min_extent, axis2),
                                 INDEX_AXIS(&CEX_this_cell->max_extent, axis2));
        ARR_APPEND(line_junction_t, arr, lj);
}

static void
read_point_junction_el(msg_t *msg, array_t *arr)
{
        point_junction_t pj;
        pj.cell = read_cell(msg);
        pj.offset = read_vec(msg, "offset", 
                             CEX_this_cell->min_extent,
                             CEX_this_cell->max_extent);
        ARR_APPEND(point_junction_t, arr, pj);
}


void
CEX_initialize_cell_junctions(msg_t *msg)
{
        REQ_STATE("cell-comm");
        CEX_jcells = read_array(msg, "jcells",
                                sizeof(cell_t), &read_jcell_el);
        CEX_surface_junctions = read_array(msg, "surface_junctions", 
                                           sizeof(surface_junction_t),
                                           &read_surface_junction_el);
        CEX_line_junctions = read_array(msg, "line_junctions", 
                                           sizeof(line_junction_t),
                                           &read_line_junction_el);
        CEX_point_junctions = read_array(msg, "point_junctions",
                                         sizeof(point_junction_t),
                                         &read_point_junction_el);
        xprintf("initialized %lu junctions: surface=%lu line=%lu point=%lu",
                ARR_LENGTH(CEX_jcells),
                ARR_LENGTH(CEX_surface_junctions), 
                ARR_LENGTH(CEX_line_junctions),
                ARR_LENGTH(CEX_point_junctions));
        init_state = "initialized";
}

int
CEX_is_initialized(void)
{
        return strcmp(init_state, "initialized")==0;
}
