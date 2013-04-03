/* -*- Mode: c -*-
 * bd.h - Brownian Dynamics Integrator
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

#ifndef _BD_H
#define _BD_H

#include "array.h"
#include "msg.h"

/* Global Information, Non-Thread Specific
 *--------------------------------------------*/
extern double CEX_T; /* absoulte temperature (kelvin) */
extern double CEX_kT; /* kB * CEX_T */
extern double CEX_dt;/* duration of 1 numerical integration (seconds)*/
extern double CEX_fric_gamma; /* solvent/particle friction (pascal*seconds*meters)*/
/* maximum change in energy for a single particle integration (joules)
 * larger changes are done as a series of sub-integrations */
extern double CEX_dU_max;

/* Force Field
 *--------------*/
/* linear interpolation table */
typedef struct {
        double x_min, x_prec;
        array_t *table;
} linterp_table_t;

/* cutoff distance of pair potential (and thereby force) (meters) */
extern double CEX_r_pair_cutoff;

/* pair potential (currently not used)
 * x_min,x_prec (meters); table (joules) */
extern linterp_table_t CEX_pair_potential ;

/* pair force
 * in this table, the force is scaled by the distance s.t. f(|r|)*r
 * is the force vector between two particles when r is the separation
 * vector.  
 * x_min,x_prec (meters); table (newton/meters) */
extern linterp_table_t CEX_pair_force ;

/* revaluate pair forces every `x cycles */
extern int CEX_force_update_rate;

/* maximum distance between 2 particles to consider them neighbors (meters) 
 * this controls how long the neighbor lists are valid, larger values
 * allows them to be valid for long times while resulting in more neighbors pairs
 * and more unnessary pairs to check durring force evaluation*/
extern double CEX_r_neighbor;
/* square of above (meters^2) */
extern double CEX_r_neighbor_sqr;

/* State for This Thread
 *-----------------------------------------------------------------*/
/* positions; includes both particles that belong to this
 * cell as well as externally updated neighbor positions. (meters)
 * First CEX_N_internal_particles are internal, and rest are
 * exteranl*/
extern array_t *CEX_positions;
extern array_t *CEX_new_positions;

/* intenral (particle-wise data structures) */
extern int CEX_N_internal_particles;
extern array_t *CEX_tags; /* unique tags for tracing particle trajectories */
extern array_t *CEX_forces;
extern array_t *CEX_random_vectors;
/* forces on internal particles (newtons) */
/* absolute displacement of internal particles since
 * neighbor lists were last rebuilt (meters).  we use this to determine
 * when the neighbor lists could potentially be invalid and 
 * need rebuilt*/
extern array_t *CEX_nl_displace;
/* neighbor list of important neighbor pairs for forces, 
 * stored as pairs of indices in CEX_positions s.t. that the
 * indices of the two particles in this i-th neighbor are
 * i*2 and i*2 + 1 */
extern array_t *CEX_internal_neighbors;
/* first indice is that of the internal particle,
 * second is that of the external particle */
extern array_t *CEX_external_neighbors;

/* auxillay data sturctures for each communicator */
extern array_t *CEX_send_indices,
               *CEX_recv_lengths,
               *CEX_tmp_recv_positions,
               *CEX_tmp_recv_tags,
               *CEX_ext_positions_offset,
               *CEX_remove_indices,
               *CEX_send_positions_buffers;

void CEX_thread_update_neighbors(void);
void CEX_thread_update_forces(void);
void CEX_slave_simulation_loop(void);
void CEX_master_simulate_cycles(int cycles);


#endif /* _BD_H */

