/* -*- Mode: c -*-
 * comm.h - Pair-wise Communications Utilities
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


/* The majority of communication that occurs in this simulation is 
 * pair-wise communcation between spacially ajacent threads.  
 * In such communication, each thread must send and recieve arbitrary
 * data with its adjacent threads.  To prevent dead lock, we 
 * have to perform this communication in synchronized order.  This is
 * done with comm_rule_t's which define when each thread sends and
 * recieves with its neighboring threads.  These are computed for specific
 * cell configurations in `higher` level code. 
 */

#ifndef _COMM_H
#define _COMM_H

#include <mpi.h>
#include "array.h"
#include "opt.h"

#define COMM_INST_SEND 1
#define COMM_INST_RECV 2

extern int CEX_rank;
#define IS_MASTER() (CEX_rank==0)
#define IS_SLAVE() (CEX_rank!=0)
extern int CEX_size;

#define REQ_MASTER() do {                    \
      if (unlikely(!IS_MASTER())) {          \
              Fatal("thread %d executed a master only segment as slave", \
                    CEX_rank);                                          \
      }                                                                 \
} while (0)                

#define REQ_SLAVE() do {                    \
      if (unlikely(!IS_SLAVE())) {          \
              Fatal("thread %d executed a slave only segment as master", \
                    CEX_rank);                                          \
      }                                                                 \
} while (0)                


typedef struct comm_t comm_t;
typedef struct comm_rule_t comm_rule_t;

struct comm_rule_t {
        int inst; /* either COMM_INST_SEND or COMM_INST_RECV */
        comm_t *comm; /* thread we're sending or recieving to */
        int tag; /* unique tag describing message */
};

/* comm rules for this thread */
extern array_t * CEX_comm_rules;

#define COMM_RULE_FOREACH(ptr, counter)                         \
        ARR_FOREACH(comm_rule_t, CEX_comm_rules, ptr, counter)

struct comm_t {
        int comm_rank; /*rank of the thread we're communicating with*/
        int arr_inx; /* auxillary data for communicators is stored in arrays.
                      * each communicator has a unique index to access 
                      * this meta data */
        comm_rule_t *current_rule; /* rule we're currently executing */
};

extern array_t * CEX_comms;

#define COMM_FOREACH(ptr, counter)                         \
        ARR_FOREACH(comm_t, CEX_comms, ptr, counter)

#define GET_COMM_DATA(tp, arr, comm) ARR_INDEX_AS(tp, arr, (comm)->arr_inx)


/* * * * * * * * * * * * * * * * * * * * *
 * inlined communicator helper functions *
 * * * * * * * * * * * * * * * * * * * * */

/* communcation routines  */
#define COMM_HELPER_ATTRS              \
        GCC_ATTRIBUTE((always_inline))

static inline void comm_send(comm_t *, void *, int, MPI_Datatype) COMM_HELPER_ATTRS;
static inline void comm_recv(comm_t *, void *, int, MPI_Datatype) COMM_HELPER_ATTRS;
static inline void comm_send_int(comm_t *, int) COMM_HELPER_ATTRS;
static inline int comm_recv_int(comm_t *) COMM_HELPER_ATTRS;
static inline int msg_bytes(comm_t *) COMM_HELPER_ATTRS;
#define MSG_SIZE(comm, dtype) (msg_bytes(comm) / sizeof(dtype))
static inline void comm_send_ints(comm_t *, array_t *) COMM_HELPER_ATTRS;
static inline void comm_send_vecs(comm_t *, array_t *) COMM_HELPER_ATTRS;
static inline void comm_send_ints_by_index(comm_t *, array_t *ints, 
                                           array_t *indices) COMM_HELPER_ATTRS;
static inline void comm_send_vecs_by_index(comm_t *, array_t *vecs, 
                                           array_t *indices) COMM_HELPER_ATTRS;
static inline void comm_recv_ints(comm_t *, array_t *dst) COMM_HELPER_ATTRS;
static inline void comm_recv_vecs(comm_t *, array_t *dst) COMM_HELPER_ATTRS;
static inline void comm_recv_extend_ints(comm_t *, array_t *dst) COMM_HELPER_ATTRS;
static inline void comm_recv_extend_vecs(comm_t *, array_t *dst) COMM_HELPER_ATTRS;

#include "comm-inline.h"

/* this hideous construct organizes pairwise communcation stages
 * by executing either SEND_BODY or RECV_BODY at the proper
 * times and binding PC_VAR to the communicator for this stage
 */

#define DO_COMM(PC_VAR, SEND_BODY, RECV_BODY) do {          \
        int _cr_c;                                          \
        comm_rule_t *_cr_p;                                 \
        comm_t *PC_VAR;                                     \
        MPI_Barrier(MPI_COMM_WORLD);                        \
        COMM_RULE_FOREACH(_cr_p, _cr_c) {                   \
            PC_VAR = _cr_p->comm;                           \
            PC_VAR->current_rule = _cr_p;                   \
            switch (_cr_p->inst) {                          \
            case COMM_INST_SEND:                            \
                    { SEND_BODY; }                          \
                    break;                                  \
            case COMM_INST_RECV:                            \
                    { RECV_BODY; }                          \
                    break;                                  \
            default:                                        \
                   Fatal("unkown comm instruction %d",      \
                         _cr_p->inst);                      \
            }                                               \
        }                                                   \
} while (0)

#endif /* _COMM_H */
