/* -*- Mode: c -*-
 * init.h - Routines to initialize simulation from control process
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

#ifndef _INIT_H
#define _INIT_H

/* Setup global parameters that are indepdent of cell*/
void CEX_initialize_system(msg_t *msg);

/* Setup data structures unique to this cell */
void CEX_initialize_random(msg_t *msg)
        GCC_ATTRIBUTE((noinline));
/* I havn't the slightest clue why, bt not inlining this function
 * prevents a fatal error (SIGKILL) when using MKL random number
 * generator.  Probaly an icc compiler bug.
 */

void CEX_initialize_cell_state(msg_t *msg);
void CEX_initialize_cell_comm(msg_t *msg);
void CEX_initialize_cell_junctions(msg_t *msg);

int CEX_is_initialized(void);
#define REQ_INIT() do {                                                 \
        if (!CEX_is_initialized()) {                                    \
                Fatal("cannot procede; require simulation initialization"); \
        }                                                               \
} while (0)

#endif /* _INIT_H */
