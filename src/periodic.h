/* -*- Mode: c -*-
 * xperiodic.h - Periodic cartesian space for this model
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

#ifndef _PERIODIC_H
#define _PERIODIC_H

#include "xperiodic.h"
#include "vector.h"

/* dimensions of carteisan space (meters) */
extern vec_t CEX_box_size;
/* half the dimension (meters). true length scale of periodic cell as
 * not 1-dimensional distance can be larger than this */
extern vec_t CEX_box_half;

#define PERIODIC_SEPARATION_VECTOR(r, pos_i, pos_j)             \
        XPERIODIC_SEPARATION_VECTOR(r, pos_i, pos_j,            \
                                    CEX_box_size, CEX_box_half)

#define PERIODIZE_LOCATION(v)                   \
        XPERIODIZE_LOCATION(v, CEX_box_size) 

#endif /* _PERIODIC_H */
