/* -*- Mode: c -*-
 * random.h - Generic random number interface
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

#ifndef _RANDOM_H
#define _RANDOM_H

#include "array.h"

void CEX_seed_random(unsigned int seed)
        GCC_ATTRIBUTE((noinline));
void CEX_generate_gauss(array_t *place, double sigma2)
        GCC_ATTRIBUTE((noinline));
void CEX_generate_gauss_vector(vec_t *vec, double sigma2)
        GCC_ATTRIBUTE((noinline));

#endif /* _RANDOM_H */
