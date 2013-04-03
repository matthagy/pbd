/* -*- Mode: c -*-
 * random.h - Intel Math Kernel Library Random Number Generator
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

#include <mkl.h>
#include <omp.h>

#include "debug.h"
#include "array.h"
#include "random.h"

/* 59-bit Multiplicative Congruential Generator */
#define BRNG VSL_BRNG_MCG59
#define METHOD VSL_METHOD_DGAUSSIAN_BOXMULLER2

static int initialized=0;
static VSLStreamStatePtr stream;

static inline void
check_vsl_error(int err)
{
        if (unlikely(err!=0)) {
                Fatal("vsl error %d", err);
        }
}

void 
CEX_seed_random(unsigned int seed)
{
        //xprintf("using mkl random number generator");
        if (initialized) {
                check_vsl_error(vslDeleteStream(&stream));
        } else {
                mkl_set_dynamic(0);
                mkl_set_num_threads(1);
        }
        check_vsl_error(vslNewStream(&stream, BRNG, seed));
        initialized = 1;
}

void 
CEX_generate_gauss(array_t *arr, double sigma)
{
        REQ_VARR(arr);
        if (unlikely(!initialized)) {
                Fatal("random number generator not yet initialized");
        }
        check_vsl_error(vdRngGaussian(METHOD, stream, ARR_LENGTH(arr)*3,
                                      ARR_DATA_AS(double, arr), 0.0, sigma));
}

void
CEX_generate_gauss_vector(vec_t *vec, double sigma)
{
        /* don't cast vec_t as double as we want strict aliasing */
        double buffer[3];
        check_vsl_error(vdRngGaussian(METHOD, stream, 3, buffer, 0.0, sigma));
        vec->x = buffer[0];
        vec->y = buffer[1];
        vec->z = buffer[2];
}
