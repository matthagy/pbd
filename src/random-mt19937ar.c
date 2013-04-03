/* -*- Mode: c -*-
 * random.c - Mersene Twister Random Number Generator
 * Derived from mt19937ar.c 
 * Original Copyright Notice Follows
 *--------------------------------------------------------------------------
 * Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 * All rights reserved.                          
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. The names of its contributors may not be used to endorse or promote 
 *      products derived from this software without specific prior written 
 *      permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <math.h>

#include "opt.h"
#include "vector.h"
#include "array.h"
#include "random.h"

typedef unsigned int rnd_int_t;

static array_t * CEX_rstate = NULL;
static rnd_int_t * CEX_RESTRICT CEX_rnext = NULL;
static int CEX_rleft=-1;
static int initialized = 0;

#define rnd_LIMIT 0xffffffffUL
#define rnd_N 624
#define rnd_M 397
#define rnd_MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define rnd_UMASK 0x80000000UL /* most significant w-r bits */
#define rnd_LMASK 0x7fffffffUL /* least significant r bits */
#define rnd_MIXBITS(u,v) ( ((u) & rnd_UMASK) | ((v) & rnd_LMASK) )
#define rnd_TWIST(u,v) ((rnd_MIXBITS(u,v) >> 1) ^ ((v)&1UL ? rnd_MATRIX_A : 0UL))

void 
CEX_seed_random(unsigned int seed)
{
        if (CEX_rstate==NULL) {
                CEX_rstate = CEX_make_array(sizeof(rnd_int_t), rnd_N);
                CEX_align_array(CEX_rstate, sizeof(rnd_int_t));
        }
        rnd_int_t *rstate = ARR_DATA_AS(rnd_int_t, CEX_rstate);
	rstate[0]= seed & 0xffffffffUL;
	for (int j=1; j<rnd_N; j++) {
		rstate[j] = (1812433253UL * (rstate[j-1] ^ (rstate[j-1] >> 30)) + j); 
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array state[].                        */
		/* 2002/01/09 modified by Makoto Matsumoto             */
		rstate[j] &= 0xffffffffUL;  /* for >32 bit machines */
	}
        CEX_rleft = 1;
        CEX_rnext = rstate;
        initialized = 1;
}

static void
CEX_rnext_state(void)
{
        rnd_int_t *rstate = ARR_DATA_AS(rnd_int_t, CEX_rstate);
        rnd_int_t *p=rstate;
        CEX_rleft = rnd_N;
        CEX_rnext = rstate;
        for (int j=rnd_N-rnd_M+1; --j; p++)
                *p = p[rnd_M] ^ rnd_TWIST(p[0], p[1]);
        for (int j=rnd_M; --j; p++)
                *p = p[rnd_M-rnd_N] ^ rnd_TWIST(p[0], p[1]);
        *p = p[rnd_M-rnd_N] ^ rnd_TWIST(p[0], rstate[0]);
}

static inline rnd_int_t rnd_gen_int32(void)
        GCC_ATTRIBUTE((always_inline));

static inline rnd_int_t
rnd_gen_int32(void)
{
        rnd_int_t y;

        if (--CEX_rleft == 0) CEX_rnext_state();
        y = *CEX_rnext++;
        /* Tempering */
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680UL;
        y ^= (y << 15) & 0xefc60000UL;
        y ^= (y >> 18);
        return y;
}

static inline void rnd_gen_gauss2(double *r1, double *r2)
        GCC_ATTRIBUTE((always_inline));

static inline void
rnd_gen_gauss2(double *r1, double *r2)
{
        double x1,x2,w;
        do {
                x1 = (2.0 / (double)rnd_LIMIT) * (double)rnd_gen_int32() - 1.0;
                x2 = (2.0 / (double)rnd_LIMIT) * (double)rnd_gen_int32() - 1.0;
                w = x1*x1 + x2*x2;
        } while (unlikely(w>=1.0));
        w = sqrt(-2.0*log(w)/w);
        *r1 = x1*w;
        *r2 = x2*w;
}

void 
CEX_generate_gauss(array_t *arr, double sigma)
{
        REQ_VARR(arr);
        if (unlikely(!initialized)) {
                Fatal("random number generator not yet initialized");
        }
        double * CEX_RESTRICT place = ARR_DATA_AS(double, arr);
        int length = ARR_LENGTH(arr) * 3;
        int half_length = length >> 1;
        for (int i=0; i<half_length; i++) {
                double *ptr = place + (i<<1);
                rnd_gen_gauss2(ptr, ptr+1);
        }
        if (length & 1) {
                double holder;
                rnd_gen_gauss2(place + length - 1, &holder);
        }
        for (int i=0; i<length; i++) {
                place[i] *= sigma;
        }
}

void
CEX_generate_gauss_vector(vec_t *vec, double sigma)
{
        double holder;
        rnd_gen_gauss2(&(vec->x), &(vec->y));
        rnd_gen_gauss2(&(vec->z), &holder);
        Vec3_MULTO(*vec, sigma);
}
