/* -*- Mode: c -*-
 * opt.h - GCC compiler optimizations/options
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

#ifndef OPT_H_
#define OPT_H_

#ifdef __GNUC__
#  define GCC_ATTRIBUTE(ATTR) __attribute__(ATTR)
#else
#  define GCC_ATTRIBUTE(ATTR) 
#endif

#ifndef NO_GCC_OPTIMIZATIONS
/* load effective location before use */
#  define prefetch(p) __builtin_prefetch(p, 1, 3)
/* branch prediction hints */
#  define likely(x)   __builtin_expect(!!(x), 1)
#  define unlikely(x) __builtin_expect(!!(x), 0)
#else
#  define prefetch(p) do {} while(0)
#  define likely(x)   (x)
#  define unlikely(x) (x)
#endif

/* compiler optimization requiring unaliased arrays
 * needed for Fortran-like array processing */
#ifndef NRESTRICT
# define CEX_RESTRICT restrict
#else
# define CEX_RESTRICT
#endif

#endif /*OPT_H_*/
