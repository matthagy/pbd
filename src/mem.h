/* -*- Mode: c -*-
 * mem.h - Memory managment utilities
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


#ifndef MEM_H
#define MEM_H

#include "opt.h"

#ifndef MACOSX
#  include <alloca.h>
#endif

#ifdef MACOSX
#  include <strings.h>
#endif


#define bzero(p,sz) (memset((p), '\0', (sz)), (void)0)

void *CEX_malloc(size_t) GCC_ATTRIBUTE((malloc));
void *CEX_realloc(void *, size_t) GCC_ATTRIBUTE((malloc));
void CEX_free(void *);

#define XNEW(TP, SZ) ((TP *)CEX_malloc(sizeof(TP)*(SZ)))
#define XRESIZE(TP, PTR, SZ) ((TP *)CEX_realloc((void *)(PTR), sizeof(TP)*(SZ)))
#define XBZERO(TP, PTR, SZ) bzero((void *)PTR, sizeof(TP)*(SZ))
#define XMEMCPY(TP, DEST, SRC, SZ) memcpy((void *)(DEST), (void *)(SRC), sizeof(TP)*SZ)
#define XALLOCA(TP, SZ) ((TP *)alloca(sizeof(TP)*(SZ)))
                         
void *CEX_aligned_malloc(size_t bytes, size_t alignment) GCC_ATTRIBUTE((malloc));
void *CEX_aligned_realloc(void *, size_t bytes, size_t alignment) GCC_ATTRIBUTE((malloc));
void CEX_aligned_free(void *);


#endif /* MEM_H */
