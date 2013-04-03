/* -*- Mode: c -*-
 * xperiodic.h - Basic macros for cartesian periodic bondary conditions
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


#ifndef _XPERIODIC_H
#define _XPERIODIC_H

#include "opt.h"
#include "vector.h"

#define XPERIODIZE_SEPARATION(vp, sz, hz)    \
        if (unlikely(vp > hz)) {            \
                vp -= sz;                   \
        } else if (unlikely(vp < -hz)) {    \
                vp += sz;                   \
        }

#define XPERIODIC_SEPARATION_VECTOR(r, pos_i, pos_j, box_size, box_sz_half) do { \
        Vec3_SUB(r, pos_j, pos_i);                                               \
        XPERIODIZE_SEPARATION(r.x, box_size.x, box_sz_half.x);                   \
        XPERIODIZE_SEPARATION(r.y, box_size.y, box_sz_half.y);                   \
        XPERIODIZE_SEPARATION(r.z, box_size.z, box_sz_half.z);                   \
        } while(0)

#define XPERIODIZE_COMPONENT(vp, sz)   \
        vp -= (vp>sz) * sz,           \
        vp += (vp<0.0) * sz

#define XPERIODIZE_LOCATION(v, size) do{                 \
                XPERIODIZE_COMPONENT((v).x, (size).x);   \
                XPERIODIZE_COMPONENT((v).y, (size).y);   \
                XPERIODIZE_COMPONENT((v).z, (size).z);   \
        } while (0)

#endif /* _XPERIODIC_H */
