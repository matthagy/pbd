/* -*- Mode: c -*-
 * cell.h - Division of cartersian space into rectangular cells
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

#ifndef _CELLS_H
#define _CELLS_H

#include "opt.h"
#include "array.h"
#include "comm.h"

typedef struct {
        comm_t *comm;
        vec_t min_extent;
        vec_t max_extent;
} cell_t;

extern cell_t *CEX_this_cell; /* information explaining this cell */

extern array_t *CEX_jcells; /* array of junctioning cells */
#define HAVE_JUNCTIONS() likely(ARR_LENGTH(CEX_jcells)!=0)


#define JCELL_FOREACH(ptr, counter)                     \
        ARR_FOREACH(cell_t, CEX_jcells, ptr, counter)

#define CHK_1EXTENT(POS, AX, MN, MX)                 \
        ((POS).  AX  >= (MN).  AX &&                 \
         (POS).  AX < (MX).  AX)

#define CHK_EXTENT(POS, MN, MX)                         \
        (CHK_1EXTENT(POS, x, MN, MX) &&                 \
         CHK_1EXTENT(POS, y, MN, MX)  &&                \
         CHK_1EXTENT(POS, z, MN, MX))

#define CELL_CONTAINS(CELLP, POS)               \
        CHK_EXTENT(POS, (CELLP)->min_extent, (CELLP)->max_extent)

static inline cell_t * find_jcell_containing(vec_t position) 
        GCC_ATTRIBUTE((always_inline));

static inline cell_t * 
find_jcell_containing(vec_t position)
{
        cell_t *jcellp; int counter;
        JCELL_FOREACH(jcellp, counter) {
                if (unlikely(CELL_CONTAINS(jcellp, position)))
                        return jcellp;
        }
        return NULL;
}

/* Cell Junction Information
 * There are three ways that cells can touch each other
 *
 *    o Surface Junctions (6)
 *    o Line Junctions (12)
 *    o Point Junctions (8)
 *
 * This means each cell can have a max of 26 junctions, although with periodic bondary conditions
 * cells may be junctions along multiple modes
 *
 *                        -----\           
 *                              -----\     
 *                                    ---+ 
 *                               -----/  | 
 *                          ----/        | 
 *----\                ----/             | 
 *     ----\  |  -----/                  | 
 *          --+-/                        |                      /---
 *            |                          |                  /---    \---
 *            |    Point                 |              /---            \----
 *            |    Junction             |          /---                     \---
 *            |                          |      /---                             \---
 *            |                         o|o /---       Reference Cell                \-+--
 *            |                        --+--                                        /--|  \----
 *            |                    ---/ o|o \----                               /---   |       \----
 *            |               ----/      |       \----                      /---      /+            \----
c *--\         |           ---/           |            \---              /---      /--- |                 \----
 *   ---\     |      ----/               |                \----     /---      /---     |                      \---
 *       ---\ |  ---/                    |                     \----       /--        /+                    /---|
 *           -+-/                        |                       |  \--/---       /--- |                 /--    |
 *                                       |                       | /---  \----/---    /+             /---       |
 *                                       |                       +-        /--\----/-- |          /--           |
 *                                       +--                     |     /---     /--\---+      /---              |
 *                                      ----\---                 | /---     /---     /-+--- +-                  |
 *                                  ---/---\----\----            +-      /--     /---     \-+--                 |
 *                              ---/       \---\---- \----       |    /--     /--           |  \-----           |
 *                          ---/               \--------  \---   | /--    /---              |        \----      |
 *                      ---/     Line Junction      \-------  \--+-    /--                  |             \---- |
 *                   +-/                                 \------ | /---                     |                  X+-
 *                   |\--                                    \---+--         Surface        |                /-
 *                   |   \---                               ---/ |  \---     Junction       |             /--
 *                   |       \--                         --/     |      \---                |           /-
 *                              \--                  ---/        |          \---            |        /--
 *                                 \--            --/                           \---        |      /-
 *                                    \---    ---/                                  \---    |   /--
 *                                        \--/                                          \-- | /-
 *                                          |                                               --
 *                                          |                     
 *                                          |                     
 */

/* Xxx These are used to index elements of vec_t struct
 * have to correspond to write struct elements
 */
#define AXIS_X 0
#define AXIS_Y 1
#define AXIS_Z 2
#define INDEX_AXIS(vec_p, axis) (((double *)(vec_p))[axis])

typedef struct {
        cell_t *cell;
        int axis; 
        int dir;
} surface_junction_t;

typedef struct {
        cell_t *cell;
        int axis;
        double offset1, offset2;
} line_junction_t;

typedef struct {
        cell_t *cell;
        vec_t offset;
} point_junction_t;

extern array_t *CEX_surface_junctions;
extern array_t *CEX_line_junctions;
extern array_t *CEX_point_junctions;


#endif /* _CELLS_H */
