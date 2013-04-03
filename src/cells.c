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

#include "cells.h"

array_t *CEX_jcells=NULL;
array_t *CEX_surface_junctions=NULL;
array_t *CEX_line_junctions=NULL;
array_t *CEX_point_junctions=NULL;
cell_t *CEX_this_cell=NULL;
