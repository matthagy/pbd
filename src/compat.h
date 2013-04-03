/* -*- Mode: c -*-
 * compat.h - Compatibility helpers
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

#ifndef _COMPAT_H
#define _COMPAT_H

#ifdef USING_MPI_LAM
# define GET_MPI_STATUS_BYTES(stat) ((stat).st_length)
#endif

#ifdef USING_MPI_MVAPICH2
# define GET_MPI_STATUS_BYTES(stat) ((stat).count)
#endif

#endif /* _COMPAT_H */
