/* -*- Mode: c -*-
 * constants.h - Model and Physical Constants
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

#ifndef _CONSTANTS_H
#define _CONSTANTS_H

/* Units: Meter-Kilogram-Second (MKS) */
#define CEX_kB 1.38e-23 /* J/K */

/* meters */
#define CEX_nm 1e-9
#define CEX_mcm 1e-6
#define CEX_mm 1e-3

/* seconds */
#define CEX_ps 1e-12
#define CEX_ns 1e-9
#define CEX_mcs 1e-6
#define CEX_ms 1e-3
#define CEX_m 1

/* force */
#define CEX_pN 1e-12

#define CEX_R_particle (0.5 * 135 * CEX_nm)

#endif /* _CONSTANTS_H */
