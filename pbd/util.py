##-*- Mode: python -*-
## util.py - Assorted utilities
## --------------------------------------------------------------------------
## Copyright (C) 2009, Matthew Hagy (hagy@gatech.edu)
## All rights reserved.
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY# without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


from numpy import *

def tessellate(size, positions, depth=1,
               cube_offsets=mgrid[0:2:1,0:2:1,0:2:1].swapaxes(0,3).reshape(8,3)):
    if not depth:
        return size, positions
    return tessellate(2*size,
                      concatenate(list(positions + offset * size
                                       for offset in cube_offsets)),
                      depth-1)
