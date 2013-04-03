##-*- Mode: python -*-
## xtypes.py - Various algebric types
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


from numpy import ndarray, integer as numpy_integer, floating as numpy_float

from jamenson.runtime.atypes import (typep, as_optimized_type, anytype,
                                     intersection,
                                     IsInstanceType, IsType)
from jamenson.runtime.atypes.ptypes import integer_type as pinteger_type

integer_type = as_optimized_type((pinteger_type, numpy_integer))
positive_integer_type = as_optimized_type(intersection(integer_type, lambda x : x >= 0))
single_char_type = as_optimized_type(intersection(str, lambda x: len(x)==1))
seq_type = as_optimized_type((list,tuple,ndarray))
string_or_seq_type = as_optimized_type((seq_type, str))
float_type = as_optimized_type((float, numpy_float))
number_type = as_optimized_type((integer_type, float_type))
