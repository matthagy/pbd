##-*- Mode: python -*-
## debug.py - Simple debugging helpers
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

from __future__ import absolute_import

import sys

level = 1

def say(msg='', *args):
    msg = str(msg)
    msg = msg%args if args else msg
    for line in msg.split('\n'):
        print >>sys.stderr, line
    sys.stderr.flush()

def set_level(new_level=None):
    global level
    if new_level is None:
        new_level = 1
    assert isinstance(new_level, (int,long))
    level = new_level

def msg(msg='', *args):
    if level >= 1:
        say(msg, *args)

def warn(msg='', *args):
    if level >= 0:
        say(msg, *args)

def error(_msg='', *args):
    say(_msg, *args)
    sys.exit(1)


