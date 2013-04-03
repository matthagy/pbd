##-*- Mode: python -*-
## msg.py - Serialization for invoking commands in simulation process
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

from numpy import array

from jamenson.runtime.multimethod import MultiMethod, defmethod
from jamenson.runtime.atypes import anytype, typep
from jamenson.runtime.atypes.ptypes import integer_type

from .xtypes import (integer_type, positive_integer_type, single_char_type, seq_type,
                     string_or_seq_type, number_type)


def make_writing_message(command_name, frmt='', *args):
    return WritingMessage(command_name).write_frmt(frmt, *args)

def uint2str(value):
    assert typep(value, positive_integer_type)
    assert 0<=value<=0x100000000L
    return ''.join(chr((value>>(8*i)) & 0xFF) for i in xrange(3,-1,-1))

def str2uint(s):
    assert len(s)==4
    return sum(ord(c)<<(8*i) for i,c in zip(xrange(3,-1,-1), s))


class WritingMessage(object):
    '''serialize arguements for a command invocation
    '''

    def __init__(self, command_name=None):
        self.buffer = []
        if command_name:
            assert isinstance(command_name, str)
            self.write_char_array(command_name)

    def prepare(self):
        bytes = ''.join(self.buffer)
        self.buffer = [bytes]
        return bytes

    def write_char(self, c):
        assert isinstance(c, str)
        assert len(c)==1
        self.buffer.append(c)
        return self

    def write_uint(self, op):
        self.buffer.append(uint2str(op))
        return self

    def write_int(self, op):
        assert abs(op) <= sys.maxint
        sign = 1 if op>=0 else 0
        return self.write_uint(abs(op) | (sign<<31))

    def write_array(self, seq, el_writer):
        if not typep(seq, string_or_seq_type):
            seq = list(seq)
        self.write_uint(len(seq))
        for el in seq:
            el_writer(self, el)
        return self

    def write_char_array(self, s):
        return self.write_array(s, self.__class__.write_char)

    def write_int_array(self, arr):
        return self.write_array(arr, self.__class__.write_int)

    def write_uint_array(self, arr):
        return self.write_array(arr, self.__class__.write_uint)

    def write_double(self, op):
        return self.write_char_array('%.10e' % op)

    def write_double_array(self, arr):
        return self.write_array(arr, self.__class__.write_double)

    def write_vec(self, v):
        x,y,z = v
        self.write_double(x)
        self.write_double(y)
        self.write_double(z)
        return self

    def write_vec_array(self, arr):
        return self.write_array(arr, self.__class__.write_vec)

    def write_submsg(self, msg):
        return self.write_char_array(''.join(msg.buffer))

    def write_frmt(self, frmt, *args):
        assert isinstance(frmt, str)
        for frmt,arg in zip(frmt, args):
            write_frmt(self, frmt, arg)
        return self

    def write_msg(self, msg):
        msg.buffer.append(self.prepare())

write_frmt = MultiMethod(name='write_frmt',
                         signature='writer,frmt,value',
                         doc='''direct formatted writing
                         ''')

@defmethod(write_frmt, [anytype, "u", positive_integer_type])
def meth(writer, frmt, i):
    writer.write_uint(i)

@defmethod(write_frmt, [anytype, "i", integer_type])
def meth(writer, frmt, i):
    writer.write_int(i)

@defmethod(write_frmt, [anytype, "U", seq_type])
def meth(writer, frmt, seq):
    writer.write_uint_array(seq)

@defmethod(write_frmt, [anytype, "I", seq_type])
def meth(writer, frmt, seq):
    writer.write_int_array(seq)

@defmethod(write_frmt, [anytype, "s", str])
def meth(writer, frmt, s):
    writer.write_char_array(s)

@defmethod(write_frmt, [anytype, "f", number_type])
def meth(writer, frmt, seq):
    writer.write_double(seq)

@defmethod(write_frmt, [anytype, "F", seq_type])
def meth(writer, frmt, seq):
    writer.write_double_array(seq)

@defmethod(write_frmt, [anytype, "v", seq_type])
def meth(writer, frmt, seq):
    writer.write_vec(seq)

@defmethod(write_frmt, [anytype, "V", seq_type])
def meth(writer, frmt, seq):
    writer.write_vec_array(seq)

@defmethod(write_frmt, [anytype, "m", WritingMessage])
def meth(writer, frmt, msg):
    writer.write_submsg(msg)

@defmethod(write_frmt, [anytype, "o", anytype])
def meth(writer, frmt, op):
    op.write_msg(writer)

class NamedItems(object):
    '''Data structure for writing named elements in WritingMessage
    '''

    def __init__(self, items):
        self.items = items

    def write_msg(self, msg):
        for name,frmt,op in self.items:
            msg.write_char_array(name)
            msg.write_frmt(frmt, op)

class StructArray(object):
    '''Data structure for writing arrays of NamedItems
    '''

    def __init__(self, elements):
        self.elements = elements

    def write_msg(self, msg):
        msg.write_array(self.elements, self.write_element)

    @staticmethod
    def write_element(msg, el):
        if not isinstance(el, NamedItems):
            el = NamedItems(el)
        el.write_msg(msg)


class ReadingMessage(object):
    '''interpret bytes returned from a command invocation
       in wrapped process
    '''

    def __init__(self, bytes):
        self.buffer = list(bytes)
        self.buffer.reverse()

    def read_char(self):
        return self.buffer.pop()

    def read_uint(self):
        return str2uint(''.join(self.read_char() for _ in '1234'))

    def read_int(self):
        base = self.read_uint()
        sign = (base >> 31) & 1
        ibase = base & ((1<<31)-1)
        return ibase if sign==1 else -ibase

    def read_array(self, el_reader):
        return [el_reader(self) for i in xrange(self.read_uint())]

    def read_char_array(self):
        return ''.join(self.read_array(self.__class__.read_char))

    def read_int_array(self):
        return array(self.read_array(self.__class__.read_int))

    def read_double(self):
        return float(self.read_char_array())

    def read_double_array(self):
        return array(self.read_array(self.__class__.read_double))

    def read_vec(self):
        return array([self.read_double() for _ in 'xyz'])

    def read_vec_array(self):
        return array(self.read_array(self.__class__.read_vec))

    null = object()
    def req_eofp(self):
        if self.buffer:
            raise RuntimeError("expected end of message with %d character remaining" %
                               (len(self.buffer),))
        return self.null

    def read_frmt(self, frmt):
        return [op for op in [read_frmt(self,f) for f in frmt]
                if op is not self.null]

read_frmt = MultiMethod(name='read_frmt',
                        signature='reader,frmt',
                        doc='''direct formatted reading
                         ''')

@defmethod(read_frmt, [anytype, "u"])
def meth(reader, frmt):
    return reader.read_uint()

@defmethod(read_frmt, [anytype, "i"])
def meth(reader, frmt):
    return reader.read_int()

@defmethod(read_frmt, [anytype, "s"])
def meth(reader, frmt):
    return reader.read_char_array()

@defmethod(read_frmt, [anytype, "I"])
def meth(reader, frmt):
    return reader.read_int_array()

@defmethod(read_frmt, [anytype, "f"])
def meth(reader, frmt):
    return reader.read_double_array()

@defmethod(read_frmt, [anytype, "v"])
def meth(reader, frmt):
    return reader.read_vec()

@defmethod(read_frmt, [anytype, "V"])
def meth(reader, frmt):
    return reader.read_vec_array()

@defmethod(read_frmt, [anytype, "x"])
def meth(reader, frmt):
    return reader.req_eofp()

class ReadingMessageList(object):

    def __init__(self, msgs):
        self.msgs = list(msgs)

    def read_frmt(self, frmt):
        return [msg.read_frmt(frmt) for msg in self.msgs]
