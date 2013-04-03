##-*- Mode: python -*-
## cex.py - Interface to cex (C-Extension) program
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
from __future__ import with_statement

import os
import signal
from subprocess import Popen
import random
import atexit
import time

from HH2.pathutils import DirPath
from jamenson.runtime.atypes import typep

import setup
from .msg import (make_writing_message, ReadingMessage, ReadingMessageList,
                  uint2str, str2uint)
from .xtypes import positive_integer_type


class CexInterface(object):
    '''Interface to cex (C-extension) process for parallel Brownian Dynamics
       simulation
    '''

    # # # #
    # API #
    # # # #

    @classmethod
    def create(self):
        '''create an interface to a single threaded simulation process
        '''
        return self.create_ex([setup.paths.cex.abspath()])

    @classmethod
    def create_mpi(self, nproc=1, mpiargs=[], mpirun=None, nodes=None):
        '''create an interface to a process launched over MPI
        '''
        assert typep(nproc, positive_integer_type)
        if mpirun is None:
            mpirun = 'mpirun'
        if not mpirun.startswith('./') or mpirun.startswith('/'):
            mpirun = which(mpirun)
        if isinstance(mpiargs, str):
            mpiargs = mpiargs.split()
#         if nodes is None:
#             nodes = get_pbs_nodes()
        basedir = make_tmp_directory()
#         machine_path = basedir.child('machines')
#         assert len(nodes) >= nproc
#         with open(machine_path, 'w') as fp:
#             for node in nodes[:nproc]:
#                 print >>fp, node
        return self.create_ex([mpirun] + list(mpiargs) +
                              [#'-machinefile', machine_path,
                               '-np', nproc, setup.paths.cex.abspath()],
                              basedir=basedir)

    shutting_down = False
    def shutdown(self):
        '''shutdown wrapped cex process
        '''
        if not self.check_active():
            return
        if self.shutting_down:
            return
        self.shutting_down = True
        self.on_each_slave_async(make_writing_message('exit')).read_frmt('x')
        self.do_command(0, make_writing_message('exit')).read_frmt('x')
        self.active = False
        force_proc_exit(self.proc)

    size = None
    def get_size(self):
        '''retrieve the number of threads that compose simulation
        '''
        if self.size is None:
            self.size, = self.perform_command(0, make_writing_message('poll_size', 'u'),
                                              ).read_frmt("ux")
        return self.size

    def perform_command(self, rank, msg):
        '''send a writing msg to the process of specified rank and
           returns the reading msg
        '''
        self.req_active()
        return self.do_command(rank, msg)

    def map_slave(self, msgs):
        '''map a sequence of writing msgs across all slave
           and return the resulting reading msgs
        '''
        self.req_active()
        return ReadingMessageList([self.do_command(rank, imsg)
                for rank,imsg in self.iter_slave_msgs(msgs)])

    def on_each_slave(self, msg):
        return self.map_slave([msg] * (self.get_size()-1))

    def map_all(self, msgs):
        self.req_active()
        if not isinstance(msgs, list):
            msgs = list(msgs)
        return ReadingMessageList([self.do_command(0, msgs[0])] +
                                  self.map_slave(msgs[1:]).msgs)

    def on_each(self, msg):
        return self.map_all([msg] * self.get_size())

    def map_slave_async(self, msgs):
        self.req_active()
        self.map_slave_async_send(msgs)
        return self.map_slave_async_recv()

    def on_each_slave_async(self, msg):
        return self.map_slave_async([msg] * (self.get_size()-1))

    def map_all_async(self, msgs):
        self.req_active()
        if not isinstance(msgs, list):
            msgs = list(msgs)
        self.map_slave_async_send(msgs[1:])
        return ReadingMessageList([self.do_command(0,msgs[0])] +
                                  self.map_slave_async_recv().msgs)

    def on_each_async(self, msg):
        return self.map_all_async([msg] * self.get_size())


    # # # # # # #
    # Internals #
    # # # # # # #

    # Process Management
    @classmethod
    def create_ex(cls, launch_args, basedir=None):
        #print launch_args
        if basedir is None:
            basedir = make_tmp_directory()
        try:
            write_fifo_path = basedir.child('python2c-fifo')
            read_fifo_path = basedir.child('c2python-fifo')
            os.mkfifo(write_fifo_path)
            os.mkfifo(read_fifo_path)
            proc = Popen(map(str, launch_args + [write_fifo_path, read_fifo_path]))
            try:
                write_fifo, read_fifo = open_fifos(write_fifo_path, read_fifo_path)
            except IOError:
                force_proc_exit(proc)
                raise
            else:
                return cls(proc, write_fifo, read_fifo)
        except:
            remove_directory(basedir)
            raise

    def __init__(self, proc, write_fp, read_fp):
        self.proc = proc
        self.write_fp = write_fp
        self.read_fp = read_fp
        atexit.register(self.shutdown)

    active = True
    def req_active(self):
        if not self.check_active():
            raise RuntimeError('proc exited early with code %r' %
                                (self.proc.poll(),))

    def check_active(self):
        if self.active:
            self.active = self.proc.poll() is None
        return self.active

    # IO
    def write(self, bytes):
        try:
            self.write_fp.write(bytes)
            self.write_fp.flush()
        except IOError:
            self.req_active()
            raise

    def read(self, nbytes):
        try:
            bytes = self.read_fp.read(nbytes)
            if len(bytes) != nbytes:
                raise EOFError
            return bytes
        except (IOError, EOFError):
            self.req_active()
            raise

    def write_uint(self, value):
        self.write(uint2str(value))

    def read_uint(self):
        return str2uint(self.read(4))

    def do_command(self, rank, msg):
        self.write_uint(rank)
        bytes = msg.prepare()
        self.write_uint(len(bytes))
        self.write(bytes)
        return ReadingMessage(self.read(self.read_uint()))

    def check_slave_msgs(self, msgs):
        if not isinstance(msgs, list):
            msgs = list(msgs)
        if len(msgs) != self.get_size()-1:
            raise ValueError("need %d msgs; given %d" % (self.get_size()-1, len(msgs)))
        return msgs

    def iter_slave_msgs(self, msgs):
        for i,msg in enumerate(self.check_slave_msgs(msgs)):
            yield i+1,msg

    def map_slave_async_send(self, msgs):
        for rank,msg in self.iter_slave_msgs(msgs):
            self.do_command(0, make_writing_message('send_msg', 'im', rank, msg))

    def map_slave_async_recv(self):
        return ReadingMessageList(self.do_command(0, make_writing_message('recv_msg', 'i', rank))
                                  for rank in xrange(1, self.get_size()))



# # # # # # # # # # #
# Helper Functions  #
# # # # # # # # # # #

def get_path():
    return map(DirPath, os.environ['PATH'].split(':'))

def which_iter(name, paths=None):
    for path in paths or get_path():
        path = path.child(name)
        if path.exists():
            yield path

def which(name, paths=None):
    paths = list(which_iter(name, paths))
    if not paths:
        raise RuntimeError("couldn't find %s in path" % (name,))
    return paths[0]

def make_tmp_directory():
    path = DirPath('/tmp/%sP%dS%X-c-python-bridge' %
                   (os.environ['USER'], os.getpid(),
                    random.randrange(0xffffffff)))
    atexit.register(lambda : remove_directory(path))
    path.reqdir()
    return path

def remove_directory(path):
    try:
        path.rmdir(recursive=True)
    except OSError:
        pass

def force_proc_exit(proc):
    exited = lambda : proc.poll() is not None
    if exited():
        return
    for i in xrange(5):
        if exited():
            return
        time.sleep(0.05)
    if exited():
        return

    try:
        os.kill(proc.pid, signal.SIGTERM)
    except OSError:
        pass
    else:
        for i in xrange(5):
            if exited():
                return
            time.sleep(0.02)
        if exited():
            return
        try:
            os.kill(proc.pid, signal.SIGKILL)
        except OSError:
            pass

def open_fifos(write_fifo_path, read_fifo_path):
    #use alarm to abort blocking IO
    old_handler = signal.signal(signal.SIGALRM, lambda sig, frame: None)
    try:
        signal.alarm(200)
        try:
            #IOError will be raised from interrupt if we recieve SIGALRM
            return open(write_fifo_path,'wb'), open(read_fifo_path,'rb')
        finally:
            signal.alarm(0)
    finally:
        signal.signal(signal.SIGALRM, old_handler)

def get_pbs_nodes(ncpus=8):
    return reduce(lambda a,b: a+b, ([node]*ncpus for node in get_pbs_machines()))

def get_pbs_machines():
    return filter(None, (line.strip() for line in open(os.environ['PBS_NODEFILE'])))
