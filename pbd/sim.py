##-*- Mode: python -*-
## cex.py - High Level Interface to Simulation
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
from __future__ import division

import sys
import os
import time
import struct
import itertools

from numpy import *
from numpy.random import RandomState

from .debug import msg, warn, error
from . import state
from . import cells
from . import constants
from . import forcefield
from . msg import make_writing_message, NamedItems, StructArray
from . import cex


class Simulator(object):
    '''abstracts away interface to simulation process
    '''

    # # # #
    # API #
    # # # #

    @classmethod
    def create(cls, cexinf, parameters=None, configuration=None,
               divisions=None, random_seed=None):
        '''create a Simulator from an uninitialized CexInterface
        '''
        if parameters is None:
            parameters = state.Parameters()
        if configuration is None:
            configuration = state.Configuration(positions=array([]).reshape(0,3))
        assert isinstance(parameters, state.Parameters)
        assert isinstance(configuration, state.Configuration)
        initialize(cexinf, parameters, configuration, divisions, random_seed)
        return cls(cexinf, parameters.time_step, configuration.time, parameters)

    def simulate(self, n_cycles, max_c_cycles=2500):
        '''simulate n_cycles integration cycles. max_c_cycles specifies the number
           of steps to perform in a single command invocation the childr process.
        '''
        n_max,extra = divmod(n_cycles, max_c_cycles)
        for i in xrange(n_max):
            self.simulate_cycles(max_c_cycles)
        if extra:
            self.simulate_cycles(extra)

    def get_configuration(self):
        wall_time = time.time()
        return state.Configuration(time=self.get_time(),
                                   wall_time=wall_time,
                                   positions=self.get_positions())

    def get_time(self):
        return self.start_time + self.simulated_cycles * self.time_step

    def get_positions(self):
        '''retrieve the positions of all particles in the simulation

           return particle positions in order in which particles where
           first give.  allows for tracing position of individual positions
        '''
        def fix_array(arr, n):
            if not len(arr):
                return zeros((0,n) if n else (0,), int)
            return arr
        positions, tags = map(concatenate,
                              zip(*list((fix_array(positions, 3),fix_array(tags, 0))
                                        for positions,tags in
                                        self.cexinf.on_each_async(
                           make_writing_message('collect_thread_positions_and_tags')).read_frmt('VIx'))))
        assert len(set(tags)) == len(tags)
        tags,positions = zip(*sorted(zip(tags, positions)))
        positions = array(positions)
        box_size = array(self.parameters.box_size)
        if self.parameters is not None:
            if any(positions >= box_size):
                #print where(positions >= box_size)
                raise RuntimeError('invalid positions %s' % ', '.join(
                    '%s (%d) ' % (positions[i], i) for i in where(positions >= box_size)[0]))
        return positions

    def update_neighbors(self):
        '''update neighbor information across all threads
        '''
        self.cexinf.on_each_async(make_writing_message('thread_update_neighbors')).read_frmt('x')

    def get_state(self):
        '''retrieve the internal state of each thread.  largely only useful for
           debugging
        '''
        wall_time = time.time()
        def read_neighbors(arr):
            assert len(arr) % 2 == 0
            return arr.reshape(len(arr)>>1, 2)
        return state.SimulationState(
            time=get_time(),
            wall_time=wall_time,
            threads=[state.ThreadState(
                         positions=positions, tags=tags,
                         internal_neighbors=read_neighbors(internal_neighbors),
                         external_neighbors=read_neighbors(external_neighbors))
                     for positions,tags,internal_neighbors,external_neighbors in
                     self.cexinf.on_each_async(make_writing_message('collect_thread_state')).read_frmt('VIIIx')])

    # # # # # # #
    # Internals #
    # # # # # # #

    def __init__(self, cexinf, time_step, start_time=0, parameters=None):
        self.cexinf = cexinf
        self.time_step = time_step
        self.start_time = start_time
        self.simulated_cycles = 0
        self.parameters = parameters

    def simulate_cycles(self, steps):
        self.cexinf.map_slave_async_send([make_writing_message("slave_simulation_loop")]*(self.cexinf.get_size()-1))
        self.cexinf.perform_command(0, make_writing_message('master_simulate_cycles', 'i', steps))
        self.simulated_cycles += steps
        self.cexinf.map_slave_async_recv()


# # # # # # # # # # # # # #
# Initialization Routines #
# # # # # # # # # # # # # #

def initialize(cexinf, parameters, configuration, divisions, random_seed):
    '''initialize a cex process (through cexinf) for the
       simulation of the specified system
    '''
    initialize_thread_names(cexinf)
    initialize_system(cexinf, parameters)
    initialize_random(cexinf, random_seed)
    thread_cells = create_cells(array(parameters.box_size),
                                configuration.positions, divisions, cexinf.get_size())
    setup_comm_rules(thread_cells)
    initialize_thread_cells(cexinf, thread_cells)

def initialize_thread_names(cexinf):
    cexinf.map_all_async(make_writing_message('set_thread_name', 's',
                                              'master' if i==0 else 'slave%d' % i)
                         for i in xrange(cexinf.get_size()))


def initialize_system(cexinf, parameters):
    '''setup-thread independent state
    '''
    kT = constants.kB * parameters.temperature
    cexinf.on_each_async(make_writing_message("initialize_system", "o",
                           NamedItems([
         #Size of Full Ensemble Cartesian Space
         ["box_size", "v", parameters.box_size],
         #Integration Constants
         ["T", "f", parameters.temperature],
         ["dt", "f", parameters.time_step],
         ["dU_max", "f", parameters.dU_max * kT],
         # friction from Stokes-Einstein relationship
         ["fric_gamma", "f", 6 * pi * parameters.eta_solv * constants.R_particle],
         #Force Field
         ["force_update", "i", parameters.force_update_rate],
         ["r_pair_cutoff", "f", parameters.r_potential_cutoff],
         ["pair_potential", "o", LinterpWriter(kT * parameters.pair_potential.make_potential_table(
                                                        parameters.linterp_r_min,
                                                        parameters.r_potential_cutoff * 1.05,
                                                        parameters.linterp_size))],
         #evaluate_forces relies on force linterp being normalized
         #for vector length
         ["pair_force", "o", LinterpWriter(scale_force_table(
                              kT * parameters.pair_potential.make_force_table(
                                     parameters.linterp_r_min,
                                     parameters.r_potential_cutoff * 1.05,
                                     parameters.linterp_size)))],
         #Neighbor Lists
         ["r_neighbor", "f", parameters.r_neighbor]]))
      ).read_frmt("x")


class LinterpWriter(object):

    def __init__(self, linterp):
        self.linterp = linterp

    def write_msg(self, msg):
        msg.write_frmt("o", NamedItems([
            ["x_min", "f", self.linterp.x_min],
            ["x_prec", "f", self.linterp.x_prec],
            ["table", "F", self.linterp.y]]))

def scale_force_table(table):
    '''evaluate_forces relies on force linterp being normalized
       for vector length
    '''
    r  = table.x_min + table.x_prec * arange(len(table.y))
    r[where(r==0)] = 1 # just do this to prevent zero-division error
                       # this distance (r=0) is unimportant for the table itself
    table.y /= r
    return table


def initialize_random(cexinf, random_seed):
    if random_seed is None:
        random_seed = generate_seed()
    random_seed = int(random_seed)
    msg('initializing random state with seed=0x%X', random_seed)
    rnd = RandomState(random_seed)
    cexinf.map_all_async(make_writing_message('initialize_random', 'u', seed)
                         for seed in rnd.randint(0xfffffff, size=cexinf.get_size())
                         ).read_frmt('x')

def generate_seed():
    frmt = '@I'
    number, = struct.unpack(frmt, os.urandom(struct.calcsize(frmt)))
    return int(number & sys.maxint)

def create_cells(box_size, positions, divisions, world_size):
    '''partition simulation space into cells using divisions
       divide positions among these cells, and tag them with
       their index in the orignal positions array
    '''
    if divisions is None:
        divisions = world_size
    div_dimensions = cells.get_divs(divisions)
    n_threads = multiply.reduce(list(div_dimensions))
    if n_threads != world_size:
        error("bad number of threads %d with divisions %s for %d theads",
              n_threads, div_dimensions, world_size)
    msg('initializing cells for n=%d with dimensions %s',
        n_threads, div_dimensions)
    thread_cells = cells.partition_positions(box_size, positions, div_dimensions)
    for cell in thread_cells:
        cell.junctioned_cells = list(set(junction.cell for junction in cell.junctions))
        cell.jcell_indices = dict((jcell,i) for i,jcell in enumerate(cell.junctioned_cells))
    N_particles = array(list(len(cell.positions) for cell in thread_cells))
    msg('particle distribution min=%d, max=%d, mean=%.1f, std=%.1f',
        N_particles.min(), N_particles.max(),
        N_particles.mean(), N_particles.std())
    #relate tags to positions
    tag_map = dict(zip(map(tuple, positions), xrange(len(positions))))
    assert len(tag_map) == len(positions)
    for cell in thread_cells:
        cell.tags = array(list(tag_map[pos] for pos in map(tuple, cell.positions)))
    all_tags = set(tag_map.itervalues())
    for cell in thread_cells:
        stags = set(cell.tags)
        assert all_tags >= stags
        all_tags -= stags
    assert not all_tags
    return thread_cells

def setup_comm_rules(thread_cells):
    '''determine order in which threads send and recieve
       data between each other to prevent dead lock without
       global synchronization
    '''
    for rank_i,cell in enumerate(thread_cells):
        cell.rank = rank_i
    rounds = cells.calculate_communication_rules(thread_cells)
    for cell in thread_cells:
        cell.comm_rules = []
    next_tag = iter(itertools.count(1)).next
    for round in rounds:
        for link in round:
            msg_tag = next_tag()
            link.cell_i.comm_rules.append(['send', link.cell_j, msg_tag])
            link.cell_j.comm_rules.append(['recv', link.cell_i, msg_tag])
            msg_tag = next_tag()
            link.cell_j.comm_rules.append(['send', link.cell_i, msg_tag])
            link.cell_i.comm_rules.append(['recv', link.cell_j, msg_tag])

def initialize_thread_cells(cexinf, thread_cells):
    '''setup up thread specific information for the cell
       each thread simulates and the junctions between cells.
    '''
    cexinf.map_all_async(make_writing_message('initialize_cell_state', 'o',
                           NamedItems([
                               ['min_extent', 'v', cell.extent.min_extent],
                               ['max_extent', 'v', cell.extent.max_extent],
                               ['positions', 'V', cell.positions],
                               ['tags', 'I', cell.tags]]))
                      for cell in thread_cells).read_frmt('x')
    inst_map = dict(send=1, recv=2)
    cexinf.map_all_async(make_writing_message('initialize_cell_comm', 'o',
                                  NamedItems([
                                    ['comms', 'o', StructArray(
                                       [['comm_rank', 'i', jcell.rank]]
                                       for jcell in  cell.junctioned_cells)],
                                    ['comm_rules', 'o', StructArray(
                                      [['inst', 'i', inst_map[inst]],
                                       ['comm_index', 'i', cell.jcell_indices[jcell]],
                                       ['tag', 'i', msg_tag]]
                                      for inst,jcell,msg_tag in cell.comm_rules)]]))
                         for cell in thread_cells).read_frmt('x')
    cexinf.map_all_async(make_writing_message('initialize_cell_junctions', 'o',
                                      create_cell_junction_msg(cell))
                               for cell in thread_cells).read_frmt('x')

def create_cell_junction_msg(cell):
    surface_junctions = []
    line_junctions = []
    point_junctions = []
    index_axis = 'xyz'.index
    direction_map = {'+':1, '-':-1}
    for junction in cell.junctions:
        jcell_item = ['cell_index','i', cell.jcell_indices[junction.cell]]
        name = junction.df[0]
        if name=='surface':
            axis,direct = junction.df[1:]
            surface_junctions.append(NamedItems([
                         jcell_item,
                         ['axis', 'i', index_axis(axis)],
                         ['dir', 'i', direction_map[direct]]]))
        elif name=='line':
            axis,off1,off2 = junction.df[1:]
            line_junctions.append(NamedItems([
                         jcell_item,
                         ['axis', 'i', index_axis(axis)],
                         ['offset1', 'f', off1],
                         ['offset2', 'f', off2]]))
        elif name=='point':
            point_junctions.append(NamedItems([
                         jcell_item,
                         ['offset', 'v', junction.df[1:]]]))

        else:
            raise RuntimeError('bad name %r' % (name,))
    return NamedItems([
             ['jcells', 'o', StructArray(
                   [['comm_index', 'i', cell.jcell_indices[jcell]],
                    ['min_extent', 'v', jcell.extent.min_extent],
                    ['max_extent', 'v', jcell.extent.max_extent]]
                for jcell in cell.junctioned_cells)],
             ['surface_junctions', 'o', StructArray(surface_junctions)],
             ['line_junctions', 'o', StructArray(line_junctions)],
             ['point_junctions', 'o', StructArray(point_junctions)]])


