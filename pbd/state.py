##-*- Mode: python -*-
## state.py - Persistance State of Simulation
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

import numpy
from enthought.traits.api import (HasStrictTraits, Trait, Float,
                                  Range, List)

from .constants import R_particle, nm, mcm, mm, ps,ns
from .forcefield import PotentialBase, get_force_field


def validate_pair_potential(op, name, value):
    if isinstance(value, str):
        value = get_force_field(value)
    elif value is None:
        value = get_force_field("zero")
    elif not isinstance(value, PotentialBase):
        raise TypeError
    return value

validate_pair_potential.info = 'particle pairwise potential'

class Parameters(HasStrictTraits):
    '''configuration-independent description of simulation
    '''

    box_size = List(desc='''dimensions of cartesian space
                     units: meters
                     ''',
                    minlen=3, maxlen=3,
                    trait=Range(low=1*nm, high=1*mm),
                    value=list(mcm * numpy.ones(3)))

    temperature = Range(desc='''temperature of thermal bath
                        units: kelvin
                        ''',
                        low=1, high=2000, value=300)

    time_step = Range(desc='''brownian dynamics integration time step
                              units: seconds
                              ''',
                      low=1*ps, high=1000*ns, value=0.1*ns)

    dU_max = Range(desc='''maximum change in energy that can occur durring
                           any single particle integration.  subcycles are
                           used to ensure that each individual subcycle is
                           has a smaller change in energy than this limit
                           units: kT
                           ''',
                   low=1e-3, high=100, value=0.5)

    eta_solv = Range(desc='''viscosity of solvent in which particles are
                           emerged
                           units: pascal*second
                        ''',
                   low=1e-6, #dilute gas
                   high=1.0, #thick corn syrup
                   value=8.94e-4 #viscosity of water at 25C
                   )

    pair_potential = Trait(get_force_field("zero"),
                      validate_pair_potential,
                      desc='''particle pairwise potential
                      ''')

    r_potential_cutoff = Range(desc='''separation distance at which
                                       pair potential becomes zero
                                       units: meters
                                       ''',
                               low=2*R_particle, high=5*R_particle,
                               value=2*R_particle + 20*nm)

    force_update_rate = Range(desc='''reevaluate pairwise forces every
                                      n-integration cycles.  each force
                                      evaluation requires communicating
                                      the position of neighboring particles
                                      between threads
                                   ''',
                              low=1, high=10, value=1)

    linterp_size = Range(desc='''number of points to include in the linear
                                 interpolation table for pairwsie for potential
                                 and force evaluation durring simulation.  larger
                                 tables more accuratley represent details of functions
                                 while using more memory (and potentially smashing
                                 the L2 cache)
                                 ''',
                         low=2, #used for flat potentials w/ ideal gas
                         high=10000,
                         value=250)

    linterp_r_min = Range(desc='''minimum separation distance to model in
                                  linear interpolation tables.  larger
                                  values improve table resolution (higher
                                  point density), but simulation behavior
                                  is not defined should two particles come
                                  closer than this distance (likely SEGFAULT)
                                  units: meters
                                  ''',
                          low=0, high=2*R_particle, value=0)


    r_neighbor = Range(desc='''all pairs of particles within this
                               distance are included in the pairwise
                               particle neighbor lists.  along with
                               r_potential_cutoff, controls lifetime
                               of neighbor lists.
                               units: meter
                            ''',
                       low=2*R_particle, high=10*R_particle,
                       value=2*R_particle + 30*nm)



def validate_positions(op, name, value):
    value = numpy.asarray(value, float)
    if not (len(value.shape)==2 and value.shape[1]==3):
        raise TypeError
    return value

validate_positions.info = 'a float array of shape (N,3) where N is the number of particles'

class Configuration(HasStrictTraits):
    '''single snapshot of particles
    '''

    time = Range(desc='''time associated with this configuration, relative
                         to some initial reference configuration.
                         units: seconds
                      ''',
                 low=0.0, high=None, value=0.0)

    positions = Trait(None, validate_positions,
                      desc='''location of all particles in configuration
                      ''')

    wall_time = Range(desc='''time at which state was generated
                      deltas between states can be used for rough
                      bench marking
                      ''',
                      low=0.0, high=None, value=0.0)


def validate_tags(op, name, value):
    value = numpy.asarray(value, int)
    if len(value.shape)!=1:
        raise TypeError
    return value

validate_tags.info = 'a 1 dimensional integer array specifying particle tags'

def validate_neighbors(op, name, value):
    value = numpy.asarray(value, int)
    if not (len(value.shape)==2 and value.shape[1]==2):
        raise TypeError
    return value

validate_neighbors.info = 'an integer array of shape (N,2) where N is number of neighbor pairs'

class ThreadState(HasStrictTraits):
    '''Snapshot of internal state in a thread
       Mainly used for debugging
    '''

    positions = Trait(None, validate_positions,
                      desc='''location of all particles in thread
                              including both internal and external particles
                      ''')

    tags = Trait(None, validate_tags,
                 desc='''unique tags of internal particles in thread
                 ''')

    internal_neighbors = Trait(None, validate_neighbors,
                               desc='''internal neighbor table
                               ''')

    external_neighbors = Trait(None, validate_neighbors,
                               desc='''external neighbor table
                               ''')


class SimulationState(HasStrictTraits):
    '''State of all threads
    '''

    time = Range(desc='''time associated with this configuration, relative
                         to some initial reference configuration.
                         units: seconds
                      ''',
                 low=0.0, high=None, value=0.0)

    threads = List(trait=Trait(ThreadState),
                   desc='''all threads in simulation
                   ''')

    wall_time = Range(desc='''time at which state was generated
                      deltas between states can be used for rough
                      bench marking
                      ''',
                      low=0.0, high=None, value=0.0)
