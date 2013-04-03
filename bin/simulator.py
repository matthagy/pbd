#!/usr/bin/env python

'''simulate ensemble of colloid particles using Brown Dynamics simulation
'''

from __future__ import division

import sys
import os
from optparse import OptionParser
from math import ceil
from itertools import count

from hlab.pathutils import FilePath
from hlab import objstream

import autopath

from pbd.constants import nm, ns,mcs, ms
from pbd.debug import msg, error
from pbd.cex import CexInterface
from pbd.sim import Simulator
from pbd import state

def main():
    global config
    config = configure()
    parameters, configuration = load_initial_state()
    duration = config.duration
    if config.absolute_duration and duration:
        msg('offseting simulation duration of length %.3g ms by intial time %.3g ms',
            1e3 * duration, 1e3 * configuration.time)
        duration -= configuration.time
        if duration <= 0:
            msg('simulation complete')
            return
    cexinf = create_cex()
    sim = Simulator.create(cexinf,
                           parameters=parameters,
                           configuration=configuration,
                           random_seed=config.random_seed)
    outstream = initialize_output_stream(parameters, configuration)
    save_cycles = (duration and
                   int(ceil(duration / config.save_rate)))
    integration_cycles = int(ceil(config.save_rate / parameters.time_step))
    msg("simulating %s cycles of length %.3g mcs (%d integrations)",
        save_cycles or 'oo', config.save_rate / mcs, integration_cycles)
    for save_i in xrange(1,1+save_cycles) if save_cycles else count(1):
        sim.simulate(integration_cycles,
                     max_c_cycles=config.max_c_integrations)
        msg("saving cycle %d of %s", save_i, save_cycles or 'oo')
        outstream.write(sim.get_configuration()
                        if not config.thread_dump else
                        sim.get_state())
        outstream.flush()

def configure():
    parser = OptionParser(usage='%prog [OPTIONS] <filename>',
                          epilog=__doc__,
                          add_help_option=False)
    parser.add_option('-?','--help',
                      action='help',
                      help='show this message and exit')
    parser.add_option('--init-file',
                      dest='initfile',
                      action='store',
                      default=None,
                      metavar='FILE',
                      help='specify initial configuration file FILE; ' +
                           'otherwise load from simulation file')
    parser.add_option('--save-rate',
                      dest='save_rate',
                      action='store',
                      type='float',
                      default=10*mcs,
                      metavar='TIME',
                      help='specify how often to save configurations in seconds')
    parser.add_option('--duration',
                      dest='duration',
                      action='store',
                      type='float',
                      default=0,
                      metavar='TIME',
                      help='total time to simulate in seconds')
    parser.add_option('--absolute-duration',
                      dest='absolute_duration',
                      action='store_true',
                      default=False,
                      help='offset duration by time of initial configuration')
    parser.add_option('--max-c-integrations',
                      dest='max_c_integrations',
                      action='store',
                      type='int',
                      default=4000,
                      metavar='N',
                      help='perform at most N integrations in one C-call')
    parser.add_option('--clobber',
                      dest='clobber',
                      default=False,
                      action='store_true',
                      help='erase existing trajectory as opposed to appending')
    parser.add_option('--thread-dump',
                      dest='thread_dump',
                      default=False,
                      action='store_true',
                      help='dump state of individual threads instead of just positions;' +
                           'just for debugging')
    parser.add_option('--random-seed',
                      dest='random_seed',
                      default=None,
                      type='int',
                      action='store',
                      metavar='SEED',
                      help='specify random seed; otherwise randomly seeded from system entropy')
    parser.add_option('--nproc',
                      dest='nproc',
                      default=1,
                      type='int',
                      action='store',
                      metavar='N',
                      help='specify number of threads to use in simulation')
    parser.add_option('-m','--mpi',
                      dest='mpiargs',
                      default=[],
                      action='append',
                      metavar='ARGS',
                      help='pass argument onto mpirun')
    parser.add_option('--mpirun',
                      dest='mpirun',
                      default=None,
                      action='store',
                      metavar='PATH',
                      help='specify path to mpirun program')
    parser.add_option
    config,args = parser.parse_args()
    if not args:
        outfile = '-' #write to stdout by default
    elif len(args)>1:
        parser.error("takes exactly one argument; given %d" % (len(args),))
    else:
        outfile, = args
    config.outfile = outfile
    return config

def load_initial_state():
    global config
    if config.initfile is None:
        instream = create_objstream(config.outfile, mode='r', stdfileobj=None)
        if instream is None:
            error("no initial configuration file and output is pipe")
    else:
        instream = create_objstream(config.initfile, mode='r', stdfileobj=sys.stdin)
    instream.ignore_corrupt_entries = True
    try:
        parameters = instream[0]
        configuration = instream[-1]
    except EOFError:
        error("EOF on input file %s", instream.name)
    except IOError,e:
        error("error loading input file %s: %s", instream.name, e.strerror)
    if not isinstance(parameters, state.Parameters):
        error("in input file %s: bad first object %r; expected Parameters",
              instream.name, parameters)
    if not isinstance(configuration, state.Configuration):
        error("in input file %s: bad last object %r; expected Configuration",
              instream.name, config)
    return parameters, configuration

def initialize_output_stream(parameters, configuration):
    if config.clobber:
        obs =  create_objstream(config.outfile, "w", sys.stdout)
        #initialize new stream
        #due to issues with sympy's pickling, protocol 2 leads to errors
        #write parameters with basic protocol
        obs.pickle_protocol = 0
        obs.write(parameters)
        obs.pickle_protocol = 2
        obs.write(configuration)
        return obs
    #don't append to file unless we know the file is free of defects
    if config.initfile and os.path.realpath(config.outfile) != os.path.realpath(config.initfile):
        msg("warning, appending to different file than initialization file")
    if os.path.exists(config.outfile):
        obs = create_objstream(config.outfile, "r", None)
        if obs is None:
            error("can't append to standard stream")
        try:
            obs.read_all_locators()
        except objstream.CorruptFile, e:
            error("can't append to corrupt file %s; %s",
                  obs.name, e)
    return create_objstream(config.outfile, "a", None)

def create_objstream(s, mode='r', stdfileobj=None):
    assert mode in 'rwa'
    if s=='-':
        if mode=='a':
            error("can't append to standard stream")
        if stdfileobj is None:
            return None
        fileobj = stdfileobj
    else:
        try:
            fileobj = open(s, mode+'b')
        except IOError,e:
            error("couldn't open %s for %s; %s", s,
                  dict(r='reading', w='writing', a='appending')[mode],
                  e.strerror)
    stream = (objstream.Reader if mode=='r' else objstream.Writer)(fileobj)
    if mode in 'wa':
        stream.pickle_protocol = 2
    stream.name = fileobj.name
    return stream

def create_cex():
    if config.nproc == 1 and not config.mpiargs and config.mpirun is None:
        return CexInterface.create()
    return CexInterface.create_mpi(nproc=config.nproc,
                                   mpiargs=config.mpiargs,
                                   mpirun=config.mpirun)


__name__ == '__main__' and main()
