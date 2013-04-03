'''Routines for dividing cartessian space into rectangular cells
   and determining the junctions between these cells
'''

from __future__ import division
from __future__ import with_statement
from __future__ import absolute_import

from numpy import *
from collections import defaultdict

from jamenson.runtime.atypes import typep, as_optimized_type, anytype
from jamenson.runtime.atypes.ptypes import integer_type
from jamenson.runtime.multimethod import MultiMethod, defmethod

class Extent(object):

    def __init__(self, min_extent, max_extent):
        self.min_extent = min_extent
        self.max_extent = max_extent


class Cell(object):

    def __init__(self, index, extent, positions):
        self.index = index
        self.extent = extent
        self.positions = positions
        self.junctions = []

    def __hash__(self):
        return hash(tuple(self.index))

    def __eq__(self, other):
        return self is other


def divs():
    """We divide cartessian space by cutting each axis into segments and
       different divisions to get a specific number of cells.  For instance,
       the division (2,1,3) divides the x-axis in half, dosn't divide the y-axis,
       and divides the z-axis into thirds.  This gives a total of 6 cells, all of
       uniform dimensions.

       Here we calculate the ideal divisions to get a specific number of cells.  We
       define ideal as giving cells with most uniform dimensions, i.e. more
       squarish and lest rectangularish.
    """
    MAX_DIMENSION_DIVIDES = 10
    m = defaultdict(list)
    xr = 1+arange(MAX_DIMENSION_DIVIDES)
    for x in xr:
        for y in xr:
            for z in xr:
                div = array([x,y,z])
                m[multiply.reduce(div)].append(div)
    return dict((size,sorted(divs, key=lambda divs: (divs.std(), tuple(divs)))[0])
                for size,divs in m.iteritems())

divs = divs()

def get_divs(size):
    return divs[size]

# Calculate the offset that describe the junction data
# structures in md.c.  See the section on junctions in
# this source file for an explaination.
offsets = [
    ['surface', [
     [axis,i]
      for axis in 'xyz' for i in [-1,1]]],
    ['line', [
     (x,y,z)
     for x in [-1,0,1] for y in [-1,0,1] for z in [-1,0,1]
     if sum(abs(array([x,y,z])))==2
     ]],
    ['point', [
     [x,y,z]
     for x in [-1,1] for y in [-1,1] for z in [-1,1]]]]

calculate_junction_offset = MultiMethod('calculate_junction_offset',
                                        signature='name,cell,box_size,data')

def offset_shift(off, cell):
    one = ones_like(off)
    return ((cell.extent.max_extent - cell.extent.min_extent) *
            (one + (off - one) * 0.5) +
            cell.extent.min_extent)

@defmethod(calculate_junction_offset, '"surface",Cell,ndarray,anytype')
def meth(name, cell, box_size, (axis, axis_off)):
    offset = array(dict(x=[1,0,0], y=[0,1,0], z=[0,0,1])[axis]) * axis_off
    return offset, (axis, '-' if axis_off==-1 else '+')

@defmethod(calculate_junction_offset, '"line",Cell,ndarray,anytype')
def meth(name, cell, box_size, (x,y,z)):
    v = array([x,y,z])
    axis = 'xyz'[where(v==0)[0]]
    off1,off2 = [i for i in v if i!=0]
    if axis=='x':
        offset = [0,off1,off2]
    elif axis=='y':
        offset = [off1,0,off2]
    elif axis=='z':
        offset = [off1,off2,0]
    else:
        raise 'hell'
    offset =  array(offset)
    return offset, (axis,) + tuple(offset_shift(offset, cell)[where(offset!=0)])

@defmethod(calculate_junction_offset, '"point",Cell,ndarray,anytype')
def meth(name, cell, box_size, offset):
    offset = array(offset)
    return offset, tuple(offset_shift(offset, cell))


class Junction(object):

    def __init__(self, cell, df):
        self.cell = cell
        self.df = df

def partition_positions(size, positions, divs):
    if typep(divs, integer_type):
        divs = get_divs(divs)
    divs = asarray(divs, int)
    indices = floor(((positions / size) * divs)).astype(int)
    cells = dict(((x,y,z),Cell(array([x,y,z]), None, []))
                 for x in xrange(divs[0])
                 for y in xrange(divs[1])
                 for z in xrange(divs[2]))
    for index,position in zip(map(tuple, indices), positions):
        cells[index].positions.append(position)
    adj = ones(3, int)
    for cell in cells.itervalues():
        cell.positions = array(cell.positions).reshape(len(cell.positions), 3)
        min_extent = (cell.index / divs) * size
        max_extent = ((cell.index + adj) / divs) * size
        cell.extent = Extent(min_extent, max_extent)
        for name,defs in offsets:
            for df in defs:
                off, junction_df = calculate_junction_offset(name, cell, size, df)
                junction_cell = cells[tuple((cell.index + off) % divs)]
                if junction_cell is not cell:
                    cell.junctions.append(Junction(junction_cell, (name,) + tuple(junction_df)))
    return sorted(cells.itervalues(), key=lambda cell: tuple(cell.index))

# Calculation of communication rules between cells
# TODO: Explain this

class Link(object):

    def __init__(self, cell_i, cell_j):
        key_i = tuple(cell_i.index)
        key_j = tuple(cell_j.index)
        assert key_i != key_j
        if key_i > key_j:
            key_i, key_j = key_j, key_i
            cell_i, cell_j = cell_j, cell_i
        self.key = key_i, key_j
        self.cell_i = cell_i
        self.cell_j = cell_j

    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        if type(self) is not type(other):
            return NotImplemented
        return self.key == other.key


def group_junctions(cells):
    precedence_map = {None: 0,
                      'point': -1,
                      'line': -2,
                      'surface': -3}
    comms = defaultdict(int)
    for cell in cells:
        for junction in cell.junctions:
            link = Link(cell, junction.cell)
            prec = precedence_map[junction.df[0]]
            comms[link] = min(prec, comms[link])
    groups = defaultdict(list)
    for link,prec in comms.iteritems():
        groups[prec].append(link)
    if not groups:
        return []
    _,groups = zip(*sorted(groups.iteritems()))
    return groups

def calculate_communication_rules(cells):
    groups = group_junctions(cells)
    rounds = []
    round_comms = []
    comms_in_round = set()
    #XXX group rounds by type of junctions
    #perform surface, then line, then point communications
    for comms in groups:
        while comms:
            for link in comms:
                if not (link.cell_i in comms_in_round or
                        link.cell_j in comms_in_round):
                    break
            else:
                #start another round, everyone is busy
                rounds.append(round_comms)
                round_comms = []
                comms_in_round.clear()
                continue #while
            comms.remove(link)
            comms_in_round.add(link.cell_i)
            comms_in_round.add(link.cell_j)
            round_comms.append(link)
    if round_comms:
        rounds.append(round_comms)
    rounds.sort(key=len, reverse=True)
    return rounds


