##-*- Mode: python -*-
## forcefield.py - Defines pair-potentials for forcefield
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

from __future__ import division

import sympy as S
import numpy as N
from scipy_optimize import fminbound

from jamenson.runtime.multimethod import MultiMethod, defmethod

from hlab.ihlab import IDomain
from hlab.infinity import oo
from hlab.domain import ds, Boundary
from hlab.calculated import PickelingBase, calculated
from hlab.prof import Profile
from hlab.xalgebra import (AlgebraBase, mm_neg, defboth_mm_eq,
                           defboth_mm_add, defboth_mm_mul,
                           scalar_number_type)
from hlab.xlinterp import HomogenousTable
from hlab.piecewise import PieceWiseCollection, PieceWiseFunction

try:
    from Gnuplot import Gnuplot, Data
except ImportError:
    have_gnuplot = False
else:
    have_gnuplot = True

from constants import R_particle, nm, kB

r = S.Symbol('r')
h = r - 2*R_particle

def make_r_evaluator(expr):
    '''from a sympy expression, where r is the only free variable,
       create a function of r
    '''
    exprs = str(expr)
    f =  eval('lambda r: %s' % (exprs,), vars(N))
    f.func_name = 'eval: %s' % (exprs,)
    return f

class PotentialBase(PickelingBase, AlgebraBase):

    def __init__(self, name):
        self.name = name

    def __str__(self):
        return '<%s %s>' % (self.__class__.__name__, self.name)

    def make_potential_table(self, r_min, r_max, n_points):
        return HomogenousTable.fromfunc(self.potential_eval, r_min, r_max, n_points)

    def make_force_table(self, r_min, r_max, n_points):
        return HomogenousTable.fromfunc(self.force_eval, r_min, r_max, n_points)

    if have_gnuplot:
        def make_potential_plot(self, **kwds):
            return self.make_plot(self.potential_eval, _base_="potential", **kwds)

        def make_force_plot(self, **kwds):
            return self.make_plot(self.force_eval, _base_="force", **kwds)

        def plot_potential(self, **kwds):
            return plot_potentials(self, **kwds)

        def plot_force(self, **kwds):
            return plot_forces(self, **kwds)

        def make_plot(self, func, _base_='func', with_='l', rs=None, title=None):
            if rs is None:
                rs = N.linspace(135*nm, 155*nm, 100)
            y = map(func, rs)
            if title is None:
                title = '%s:%s' % (self.name, _base_)
            return Data(N.asarray(rs)/nm, y, title=title, with_=with_)


if have_gnuplot:
    def do_plots(func, *pots, **kwds):
        try:
            persist = kwds.pop('kwds')
        except KeyError:
            persist = True
        gp = Gnuplot(persist=persist)
        gp.xlabel('Separation Distance (nm)')
        gp.plot(*list(func(pot, **kwds) for pot in pots))
        return gp
    def plot_potentials(*pots, **kwds):
        gp = do_plots(lambda p, **kwds: p.make_potential_plot(**kwds), *pots, **kwds)
        gp.ylabel('Energy (kT)')
        gp.replot()
        return gp
    def plot_forces(*pots, **kwds):
        gp = do_plots(lambda p, **kwds: p.make_force_plot(**kwds), *pots, **kwds)
        gp.ylabel('Force (kT/nm)')
        gp.replot()
        return gp


class AnalyticPotential(PotentialBase):
    '''define in terms of a sympy expression
    '''

    def __init__(self, name, potential_expr):
        self.name = name
        self.potential_expr = potential_expr

    @calculated
    def force_expr(self):
        return -self.potential_expr.diff(r)

    @calculated
    def potential_eval(self):
        return make_r_evaluator(self.potential_expr)

    @calculated
    def force_eval(self):
        return make_r_evaluator(self.force_expr)

@defboth_mm_add([AnalyticPotential, AnalyticPotential])
def meth(a,b):
    return AnalyticPotential('%s+%s' (a.name,b.name),
                             a.potential_expr + b.potential_expr)

@defboth_mm_add([AnalyticPotential, scalar_number_type])
def meth(a,s):
    return AnalyticPotential('%s+%s' (a.name,s),
                             a.potential_expr + s)

@defboth_mm_mul([AnalyticPotential, AnalyticPotential])
def meth(a,b):
    return AnalyticPotential('%s*%s' (a.name,b.name),
                             a.potential_expr * b.potential_expr)

@defboth_mm_mul([AnalyticPotential, scalar_number_type])
def meth(a,s):
    return AnalyticPotential('%s*%s' (a.name+s),
                             a.potential_expr * s)


class InterpolatedPotential(PotentialBase):
    '''defined in terms of a HomogenousTable
    '''

    def __init__(self, name, potential_linterp):
        assert isinstance(potential_linterp, HomogenousTable)
        self.name = name
        self.potential_linterp = potential_linterp

    @calculated
    def force_linterp(self):
        return -self.potential_linterp.diff()

    @calculated
    def potential_eval(self):
        return self.potential_linterp.interpolate

    @calculated
    def force_eval(self):
        return self.force_linterp.interpolate


class PieceWisePotential(PotentialBase):
    '''defined as a collection of other potentials
    '''

    def __init__(self, name, pieces=None):
        self.name = name
        self.pieces = PieceWiseCollection(pieces)

    def add_piece(self, domain, piece):
        self.pieces.add(domain, piece)

    @calculated
    def potential_eval(self):
        return PieceWiseFunction(self.pieces.map(lambda piece: piece.potential_eval))

    @calculated
    def force_eval(self):
        return PieceWiseFunction(self.pieces.map(lambda piece: piece.force_eval))


force_field_cache = {}

def get_force_field(name):
    try:
        return force_field_cache[name]
    except KeyError:
        ff = force_field_cache[name] = make_force_field(name)
        return ff

make_force_field = MultiMethod('make_force_field')

def closure():
    '''original potentials
    '''
    u0 = 20
    u1 = 10
    l0 = 3*nm
    l1 = 0.5*nm
    l2 = 1*nm

    @defmethod(make_force_field, ["zero"])
    def meth(name):
        return InterpolatedPotential(name, HomogenousTable.fromfunc(lambda r: 0, -1e10, 1e10, 20))

    @defmethod(make_force_field, ["repulsive"])
    def meth(name):
        U_rep = u0*S.exp(-h/l0) + u1*S.exp(-h/l1)
        return AnalyticPotential(name, U_rep)

    @defmethod(make_force_field, ["attr-low-barrier"])
    def meth(name):
        U_lb = (u1 + u0*S.exp(-h/l0)) / (1 + u0*S.exp(-h/l2)) - u1 * (1 - S.exp(-h/l1))
        return AnalyticPotential(name, U_lb)

    @defmethod(make_force_field, ["attr-high-barrier"])
    def meth(name):
        return AnalyticPotential(name, 5 * get_force_field('attr-low-barrier').potential_expr)

closure()
del closure

def closure():
    '''kostansek potentials
    '''
    # potential from kostansek1996
    A = 1e-13 * 1e-7 #1e-13 ergs (1erg = 1pe-7 J)
    kappa = 1 / (2*nm) #2-5 nm
    F = 5e-11 # 5e-11 - 200e-11 N
    L = 1*nm #posibly larger
    epsilon_r = 80 #relative permitivity of H20 at ~300K
    epsilon_0 = 8.85e-12 #permitivity of vaccum (A^2s^4/Kg/m^3)
    psi = -55e-3 #55mV zeta potential

    @defmethod(make_force_field, ["kostansek"])
    def meth(name):
        U_attr = -A*R_particle/(12*h)
        U_rep = 2*R_particle*S.pi*epsilon_r*epsilon_0*psi**2*S.log(1+S.exp(-kappa*h))
        U_hyd = -F*R_particle*S.exp(-h/L)
        U_tot = U_attr + U_rep + U_hyd
        kT = kB*298
        U_tot = U_tot / kT
        return AnalyticPotential(name, U_tot)

    @defmethod(make_force_field, ["kostansek-shift"])
    def meth(name):
        U_tot = get_force_field("kostansek").potential_expr
        U_shift = U_tot - U_tot.subs(r, 155*nm)
        return AnalyticPotential(name, U_shift)

    def fudge_kostansek(name, hc):
        rc = 2*R_particle + hc * nm
        kos = get_force_field("kostansek-shift")
        U_tot = kos.potential_expr
        U_k = U_tot
        dU_k = U_k.diff(r)
        A,B,C,U,dU,R = map(S.Symbol, 'A B C U dU R'.split())
        diam = 2*R
        sol = S.solve_poly_system([
              A*diam**2 +  B*diam + C,
              A*r**2 +  B*r + C - U,
              2*A*r +   B       - dU],
            A, B, C)
        sol, = sol
        A,B,C = [S.simplify(expr.subs(U,U_k).subs(dU, dU_k).subs(r,rc).subs(R,R_particle))
                 for expr in sol]
        para = A*r**2 + B*r + C
        para = S.simplify(para)
        para = AnalyticPotential(name + '-para', para)
        return PieceWisePotential(name, [[ds[-oo:rc], para],
                                         [ds[Boundary(rc, inclusive=False):oo], kos]])

    @defmethod(make_force_field, ["kostansek-cut3.4"])
    def meth(name):
        return fudge_kostansek(name, 3.4)
    @defmethod(make_force_field, ["kostansek-cut3.0"])
    def meth(name):
        return fudge_kostansek(name, 3.0)

closure()
del closure

def closure():
    #from decimal import Decimal as dc
    dc = float
    n = dc("36")
    k = dc("0.1")
    Phi_p = dc("0.4")
    Phi = dc("0.3")
    L = 1 + k
    A = dc("2") / dc("3") * L ** dc("3")
    Q = dc("3") * Phi_p / (2 * k ** dc("3"))

    def U_HC(r_d):
        return r_d ** -n
    def H(r_d):
        return 1 if r_d < L else 0
    def U_AO(r_d):
        return Q * H(r_d) * (r_d * L**dc("2") -
                      (r_d ** dc("3")) / dc("3") -
                      A)
    d = 2 * R_particle
    def U(r):
        r_d = dc(repr(r / d))
        if r_d < 0.05:
            r_d = 0.05
        return U_HC(r_d) + U_AO(r_d)
    tbl = HomogenousTable.fromfunc(U, -10*nm, 200*nm, 10000)
    @defmethod(make_force_field, ["asakura-oosawa"])
    def meth(name):
        return InterpolatedPotential(name, tbl)

closure()
del closure

def get_scaled_barrier_kostansek(scale, base=None):
    if base is None:
        base = get_force_field('kostansek-cut3.4')
    U_eval = base.potential_eval
    r_max = fminbound(lambda r: -U_eval(r), 2*R_particle, 2*R_particle + 20*nm, xtol=0.01*nm)
    r_min = fminbound(U_eval, 2*R_particle, r_max, xtol=0.01*nm)
    r_0 = fminbound(lambda r: abs(U_eval(r)), 2*R_particle, r_max, xtol=0.01*nm)

    r_basis = N.linspace(0, 2*R_particle + 30*nm, 2000)
    r_sig = r_basis[N.where(r_basis>r_0)]
    r_flip = r_0 + 0.05 * (r_max - r_0)
    kappa = 100*1e10
    trans_sig = N.ones_like(r_sig) - (1-scale) * 1.0 / (1.0 + N.exp(-kappa*(r_sig-r_flip)))
    transformation = N.concatenate([N.ones_like(r_basis[N.where(r_basis<=r_0)]),
                                    trans_sig])
    U_basis = N.array(map(U_eval, r_basis))
    U = U_basis * transformation
    U = Profile(U).smoothed_guass(window_size=5, falloff=0.01).data
    return InterpolatedPotential(base.name + 'scaled%.2f' % (scale,),
                                 HomogenousTable(U, r_basis[0],
                                                 r_basis[1] - r_basis[0]))

def test():
    gp = plot_forces(#get_force_field('repulsive'),
                     #    get_force_field('attr-low-barrier'),
                     #    get_force_field('attr-high-barrier'),
                     #    get_force_field('kostansek'),
                         get_force_field('kostansek-cut3.4'),
                     #    get_force_field('kostansek-cut3.0'),
                         get_scaled_barrier_kostansek(1.0),
                         get_scaled_barrier_kostansek(0.8),
                         get_scaled_barrier_kostansek(0.5),
                         get_scaled_barrier_kostansek(0.1)
                         )
    gp('set yrange[-3e10:1e10]')
    #gp('set yrange[-20:20]')
    gp('set xrange[134:155]')
    gp.replot()
    from time import sleep
    sleep(1)


#__name__ == '__main__' and test()
