#!/usr/bin/python

from __future__ import print_function

SPE = 'Na'
DEG = 'I'
T = 10000
t_list = [1000, 2000, 4000, 6000, 8000, 10000,12000, 14000, 16000]

import pylab as pl
from mad import *
from mact_e import *

SPE = 'Na'
DEG = 'I'
T = 10000
t_list = [1000, 2000, 4000, 6000, 8000, 10000,12000, 14000, 16000]

lev1 = Level(e=0.0000, g=2, cfg='2p6.3s', term='2S', p='e')
state2 = State(e=16956.1702, g=2, cfg='2p6.3p', term='2P*', p='o')
state3 = State(e=16973.3661, g=4, cfg='2p6.3p', term='2P*', p='o')
lev2 = Level(state2, state3)

line = Line(lower=lev1, upper=lev2, f=1.)

SYM = str(SPE).lower() + str(DEG).lower()
ION = ions.get(SYM)

print('R0:', orbital_radius(ION, lev1))
print('R0:', orbital_radius(ION, lev2))
print('n*:', neff(ION, lev1))
print('n*:', neff(ION, lev2))


upse = ups_sce_vanregemorter(line.lower.g, line.f, 2.1, DEG)
ups_vr = [upse(t) for t in t_list]
print(t_list)
print('Van Regemorter:', ups_vr)

upse = ups_sce_seaton(ION, line)
ups_seaton = [upse(t) for t in t_list]
print('Seaton:        ', ups_seaton)
               
upse = ups_sce_fisher(line.lower.g, line.f, 2.1, line.dn())
ups_fisher = [upse(float(t)) for t in t_list]
print('Fisher:        ', ups_fisher)

print('x0:', 2.1 * Cst.Q / Cst.K / T)
print('Average gaunt Fisher:', average_gaunt_fisher(line.dn(), 2.1, T))
print('Average gaunt Van Regemorter:', average_gaunt_vanregemorter(DEG, 2.1, T))
print('Average gaunt Van Regemorter2:', average_gaunt_vanregemorter2(DEG, 2.1, T))
print('Average gaunt Seaton:', average_gaunt_seaton(ups_seaton[5], line))

x_vr = pl.linspace(0., 8., 100)
gaunt_vr = make_gaunt_vanregemorter('I')
gi_vr = gaunt_vr(x_vr)
gaunt_vr = make_gaunt_vanregemorter('II')
gii_vr = gaunt_vr(x_vr)

x_fisher = pl.linspace(0.,  8., 100)
gaunt_fisher = make_gaunt_fisher(line.dn())
gn0_fisher = gaunt_fisher(x_fisher)
gaunt_fisher = make_gaunt_fisher(1)
gn_fisher = gaunt_fisher(x_fisher)


pl.figure(figsize=(8, 6), dpi=100, facecolor='w', edgecolor='k')
pl.plot(x_vr, gi_vr, label='Van Regemorter (Neutral)')
pl.plot(x_vr, gii_vr, label='Van Regemorter (Positive ion)')
pl.plot(x_fisher, gn0_fisher, label='Fisher (Dn = 0)')
pl.plot(x_fisher, gn_fisher, label='Fisher (Dn > 0)')
pl.xlabel("sqrt(E/Eo)")
pl.ylabel("Gaunt factor")
pl.legend(loc=0, frameon=False)
pl.show()

pl.figure(figsize=(8, 6), dpi=100, facecolor='w', edgecolor='k')
pl.plot(t_list, ups_vr, label='Van Regemorter (1962)')
pl.plot(t_list, ups_fisher, label='Fisher (1994)')
pl.plot(t_list, ups_seaton, label='Seaton (1962)')
pl.xlabel("T [K]")
pl.ylabel("Maxwellian average Gaunt factor")
pl.legend(loc=0, frameon=False)
pl.show()


plot_ex(ION, line)


