#!/usr/bin/python

import pylab as pl
import janicki as j

E = pl.arange(0.01, 1000, 0.01)
Z = 1

def f(N, Z, E):
	gbf_av = []
	for x in E:
		dumb, av = j.gbf_ave(N, Z, x)
		gbf_av.append(av)
	return gbf_av

N = 1
gbf_av1 = f(N, Z, E)
N = 2
gbf_av2 = f(N, Z, E)
N = 3
gbf_av3 = f(N, Z, E)
N = 4
gbf_av4 = f(N, Z, E)
N = 5
gbf_av5 = f(N, Z, E)
N = 6
gbf_av6 = f(N, Z, E)
N = 7
gbf_av7 = f(N, Z, E)



pl.figure(facecolor='w', edgecolor='k')
pl.xscale('log')
E = [e/13.6057/float(Z)**2 for e in E]
pl.plot(E, gbf_av1, label="n=1")
pl.plot(E, gbf_av2, label="n=2")
pl.plot(E, gbf_av3, label="n=3")
pl.plot(E, gbf_av4, label="n=4")
pl.plot(E, gbf_av5, label="n=5")
pl.plot(E, gbf_av6, label="n=6")
pl.plot(E, gbf_av7, label="n=7")
pl.xlabel("$E/Z^2$")
pl.ylabel("$\\bar{g_{bf}}$")
pl.legend(loc='best')
pl.show() 
