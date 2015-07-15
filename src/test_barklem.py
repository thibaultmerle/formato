#!/usr/bin/python

import pylab as pl
from barklem import retcross

#Spec     Wavel      Elow     Elimit       Eupp     Elimit   L      J     N*low N*upp  Sigma    Alpha Log(Gamma/N)
#  
# Multiplet 1 of Na                                                                                  
#Na 1  5889.953      0.000  41449.650  16973.379  41449.650 0->1 0.5->1.5 1.627 2.117  406.87   0.273  -7.526

NSlow, NSupp = 1.627, 2.117
Llow, Lupp = 0, 1

# First test with neff

CROSS, ALPHA, IFAIL = retcross(NSlow,NSupp,Llow,Lupp)
print CROSS, ALPHA, IFAIL

# Second test with energy levels

Elow = 0.000
Eupp = 16973.379
Elimit = 41449.650
Zeff = 1 

RYD = 109678.758
#RYD = 109737.32
NSlow = Zeff * pl.sqrt(RYD/(Elimit-Elow))
NSupp = Zeff * pl.sqrt(RYD/(Elimit-Eupp))

print NSlow, NSupp
CROSS, ALPHA, IFAIL = retcross(NSlow,NSupp,Llow,Lupp)
print CROSS, ALPHA, IFAIL
