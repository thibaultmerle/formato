#!/usr/bin/python
# -*- coding: utf-8 -*-
''' FoRMATo Tool to create model atom for NLTE radiative transfer in late-type stars.
    Need as input binary output files of:
    mael.py    => 
    mart_bb.py => 
    mart_bf.py
    macte.py
    macth.py
'''

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
#from __future__ import unicode_literals

import argparse as ap
import pickle
import time
import sys
import os

from mad import Cst

def read_input(ifilename):
    ''' Read input binary file
    '''
    ifile = open(ifilename, 'rb')
    data = pickle.load(ifile)
    ifile.close()
    return data    

IFILE1 = 'ael'
IFILE2 = 'art_bb'
IFILE3 = 'act_e'
IFILE4 = 'act_h'
OFILE = 'atom'

DESCRIPTION = 'FoRMATo Tool to create model atom for NLTE radiative transfer in late-type stars.'
EPILOG = '2014-03-26 ThiM'

parser = ap.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
parser.add_argument('species', type=str, default=None, help='Symbol of the atomic species (e.g. "Fe" or "FE" or "fe")')
parser.add_argument('degree', nargs ='*', type=str, default=None, help='Ionization degrees (e.g. "I", "II" or "I II" for neutral, ionized or combined species)')

arguments = parser.parse_args()

SPE = arguments.species
DEG = arguments.degree

SYM_LIST = []
for deg in DEG:
    SYM_LIST.append(str(SPE).lower() + str(deg).lower())

nion = len(SYM_LIST)

print('__________________')
print('|   formato.py   |')
print('TTTTTTTTTTTTTTTTTT')

# Check the existence of mass of considered element
try: 
    mass = Cst.MASS[SPE.upper()]
except KeyError:
    print('Add the mass of this element in mad.Cst.MASS dictionnary.')
    quit(1)

IFILE1_LIST = []
IFILE2_LIST = []
IFILE3_LIST = []
IFILE4_LIST = []

for SYM in SYM_LIST:
    IFILE11 = IFILE1+'_'+SYM+'.bin'
    IFILE22 = IFILE2+'_'+SYM+'.bin'
    IFILE33 = IFILE3+'_'+SYM+'.bin'
    IFILE44 = IFILE4+'_'+SYM+'.bin'
    if os.path.isfile(IFILE11):
        IFILE1_LIST.append(IFILE11) 
    else:
        print("File", IFILE11, "does NOT exist.")
        quit(1)
    if os.path.isfile(IFILE22):
        IFILE2_LIST.append(IFILE22)
    else:
        print("File", IFILE22, "does NOT exist.")
        quit(1)
    if os.path.isfile(IFILE33):
        IFILE3_LIST.append(IFILE33)
    else: 
        print("File", IFILE33, "does NOT exist.")
        quit(1)
    if os.path.isfile(IFILE44):
        IFILE4_LIST.append(IFILE44)
    else: 
        print("File", IFILE44, "does NOT exist.")
        quit(1)    

# Read levels
levels1 = read_input(IFILE1_LIST[0])
nlev1 = len(levels1)
# Read radiative b-b transitions
lines1 = read_input(IFILE2_LIST[0])
# Read electron collision transitions
ecols1 = read_input(IFILE3_LIST[0])
# Read hydrogen collision transition
hcols1 = read_input(IFILE4_LIST[0])

if nion == 2:
    levels2 = read_input(IFILE1_LIST[1])
    nlev2 = len(levels2)
    lines2 = read_input(IFILE2_LIST[1])
    ecols2 = read_input(IFILE3_LIST[1])
    hcols2 = read_input(IFILE4_LIST[1])

if nion > 0:
    print("\nIon considered:", SPE, DEG[0])
    print(" Number of levels (including ionization stage):", str(nlev1).rjust(6))
    print(" Number of radiative bound-bound transitions:  ", str(len(lines1)).rjust(6))
    print(" Number of e⁻ collision transitions:           ", str(len(ecols1)).rjust(6))
    print(" Number of H collision transitions:            ", str(len(hcols1)).rjust(6))
if nion == 2:
    print("\nIon considered:", SPE, DEG[1])
    print(" Number of levels (including ionization stage):", str(nlev2).rjust(6))
    print(" Number of radiative bound-bound transitions:  ", str(len(lines2)).rjust(6))
    print(" Number of H collision transitions:            ", str(len(hcols2)).rjust(6))

#Remove ionization stage of first species to avoid repetition
# and combine levels1 and levels2 if necessary
if nion == 2:
    ION = levels1.pop(len(levels1)-1)
    for level in levels2:
        level.e += ION.e
        level.__format__()
    levels = levels1 + levels2
else:
    levels = levels1

print('\nGround stage:            ', end='')
levels[0].print()
print('Ionization level:        ', end='')
if nion == 2:
    ION.print()
else: 
    levels[-1].print()
if nion == 2:
    print('Second ionization level: ', end='')
    levels[-1].print()

#Combine lines1 and lines2
if nion == 2:
    for line in lines2:
        line.l += len(levels1)
        line.u += len(levels1)
    lines = lines1 + lines2
else:
    lines = lines1

#Combine ecols1 and ecols2
if nion == 2:
    for ecol in ecols2:
        ecol.l += len(levels1)
        ecol.u += len(levels1) 
    ecols = ecols1 + ecols2
else:
    ecols = ecols1

#Combine hcols1 and hcols2
if nion == 2:
    for hcol in hcols2:
        hcol.l += len(levels1)
        hcol.u += len(levels1)
    hcols = hcols1 + hcols2
else:
    hcols = hcols1

OFILE += '.'+SPE.lower()
header = SPE.upper()

for deg in DEG:
    if nion != 1:
        OFILE += '-'
    header += ' '
    OFILE += deg.lower()
    header += deg.upper()

header += "              generated by FoRMATo-2.0, " + time.asctime() 

stdout_save = sys.stdout
ofile = open(OFILE, 'w')
sys.stdout = ofile

print(header)
print('* ABUND   AWGT')
print('  X.XX  ', mass)
print('* NK NLIN NCNT NFIX')
header = str(len(levels)) +' '+ str(len(lines)) +' '+ str(len(levels)-1) +' 0'
print(header)

print('*    E[cm⁻¹]    G  CONFIGURATION & TERM                          ION NK')
for idx, level in enumerate(levels):
    if type(level.cfg) == str:
        label = (level.cfg).strip()
    elif type(level.cfg) == list:
        label = (' '.join(level.cfg)).strip()
    if type(level.term) == str:
        label += ' '+level.term+'                   '
    elif type(level.term) == list:
        label += ' '.join(level.term)+'                   '
    label = label[:20]
    towrite = level.e_p+level.g_p+" '"+label+"' "+str(level.i)+str(idx+1).rjust(4)
    print(towrite)

print('* RADIATIVE B-B TRANSITIONS')
print('*          F       NQ QMAX Q0 IW    GA       GVW      GS      LAMBDA[Å]    KR')
NQ = ' 30'
QMAX = '  50'
Q0 = ' 0.3'
IW = ' 0 '
for idx, line in enumerate(lines):
    towrite = str(line.u).rjust(4) + str(line.l).rjust(4)+ line.f_p + NQ + QMAX + Q0 + IW + line.gr_p + line.gv_p + line.gs_p + line.lbd_p + str(idx+1).rjust(6)
    print(towrite)

print('* RADIATIVE B-F TRANSITIONS')
print('* ION   I    A0   NPTS')

NQ = 10
LBD_MIN = 911.0

for idx, level in enumerate(levels[:nlev1-1]):
    if hasattr(level, 'pe') and hasattr(level, 'px'):
        print(str(nlev1).rjust(4), str(idx+1).rjust(4), format(level.px[0], '5.2e'), len(level.px), -1)
        for e, x in zip(level.pe, level.px):
            print(format(e, '8.2f'), format(x, '5.2e'))
    else:
        print(str(nlev1).rjust(4), str(idx+1).rjust(4), format(level.sbf,'5.2e'), NQ, LBD_MIN)
if nion == 2:
    for idx, level in enumerate(levels[nlev1: nlev2]):
        if hasattr(level, 'pe') and hasattr(level, 'px'):
            print(str(nlev2).rjust(4), str(idx+nlev1).rjust(4), format(level.px[0], '5.2e'), len(level.px), -1)
            for e, x in zip(level.pe, level.px):
                print(format(e, '8.2f'), format(x, '5.2e'))
        else:        
            print(str(nlev2).rjust(4), str(idx+nlev1).rjust(4), format(level.sbf,'5.2e'), NQ, LBD_MIN)

print('* COLLISION TRANSITIONS')
print('GENCOL')
print('* ELECTRONIC COLLISIONS')
print('TEMP')

temp = ''
for i in range(ecols[0].nt):
    temp += ' '+str(ecols[0].t_list[i])

print(ecols[0].nt, temp)
ecols[0].print(tot=True)

for idx, ecol in enumerate(ecols[1:]):
    if ecol.nt != ecols[idx].nt:
        temp = ''
        for i in range(ecol.nt):
            temp += ' '+str(ecol.t_list[i])
        if ecol.kw.strip() != 'CI':
            print('TEMP')
            print(ecol.nt, temp)
    ecol.print(tot=True)

print('* HYDROGEN COLLISIONS')
print('TEMP')

temp = ''
for i in range(hcols[0].nt):
    temp += ' '+str(hcols[0].t_list[i])

print(hcols[0].nt, temp)
hcols[0].print(tot=True)

for idx, hcol in enumerate(hcols[1:]):
    if hcol.nt != hcols[idx].nt:
        temp = ''
        for i in range(hcol.nt):
            temp += ' '+str(hcol.t_list[i])
#        if hcol.kw.strip() != 'CI':
        print('TEMP')
        print(hcol.nt, temp)
    hcol.print(tot=True)

print('END')
ofile.close()

sys.stdout = stdout_save

print('\nModel atom of ', SPE, DEG)
print(" Number of levels (including ionization stage):", str(len(levels)).rjust(6))
print(" Number of radiative bound-bound transitions:  ", str(len(lines)).rjust(6))
print(" Number of e⁻ collision transitions:           ", str(len(ecols)).rjust(6))
print(" Number of H collision transitions:            ", str(len(hcols)).rjust(6))

print('\nOutput ascii file:', OFILE)


