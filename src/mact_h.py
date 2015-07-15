#!/usr/bin/python
# -*- coding: utf8 -*-
''' Code to merge atomic hydrogen collisions of a given species in a given ionization degree.
'''

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import argparse as ap
import operator as op
import pickle
import sys
import os

try: 
    import pylab as pl 
except ImportError: 
    raise ImportError('pylab module not installed')

from mad import Cst, State, Collision

iplot = 0 # For numeroation of cross-section plots

def ups_sch_drawin(gi, f, Eo, mass):
    ''' Compute the effective hydrogen collision strengths
        with semi-classical formula from Drawin 1969 Zeitschrift fur Physik 225, 483
        Eo [eV] transtion energy
        mass [amu] atomic weigth in atomic mass unit
    '''
    me = Cst.ME
    mh = Cst.MASS['H']*Cst.AMU
    ma = mass*Cst.AMU

    cst = 4*pl.sqrt(2)* me * ma / mh**2 * Cst.RYD * Cst.H * Cst.C / Cst.Q 

    def func(T):
        x = Eo/(Cst.K * T/ Cst.Q)
        return cst * gi * f / Eo * (1+2/x) / x
    return func

def write_ofile(ofile, odata):
    ''' Write odata into ofile.
    '''
    head = '#' + ' '.join(sys.argv)

    stdout_save = sys.stdout
    output = open(ofile, 'w')
    sys.stdout = output
    print(head)
    print(len(t_list), t_list, file=output)
    for datum in odata:
        datum.print(tot=True)
    output.close()
    sys.stdout = stdout_save

def plot_ups(col):
    ''' Plot the hydrogen collision strength as a function of temperature
    '''
    global iplot
    t = col.t_list
    ups = col.ups_list
    ref = col.ref
    #Name of the output plot file
    ofname = format(iplot, '0>4')+'_ups_h.pdf'    
    #Figure
    pl.figure(dpi=100, facecolor='w', edgecolor='k')
    #Title of the figure
    ll = col.lower
    ul = col.upper
    title = SPE+' '+DEG+': '+str(ll.cfg)+' '+str(ll.term)+' [g='+str(ll.g)+'] '+str(ul.cfg)+' '+str(ul.term)+' [g='+str(ul.g)+'], $\lambda =$'+format(col.lbd, '8.1f')+' $\AA$'
    #Axes
    pl.xlabel('T [K]')
    pl.ylabel('$\Upsilon_H$')
    pl.semilogy()
    pl.ylim([1e-6, 1e5])
    pl.title(title)
    #Plot
    pl.plot(t, ups, 'k', label=col.ref+' ('+col.type.lower()+')')
    pl.legend()
    pl.savefig(ofname, dpi=100, format='pdf', orientation='landscape', papertype='a4')
    pl.close()

def load_act(sym, source):
    ''' Load atomic collision transition file from mechanical quantum computations
        3 ascii input files:
        ael.dat: atomic energy levels in the mael.py format
        temp.dat: list of temperatures
        ups.dat: effective collision strengths
    '''
    path = '/home/tmerle/development/formato2/ad/ct/h/qm/'+SYM+'/'+source+'/'

    dumb = pl.loadtxt(path+'ael.dat', dtype='str')
    ael_list = [State(e=e, g=g, cfg=cfg, term=term, p=p, ref=ref) for e, g, cfg, term, p, ref in dumb]
    #[ael.print() for ael in ael_list]

    t_list = pl.loadtxt(path+'temp.dat', dtype='float')
    nt = len(t_list)

    dumb = pl.loadtxt(path+'ups.dat', dtype='float')
    
    l_list, u_list, ups_list = [], [], []

    for rec in dumb:
        l_list.append(int(rec[0]))
        u_list.append(int(rec[1]))
        ups_list.append(map(float, rec[2:]))            

    col_list = []

    # Corresponding to ionization level
    umax = max(u_list)

    for i, (l, u) in enumerate(zip(l_list, u_list)):
        #print(l, u, ups_list[i])
        if u == umax:
            col = Collision(type='QM', lower=ael_list[l-1], upper=ael_list[u-1], t_list=t_list, ups_list=ups_list[i], ref=source, kw='CHCE')
        else:
            col = Collision(type='QM', lower=ael_list[l-1], upper=ael_list[u-1], t_list=t_list, ups_list=ups_list[i], ref=source, kw='UPS_H')
        col_list.append(col)

    return col_list

if __name__ == '__main__':

    t_list = [1000, 2000, 4000, 6000, 8000, 10000, 15000, 20000]

    IFILE1 = 'ael'
    IFILE2 = 'art_bb'
    OFILE0 = 'act_h'
    OFILE1 = 'mact_h' 
    QM_DATA = ['barklem2012', 'belyaev2015a', 'belyaev2015b', 'belyaev2015b_sigma', 'belyaev2015c', 'belyaev2015c2']

    DESCRIPTION = 'Tool to merge atomic collision transition with neutral hydrogen atoms'
    EPILOG = '2014-06-19 ThiM'
    parser = ap.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
    parser.add_argument('species', type=str, default=None, help='SYMbol of the atomic species (e.g. "Fe" or "FE" or "fe")')
    parser.add_argument('degree', type=str, default=None, help='Ionization degree (e.g. "I" or "i" for neutral species)')
    parser.add_argument('-f1', '--filename1', default='', type=str, help='Filename of energy levels produced by oomael.py')
    parser.add_argument('-f2', '--filename2', default='', type=str, help='Filename of radiative b-b transitions produced by oomart_bb.py')
    parser.add_argument('-lte', '--lte', default=False, action='store_true', help='Enforce collisions to ensure LTE')
    parser.add_argument('-sch', '--semi_classical_h', default=None, type=str, help='Choose semi-classical formula for electron collisions [Drawin]')
    parser.add_argument('-sh', '--scaling_factor', default=1.0, type=float, help='Scaling factor for semi-classic Drawin formula')
    parser.add_argument('-fl', '--forbidden_lines', default=False, action='store_true', help='Include forbidden lines using Drawin formula with f=1')
    parser.add_argument('-qmh', '--quantum_mechanical_h', default=None, type=str, help='Choose quantum mechanical data if available for hydrogen collisions')
    parser.add_argument('-t', '--temperature', nargs='*', default=t_list, help='Temperature list')
    parser.add_argument('-de', '--delta_e', type=float, default=50, help='allowed energy gap [cm⁻¹] for identification with quantum mechanical data (useful when dealing with fine structure')
    parser.add_argument('-p', '--plot', default=False, action='store_true', help='Plot the effective collisions strengths')

    arguments = parser.parse_args()

    SPE = arguments.species
    DEG = arguments.degree
    IF1 = arguments.filename1
    IF2 = arguments.filename2
    LTE = arguments.lte
    SCH = arguments.semi_classical_h
    S_H = arguments.scaling_factor
    FL = arguments.forbidden_lines
    QMH = arguments.quantum_mechanical_h
    t_list = arguments.temperature
    DE = arguments.delta_e
    PLOT = arguments.plot

    print('__________________')
    print('|    mact_h.py   |')
    print('TTTTTTTTTTTTTTTTTT')

    # Building of the element extension name
    SYM = str(SPE).lower() + str(DEG).lower()
    print("\nIon considered:", SPE, DEG)

    # Check the existence of inputs files
    if IF1:
        IFILE1 = IF1
    else:
        IFILE1 = IFILE1+'_'+SYM+'.bin'
    if IF2:
        IFILE2 = IF2
    else:
        IFILE2 = IFILE2+'_'+SYM+'.bin'
    if not os.path.isfile(IFILE1): 
        print("File", IFILE1, "does NOT exist.")
        quit(1)
    if IF2 and not os.path.isfile(IFILE2): 
        print("File", IFILE2, "does NOT exist.")

    # Read input level data
    ifile = open(IFILE1, 'rb')
    levels_list = pickle.load(ifile)
    ifile.close()

    # Read input radiative bb transitions
    ifile = open(IFILE2, 'rb')
    lines_list = pickle.load(ifile)
    ifile.close()

    nlev = len(levels_list)
    ncol = int((nlev-1) * ((nlev-1)-1)/2)
    nlin = len(lines_list)

    print('\nKind of H bb collisions:', QMH, SCH)
    print('Kind of H bf collisions:', 'Drawin')

    print('\nNumber of energy levels:          ', str(nlev).rjust(6))
    print('\nNumber of H collision transitions:', format(ncol, '>6'))
    print(' with oscillator strengths:', format(nlin, '>6'))
    print(' forbidden transitions:    ', format(ncol-nlin, '>6'))

    # Hydrogen semi-classical collisions
    qm_col_list = []
    lu_qm_col_list = [] # list of the level indexes
    allow_col_list = []
    forbd_col_list = []

    cst = Cst.H * Cst.C / Cst.Q 
    mass = Cst.MASS[SPE.upper()]

    if QMH:
        # Be careful
        # Treatment implemented only for qm data relative to mean levels
        # Manage distribution over fine levels included in list of energy levels selected
        # Not yet implemented for the reverse i.e. for qm fine level data 

        nqmcs = 0 # Number of qm collisions selected

        print("\nQuantum mechanical hydrogen collisions:", QMH)

        if QMH in QM_DATA:
            init_qm_col_list = load_act(SYM, QMH)

            for init_qm_col in init_qm_col_list:
                idx_l, idx_u = [], []
                # Select lower levels/states from the model related to a given QM collisions
                for ll in levels_list:
                    crit1 = abs(init_qm_col.lower.e - ll.e) <= DE
                    crit2 = init_qm_col.lower.cfg == ll.cfg 
                    crit3 = init_qm_col.lower.term == ll.term
                    crit4 = init_qm_col.lower.g >= ll.g
                    if crit1 and crit2 and crit3 and crit4:
                        idx_l.append(levels_list.index(ll))
                # Select upper levels/states from the model related to a given QM collisions
                for ul in levels_list:
                    crit1 = abs(init_qm_col.upper.e - ul.e) <= DE
                    crit2 = init_qm_col.upper.cfg == ul.cfg 
                    crit3 = init_qm_col.upper.term == ul.term
                    crit4 = init_qm_col.upper.g >= ul.g                    
                    if crit1 and crit2 and crit3 and crit4:
                        idx_u.append(levels_list.index(ul))

                # Affect a qm collision for levels/states in the energy list
                # with a weighting according statistical weight for distribution over fine levels if necessary
                if idx_l and idx_u:
                    nqmcs += 1
                    #print(idx_l, idx_u, end=' ')
                    gl = sum([levels_list[i].g for i in idx_l])                  
                    gu = sum([levels_list[i].g for i in idx_u])
                    #print(gl, gu)
                    for i in idx_l:
                        for j in idx_u:
                            ups = [levels_list[i].g*levels_list[j].g/(gl*gu) * float(upsh) for upsh in init_qm_col.ups_list]
                            col = Collision(type=init_qm_col.type, lower=levels_list[i], upper=levels_list[j], t_list=init_qm_col.t_list, ups_list=ups, kw=init_qm_col.kw, ref=QMH)
                            col.l = i+1
                            col.u = j+1
                            qm_col_list.append(col)
                            lu_qm_col_list.append((i+1, j+1))
                            iplot += 1
                            if PLOT:
                                plot_ups(col)
        else:
            print("QMH is not implemented.")
            quit(1)                  

        qm_col_list = sorted(qm_col_list, key=op.attrgetter('l', 'u', 'lbd'))

        print(" Number of initial collisions:            ", len (init_qm_col_list))
        print(" Number of initial unselected collisions: ", len(init_qm_col_list)-nqmcs)
        print(" Number of initial selected collisions:   ", nqmcs)
        print(" Number of initial splitted selected collisions:   ", len(qm_col_list))        

    if SCH:
        nsccs = 0 # Number of semi-classic collision selection
        nrfcs = 0  # Number of radiatively forbidden collision selection 
        print("\nSemi-classical hydrogen collisions:", SCH)
        print('Temperatures selected:', t_list)

        # Allowed transitions
        if SCH == 'Drawin':
            for line in lines_list:
                if (line.l, line.u) not in lu_qm_col_list:
                    Eo = cst / (line.lbd * 1e-10)
                    gi = line.lower.g
                    f = line.f
                    # upsh is now a function
                    upsh = ups_sch_drawin(gi, f, Eo, mass)
                    ups = [S_H*upsh(float(t)) for t in t_list]
                    col = Collision(type='ALLOWED', lower=line.lower, upper=line.upper, t_list=t_list, ups_list=ups, lbd=line.lbd, kw='UPS_H', ref=SCH)
                    col.l = line.l
                    col.u = line.u
                    allow_col_list.append(col)
                    iplot += 1
                    if PLOT:
                        plot_ups(col)
                    nsccs += 1
        # Forbidden transitions
        if FL:
            for i in range(1, nlev):
                # List of index (low, up) energy levels
                idx = [(line.l, line.u) for line in lines_list if line.l == i]
                for j in range(i+1, nlev):
                    if (i, j) not in idx+lu_qm_col_list:
                        ei = levels_list[i-1]
                        ej = levels_list[j-1]
                        gi = ei.g
                        Eo = cst * (ej.e-ei.e)*100.
                        upsh = ups_sch_drawin(gi, 1.0, Eo, mass)
                        ups = [S_H*upsh(float(t)) for t in t_list]
                        col = Collision(type='FORBIDDEN', lower=ei, upper=ej, t_list=t_list, ups_list=ups, kw='UPS_H', ref='DEFAULT')
                        col.l = i
                        col.u = j
                        forbd_col_list.append(col)
                        nrfcs += 1

        print(' Number of collisions a la ', SCH, ':', nsccs)
        print(' Number of radiatively forbidden collisions:', nrfcs)
    # Merge QMH and SCH
    # TODO

#    for allow_col in allow_col_list:
#       allow_col.print()
#    for forbd_col in forbd_col_list:
#        forbd_col.print()

# Hydrogen bound-free/charge-exchange collisions

    bf_col_list = []

# Merge bb allowed, forbidden and bf H collisions
    col_list = qm_col_list+allow_col_list+forbd_col_list
    col_list = sorted(col_list, key=op.attrgetter('l', 'u'))
    col_list += bf_col_list




    OFILE1 = OFILE1+'_'+SYM+'.dat'
    write_ofile(OFILE1, col_list)
    print('\nOutput ascii file: ', OFILE1)

# Write internal file
    OFILE0 = OFILE0+'_'+SYM+'.bin'
    ofile = open(OFILE0, 'wb')
    pickle.dump(col_list, ofile, protocol=1)
    ofile.close()
    print('Output binary file:', OFILE0)