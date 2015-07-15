#!/usr/bin/python
# -*- coding: utf-8 -*-
''' Code to merge atomic radiative transitions of a given species in a given ionization degree.
'''

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
#from __future__ import unicode_literals

import argparse as ap
import operator as op
import time as t
import pickle
import sys
import os

try:
    import pylab as pl
except ImportError:
    raise ImportError('pylab module is not installed.')

from mad import ions, fh, Cst, State, Line, Multiplet, none2nan
from barklem import retcross

def load_art(ifile1, source):
    ''' Load atomic radiative transition file
    '''
    dumb = pl.loadtxt(ifile1, dtype='str')

    # Case where only one line in the input lines file
    if len(dumb.shape) == 1:
        nfields = dumb.shape[0]
        dumb = [dumb]
    else:
        nfields = dumb.shape[1]

    #print('nfields:', nfields)

    if source.strip() == 'nist':
        if nfields == 12:
            dumb = [Line(lower=State(e=e_l, g=g_l, cfg=cfg_l, term=term_l), \
                              upper=State(e=e_u, g=g_u, cfg=cfg_u, term=term_u), \
                              f=10**(float(none2nan(loggf))) / float(g_l), \
                              f_acc=acc, gr=none2nan(Aki), lbd=lbd) \
                      for lbd, Aki, loggf, acc, e_l, e_u, cfg_l, term_l, cfg_u, term_u, g_l, g_u in dumb]

        elif nfields == 15:
            dumb = [Line(lower=State(e=e_l, g=g_l, cfg=cfg_l, term=term_l), \
                              upper=State(e=e_u, g=g_u, cfg=cfg_u, term=term_u), \
                              f=10**(float(none2nan(loggf))) / float(g_l), \
                              f_acc=acc, f_ref=ref2, gr=none2nan(Aki), lbd=lbd, rtype=rtype, lev_ref=ref1) \
                      for lbd, Aki, loggf, acc, e_l, e_u, cfg_l, term_l, cfg_u, term_u, g_l, g_u, rtype, ref1, ref2 in dumb]
        else:
            print('Load_art: input format problem.')
            quit(1)
    elif source.strip() in  ['vald', 'kald']:
        if nfields == 10:
            dumb = [Line(lower=State(e=e_l, g=g_l, cfg=cfg_l, term=term_l), \
                              upper=State(e=e_u, g=g_u, cfg=cfg_u, term=term_u), \
                              f=f, lbd=lbd) \
                      for e_l, g_l, cfg_l, term_l, e_u, g_u, cfg_u, term_u, f, lbd in dumb]
        # In vald lbd are already in air, so recompute them
        for line in dumb:
            setattr(line, 'lbd', line.compute_lbd(med='air'))
    elif source.strip() == 'multi':
        if nfields == 12:
            dumb = [Line(lower=levels_list[int(i)-1], upper=levels_list[int(j)-1], f=f, gr=ga, gv=gvw, gs=gs, lbd=lbd) \
                    for j, i, f, nq, qmax, q0, iw, ga, gvw, gs, lbd, kr in dumb]
    else:
        print(source, 'not implemented. Try nist | vald | kald | multi instead.')
        quit(1)

    return sorted(dumb, key=op.attrgetter('lower.e', 'upper.e', 'lbd'))

def write_ofile(ofile, ll, index=False, lev=False, tot=False):
    ''' Write linelist ll into ofile.
    '''
    head = '#' + ' '.join(sys.argv)
    if index:
        head += '\n#   J   I    F        GRAD     GVDW     GSTA      lambda[Å] Acc'
    else:
        head += '\n#E[cm⁻¹]      g  configuration                       term   ref'

    # Save standard output stream
    stdout_save = sys.stdout
    output = open(ofile, 'w')
    # Define a new standard output
    sys.stdout = output

    print(head)

    for idx, l in enumerate(ll):
        idx = str(idx+1).rjust(4)
        if index:
            print(str(l.u).rjust(4), str(l.l).rjust(4), end='')
            l.print(tot=tot)
        else:
            print(idx, end='')
            l.print(lev=lev, ret=False)

    output.close()

    # Restore standard output
    sys.stdout = stdout_save

    print('Output ascii file: ', ofile)

if __name__ == "__main__":

    IFILE0 = '/home/tmerle/development/formato2/ad/bp/hb/ionized_species.dat'
    IFILE1 = 'ael'
    OFILE0 = 'art_bb'
    OFILE1 = 'mart_bb_init'
    OFILE2 = 'mart_bb_final'
    OFILE3 = 'mart_bb_ll'
    FH_DEFAULT = 2.0
    GS_DEFAULT = 0.0

    DESCRIPTION = 'Tool to merge atomic bound-bound radiative transitions.'
    EPILOG = '2014-03-10 ThiM'

    EV_TO_CM = Cst.Q / 100. / Cst.H / Cst.C

    parser = ap.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
    parser.add_argument('species', type=str, default=None, help='Symbol of the atomic species (e.g. "Fe" or "FE" or "fe")')
    parser.add_argument('degree', type=str, default=None, help='Ionization degree (e.g. "I" or "i" for neutral species)')
    parser.add_argument('filename', default=None, help='Bound-bound radiative transitions file')
    parser.add_argument('-i', '--input', type=str, default=None, help='Atomic energy levels file (in output format of mael.py - default ael_xxx.bin)')
    parser.add_argument('-s', '--source', type=str, default='nist', help='Source of atomic line database [nist | kald | vald | multi] (default nist)')
    parser.add_argument('-f', '--f_nan', action='store_true', default=False, help='Select also lines with no f values')
    #parser.add_argument('-o', '--output', type=str, default=None, help='name of the ion for naming the output files')
    parser.add_argument('-abo', '--abo_theory', action='store_true', default=False, help='Apply ABO theory where possible.')
    parser.add_argument('-fh', '--enhancement_factor', type=float, default=None, help='enhancement factor for Van der Waals broadening according to Unsold')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Write detailed structure in output files')
    parser.add_argument('-t', '--temporary', action='store_true', default=False, help='Write temporary files (mart_bb_init, mart_bb_ll)')

    arguments = parser.parse_args()

    SPE = arguments.species
    DEG = arguments.degree
    IFILE2 = arguments.filename # Line file
    INP = arguments.input # Level file
    SOURCE = arguments.source
    #ONAME = arguments.output
    F_NAN = arguments.f_nan
    VERBOSE = arguments.verbose
    TEMP = arguments.temporary
    ABO = arguments.abo_theory
    FH = arguments.enhancement_factor

    EV_TO_CM = Cst.Q / 100. / Cst.H / Cst.C

    print('__________________')
    print('|   mart_bb.py   |')
    print('TTTTTTTTTTTTTTTTTT')

# Check the presence of the ion in the dictionnary
    SYM = str(SPE).lower() + str(DEG).lower()
    print("\nIon considered:", SPE, DEG)

    try:
        ION = ions.get(SYM)
    except KeyError:
        print('Ion not in the dictionnary mad.ions. Add it in mad.py module.')
        quit(1)
    else:
        print('Ionization stage:', ION.get(True)[0], ' cm⁻¹ (', (ION.e/EV_TO_CM).__format__('6.3f'), 'eV)', end=' ')
        print(ION.cfg, ION.term, ION.p, ION.ref)

# Input filename
    if INP:
        IFILE1 = INP
    else:
        IFILE1 = IFILE1+'_'+SYM+'.bin'

    # Check the existence of the given line command argument's files
    if not os.path.isfile(IFILE0):
        print("File", IFILE0, "does NOT exist.")
        quit(1)
    if not os.path.isfile(IFILE1): 
        print("File", IFILE1, "does NOT exist.")
        quit(1)
    if not os.path.isfile(IFILE2):
        print("File", IFILE2, "does NOT exist.")
        quit(1)

    # Read input data
    ifile = open(IFILE1, 'rb')
    levels_list = pickle.load(ifile)
    ifile.close()

    init_lines_list = load_art(IFILE2, SOURCE)

    # Species degree of the level
    if DEG.strip().upper() == 'I': i = 1
    if DEG.strip().upper() == 'II': i = 2

    for line in init_lines_list:
        line.lower.i = i
        line.upper.i = i

    # Van der Waals broadening
    if not FH:
        FH = fh.get(SYM.strip())
        if not FH:
            FH = FH_DEFAULT

    for line in init_lines_list:
        line.gv = FH
        # Stark broadening
        line.gs = GS_DEFAULT
        line.__format__()

    nabo = 0

    if ABO:

        # Compute collisions H broadening parameter according to ABO theory
        nwarn = 0
        nlwlhti = 0 # Number of lines with levels higher than ionization
        cross_list = []
        alpha_list = [] 

        if i == 1:
            # Neutral species
            print("\nInterpolate H collisions broadening parameter according to ABO theory...")
            for line in init_lines_list:
                if ION.e-line.lower.e >0. and ION.e-line.upper.e > 0.:
                    nstar_l = line.lower.i * pl.sqrt(Cst.RYD/100./(ION.e-line.lower.e))
                    nstar_u = line.upper.i * pl.sqrt(Cst.RYD/100./(ION.e-line.upper.e))
                    l_l = line.lower.get_sqn()
                    l_u = line.upper.get_sqn()
                    try:
                        cross, alpha, ifail = retcross(nstar_l, nstar_u, l_l, l_u)
                    except TypeError:
                        print("nstar_l, nstar_u, l_l, l_u:", nstar_l, nstar_u, l_l, l_u)
                    else:
                        cross = int(round(cross))
                        alpha = round(alpha, 3)
                        if not ifail:
                            line.gv = cross + alpha
                            nabo += 1
                            cross_list.append(cross)
                            alpha_list.append(alpha)
                else:
                    nlwlhti += 1
        elif i == 2:
            # First ionized species
            print("\nFind possible ABO broadening for excited species in", IFILE0,'...')
            if SPE.strip().upper() in ['BE', 'MG', 'CA', 'SR', 'BA', 'FE', 'CR']:
                eps_lbd = 0.5 # A
                eps_e = 30.0  # cm⁻¹
                abo_data = pl.loadtxt(IFILE0, dtype='S')
                abo_data = abo_data.compress(abo_data.T[1] == SPE.strip().upper(), axis=0)
                for line in init_lines_list:
                    crit1 = abs(pl.array(map(float, abo_data.T[4]))-line.lower.e) < eps_e
                    crit2 = abs(pl.array(map(float, abo_data.T[5]))-line.upper.e) < eps_e
                    dumb = abo_data.compress(pl.logical_and(crit1, crit2), axis=0)
                    #dumb = abo_data.compress(abs(pl.array(map(float, abo_data.T[3]))- line.lbd) < eps_lbd, axis=0)
                    if len(dumb) >= 1:
                        if len(dumb) == 1:
                            dumb = list(dumb.T.flatten())
                            cross = float(dumb[8])
                            alpha = float(dumb[9])
                            line.gv = cross + alpha
                            #print(line.gv)
                            nabo += 1
                        else:
                            # Compute mean of cross and alpha
                            #print("Warning:", len(dumb), dumb.T[8], dumb.T[9], end=' ')
                            cross = pl.array(map(float, dumb.T[8]))
                            cross = cross.mean()
                            alpha = pl.array(map(float, dumb.T[9]))
                            alpha = alpha.mean()
                            line.gv = cross + alpha
                            nwarn += 1
                            nabo += 1
                            #print(line.gv)
                            #raw_input()    
                        cross_list.append(cross)
                        alpha_list.append(alpha)    

    #print(cross_list, alpha_list)

    if TEMP: 
        OFILE1 = OFILE1+'_'+SYM+'.dat'
        write_ofile(OFILE1, init_lines_list, lev=True)


    print("\nAtomic linelist source: ", SOURCE.upper())
    print("\nNumber of lines:        ", len(init_lines_list))

    ll_set = set([line.lower for line in init_lines_list])
    ul_set = set([line.upper for line in init_lines_list])
    nf_tot = len(init_lines_list)
    nf_nan = len([line.f for line in init_lines_list if line.f != line.f])

    print(" Lower single levels:", len(ll_set))
    print(" Upper single levels:", len(ul_set))
    print("\n Number of lines with f values: ", nf_tot - nf_nan)
    print(" Number of lines without f values:", nf_nan)
    if ABO:
        print("\n Number of lines with ABO values:                ", nabo)
        if nwarn:
            print(" Among them with ABO warning (mean of selection):", nwarn)
        if nlwlhti:
            print(" Number of lines with levels higher than the ionization:", nlwlhti)
        if cross_list:
            print(" min cross:", format(min(cross_list), '4.0f'), " max cross:", format(max(cross_list), '4.0f'))
        if alpha_list:
            print(" min alpha:", format(min(alpha_list), '.3f'), " max alpha:", format(max(alpha_list), '.3f'))
    print(" Number of lines with Unsold enhancement factor:", nf_tot - nabo)



    # Selection of lines

    print("\nSelection of lines")

    temp_list = []
    nls = 0 # Number of lines selected

    # Selection of lines related to levels in levels_list
    nlev = len(levels_list)

    print("\n Number of levels selected in mael.py:", nlev, '\n')

    #lines_for_ll_list = [line for line in init_lines_list for ll in levels_list if line.lower == ll]
    #lines_for_ul_list = [line for line in init_lines_list for ll in levels_list if line.upper == ll]
    #temp_list = [line for line in init_lines_list if line in lines_for_ll_list and line in lines_for_ul_list]


    
    #Select lines for low levels 
    print(" Selecting lines arising from low levels selected in mael.py... ", end=' ')
    t1 = t.time()

    # Init Lines List Copy
    illc = list(init_lines_list)

    for i, ll in enumerate(levels_list[:nlev-1]):
        lines_for_ll_list = []
        ul_for_ll_list = []

        # The following list in intension is replaced by the following for loop much faster
        #lines_for_ll_list = [line for line in init_lines_list if line.lower == ll] too long!

        # Pass illc as an argument of list() is MANDATORY to do a copy of illc
        # in order to avoid troubleshooting with the remove method used after
        for line in list(illc):
            if line.lower == ll:
                lines_for_ll_list.append(line)
                illc.remove(line)

        for j, line in enumerate(lines_for_ll_list):
            for k, ul in enumerate(levels_list[i:]):
                if line.upper == ul and ll != ul:
                    ul_for_ll_list.append(ul)
                    #print(i+1, k+1, end=''); line.print()
                    temp_list.append(line)
                    #if i == k:
                    #    print(i+1, k+1, end=''); line.print()  
                    setattr(temp_list[nls], 'l', i+1)
                    setattr(temp_list[nls], 'u', i+k+1)
                    nls += 1

    t2 = t.time()-t1; print(format(t2, '8.3f'), 's')

    temp_list = sorted(temp_list, key=op.attrgetter('l', 'u', 'lbd'))



    print(' List of unselected lines...                                    ', end=' ')
    t1 = t.time()
    #unselected_list = [line for line in init_lines_list if line not in temp_list] too long!
    # or the remaning lines in init_lines_list_copy are the same as unselected_list 
    unselected_list = list(set(init_lines_list) - set(temp_list))
    t2 = t.time()-t1; print(format(t2, '8.3f'), 's')

    print(' List of unselected lines with f values...                      ', end=' ')
    t1 = t.time()
    unsel_w_f_list = [line for line in unselected_list if not pl.isnan(line.f)]
    t2 = t.time()-t1; print(format(t2, '8.3f'), 's', end='\n')
    

    print(" Number of unselected lines:             ", len(unselected_list))
    print(" Number of unselected lines with f value:", len(unsel_w_f_list))
    #for line in unsel_w_f_list:
    #    line.print(lev=True, ret=False)

    # Removed lines without f value
    if not F_NAN:
        temp_list = [line for line in temp_list if not pl.isnan(line.f)]
        nls = len(temp_list)

    #[(print(line.l, line.u, end=''), line.print()) for line in temp_list]

    # Sorted indexes of lower levels of selected lines
    idx_l = [line.l for line in temp_list]
    idx_l = sorted(list(set(idx_l)))
    # Sorted indexes of upper levels of selected lines
    idx_u = [line.u for line in temp_list]
    idx_u = sorted(list(set(idx_u)))

    select_list = []

    print(' Combining selected lines into multiplet where necessary...     ', end=' ')
    t1 = t.time()
    for l in idx_l:
        for u in idx_u:
            temp_linelist = []
            #Selection of all lines with the same lower and upper levels
            # temp_linelist = [line for line in temp_list if line.l == l and line.u == u] too long!
            # Reversed in MANDATORY to avoid troubleshooting with the remove method used 
            for line in reversed(temp_list):
                if line.l == l and line.u == u:
                    temp_linelist.append(line)
                    temp_list.remove(line)

            if temp_linelist:
                if len(temp_linelist) == 1:
                    select_list.append(*temp_linelist)
                elif len(temp_linelist) > 1:
                    multiplet = Multiplet(*temp_linelist)
                    setattr(multiplet, 'l', l)
                    setattr(multiplet, 'u', u)
                    setattr(multiplet, 'lower', levels_list[l-1])
                    setattr(multiplet, 'upper', levels_list[u-1])
                    select_list.append(multiplet)

                else:
                    print("Problem with this line:", end=''); print(l, u, end='')
                    line.print()

    t2 = t.time()-t1; print(format(t2, '8.3f'), 's', end='\n')

    nms = len(select_list)

    print(" Number of selected initial lines:", str(nls).rjust(4))
    print(" Number of selected output lines: ", str(nms).rjust(4))

    #OFILE2 = OFILE2+'_'+SYM+'.dat'
    #write_ofile(OFILE2, select_list, index=True, tot=VERBOSE)

#    for line in select_list:
#        line.print(lev=True)
#        if line.upper.cfg == None:
#            import pdb; pdb.set_trace()

#    def select_lines_for_grad(seq, ll, ul):
#        for el in seq:
#            if el.upper == ll or el.upper == ul:
#                yield el

# Compute radiative damping broadening parameters
    print("\nCompute radiative broadening parameter for each transition...   ", end=' ')
    t1 = t.time()
    for line in select_list:
        ll = line.lower
        ul = line.upper
        #ll.print(), ul.print()
        #temp_linelist = [l for l in select_list if l.upper == ll or l.upper == ul]
        select_array = pl.array(select_list)
        upper_list = pl.array([o.upper for o in select_array])
        idx1 = pl.where(upper_list == ll)[0]
        idx2 = pl.where(upper_list == ul)[0]
        idx = set(idx1).union(idx2)
        #temp_linelist = select_array[list(idx)]
        temp_linelist = select_array.take(list(idx))
        #temp_linelist = select_lines_for_grad(select_list, ll, ul)
        #temp_linelist = filter(lambda l: l.upper == ll or l.upper == ul, select_list)
        #temp_linelist = []
        #for l in select_list:
        #    if l.upper == ll or l.upper == ul:
        #        temp_linelist.append(l)
        grad = 0.
        for l in temp_linelist:
            grad += 4*pl.pi*l.f2A()            
        line.gr = grad 

    t2 = t.time()-t1; print(format(t2, '8.3f'), 's')

# Counts number of selected lines with ABO H broadening parameters
    dumb = [line.gv for line in select_list if line.gv > 10.]
    nabo = len(dumb)


    if ABO:
        print(" Number of lines with ABO H broadening parameters:        ", str(nabo).rjust(4))
    print(" Number of lines with Unsold enhancement factor FH =", FH, ":", str(len(select_list)-nabo).rjust(4),'\n')

    OFILE2 = OFILE2+'_'+SYM+'.dat'
    write_ofile(OFILE2, select_list, index=True, tot=VERBOSE)

    # Write as a linelist sorted by increasing wavelength
    linelist = []
    lbd_sort_select_list = sorted(select_list, key=op.attrgetter('lbd'))
    for line in lbd_sort_select_list:
        loggf = pl.log10(line.lower.g*line.f)
        el = line.lower.e/EV_TO_CM
        eu = line.upper.e/EV_TO_CM
        Jl = (line.lower.g-1)/2.
        Ju = (line.upper.g-1)/2.
        loggr = pl.log10(line.gr)
        if line.gs == 0.:
            loggs = 0.
        else:
            loggs = pl.log10(line.gs)
        linelist.append([line.lbd, loggf, el, Jl, eu, Ju, loggr, loggs, line.gv])
    
    linelist = pl.array(linelist)

    OFILE3 = OFILE3+'_'+SYM+'.dat'    
    head = 'lambda [Å]      loggf    El[eV]    Jl   Eu[eV]    Ju    loggr   loggs      gv'
    if len(linelist):
        pl.savetxt(OFILE3, linelist, fmt='%12.4f %11.4f %8.4f %5.1f %8.4f %5.1f %8.3f %8.3f %9.3f', header=head)
    else:
        print("Warning: no lines selected!")

# Internal output for reading by other python code
    OFILE0 = OFILE0+'_'+SYM+'.bin'
    ofile = open(OFILE0, 'wb')
    pickle.dump(select_list, ofile, protocol=1)
    ofile.close()
    print('Output binary file:', OFILE0)

