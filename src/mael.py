#!/usr/bin/python
# -*- coding: utf-8 -*-
''' Code to merge atomic energy level of a given species in a given ionization degree.
'''

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
#from __future__ import unicode_literals

import argparse as ap
import operator as op
import pickle
import sys
import os

try:
    import pylab as pl
except ImportError:
    raise ImportError('pylab module is not installed.')

from mad import Cst, ions, State, Level, SuperLevel, HyperLevel, term2num

def load_ifile(ifile):
    ''' Loading input atomic energy level file.
    '''
    try:
        dumb = pl.loadtxt(ifile, dtype='str')
    except ValueError:
        print("Problem reading input file")
        quit(1)

    if dumb.shape[1] == 3:
        dumb = [State(e=e, g=g, cfg=cfg) for e, g, cfg in dumb]
    elif dumb.shape[1] == 4:
        dumb = [State(e=e, g=g, cfg=cfg, term=term) for e, g, cfg, term in dumb]
    elif dumb.shape[1] == 5:
        dumb = [State(e=e, g=g, cfg=cfg, term=term, p=p) for e, g, cfg, term, p in dumb]
    elif dumb.shape[1] == 6:
        dumb = [State(e=e, g=g, cfg=cfg, term=term, p=p, ref=ref) for e, g, cfg, term, p, ref in dumb]
    else:
        print("Input format problem: shape = ", dumb.shape[1])
        quit(1)

    return sorted(dumb, key=op.attrgetter('e'))

def print_ifile_data(ifile, nfl, mlt_set, iterms, nsobm, ng0):
    ''' Display info about initial fine levels.
    '''
    print("\nInput file:", ifile, "\n")
    print(" Total number of levels in input file:  ", str(nfl).rjust(4))
    print(" Number of unuseful levels (with g = 0):", str(ng0).rjust(4), "\n")
    print(" Number of multiplicity systems:", len(mlt_set))
    print(" Multiplicity systems:          ")#, mlt_set)
    for idx, val in enumerate(mlt_set):
        val = str(val).__format__('4s')
        print('  system', val, ':', sorted(iterms[idx], key=term2num))
    print("\n Number of fine levels: ", nfl-ng0)
    for idx, val in enumerate(mlt_set):
        val = str(val).__format__('4s')
        ns = str(nsobm[idx]).rjust(4)
        print('  system', val, ':', ns, 'states')

def load_mean_idata(state_list, cfg_term_set):
    ''' Compute mean level from input data.
    '''

    level_list = []

    for cfg_term in cfg_term_set:
        state_sel_list = []
        cfg = cfg_term.split()[0]
        if len(cfg_term.split()) == 1:
            term = ''
        else:
            term = cfg_term.split()[1]
        for state in state_list:
            if state.get()[2] == cfg and state.get()[3] == term:
                state_sel_list.append(state)
        level_list.append(Level(*state_sel_list))

    return sorted(level_list, key=op.attrgetter('e'))

def print_mean_idata(cfg_term_set, mlt_set, nlobm):
    ''' Display info about initial mean levels.
    '''
    print("\n Number of mean levels: ", len(cfg_term_set))
    for idx, val in enumerate(mlt_set):
        val = str(val).__format__('4s')
        nl = str(nlobm[idx]).rjust(4)
        print("  system", val, ":", nl, 'mean levels')

def load_sl_idata(level_list, cfg_p_set):
    ''' Compute superlevel from input data.
    '''

    sl_list = []

    for cfg_p in cfg_p_set:
        level_sel_list = []
        cfg = cfg_p.split()[0]
        if len(cfg_p.split()) == 1:
            p = ''
        else:
            p = cfg_p.split()[1]
        for level in level_list:
            if level.get()[2] == cfg and level.get()[4] == p:
                level_sel_list.append(level)
        sl_list.append(SuperLevel(*level_sel_list))

    return sorted(sl_list, key=op.attrgetter('e'))

def print_sl_idata(nsl, nsl_ec):
    ''' Display info about initial superlevels.
    '''
    print("\n Number of superlevels: ", nsl)
    print("  with configuration:   ", str(nsl_ec).rjust(4))
    print("  without configuration:", str(nsl-nsl_ec).rjust(4))

def load_hl_idata(superl_list, DELTA_E, TCOP):
    ''' Compute hyperlevel from input data
    '''

    hyperl_list = [] # List of hyperlevels (merge of superlevels of the same parity with energy difference lower than delta_e)
    sl_temp_list = superl_list[:] # Copy of the list
    sl_select_list = []

    # Sort by parity then by energy
    sl_temp_list = sorted(sl_temp_list, key=op.attrgetter('p', 'e'))

    while True:
        n = len(sl_temp_list)
        if n == 0:
            hyperl_list.append(HyperLevel(TCOP, *sl_select_list))
            break
        if n == 1 and not sl_select_list:
            hyperl_list.append(HyperLevel(TCOP, sl_temp_list[0]))
            break
        if sl_select_list:
            e_sl_list = [sl.e for sl in sl_select_list]
            p_sl_list = [sl.p for sl in sl_select_list]
            sl1 = SuperLevel()
            sl1.e = pl.mean(e_sl_list)
            sl1.p = set(p_sl_list)
            sl1.p = list(sl1.p)[0]
            sl2 = sl_temp_list[0]
        else:
            sl1 = sl_temp_list[0]
            sl2 = sl_temp_list[1]
        if TCOP:
            criteria = abs(sl1.e-sl2.e) < DELTA_E and sl1.p == sl2.p
        else:
            criteria = abs(sl1.e-sl2.e) < DELTA_E
        if criteria:
            if sl_select_list:
                sl_temp_list.__delitem__(0)
                sl_select_list.append(sl2)
            else:
                sl_temp_list.__delitem__(0)
                sl_temp_list.__delitem__(0)
                sl_select_list.extend([sl1, sl2])
        else:
            if sl_select_list:
                hyperl_list.append(HyperLevel(TCOP, *sl_select_list))
                #[sl.print() for sl in sl_select_list]
                sl_select_list = []
            else:
                sl_temp_list.__delitem__(0)
                #print('sl1',end=''); sl1.print()
                hyperl_list.append(HyperLevel(TCOP, sl1))
                #HyperLevel(DELTA_E, sl1).print()

    return sorted(hyperl_list, key=op.attrgetter('e'))

def print_hl_idata(nhl, nhl_ec):
    ''' Display info about initial hyperlevels.
    '''
    print("\n Number of hyperlevels: ", nhl)
    print("  with different configuration:", str(nhl_ec).rjust(4))
    print("  with unique  configuration:  ", str(nhl-nhl_ec).rjust(4))

def write_ofile(ofile, e_list, tot=True):
    ''' Write e_list into ofile.
    '''
    head = '#' + ' '.join(sys.argv)
    head += '\n#E[cm⁻¹]      g  configuration                       term   ref'

    # Save standard output stream
    stdout_save = sys.stdout
    output = open(ofile, 'w')
    # Define a new standard output
    sys.stdout = output

    print(head)

    for idx, e in enumerate(e_list):
        idx = str(idx+1).rjust(4)
        print(idx, end='')
        e.print(tot)

    output.close()

    # Restore standard output
    sys.stdout = stdout_save

    print('Output ascii file: ', ofile)

def write_ofile2(ofile, odata):
    ''' Write odata (term/pcfg) into ofile.
    '''
    pl.savetxt(ofile, odata, fmt='%s')

if __name__ == '__main__':

    OFILE0 = 'ael'
    OFILE1 = 'mael_fl'
    OFILE2 = 'mael_ml'
    OFILE3 = 'mael_sl'
    OFILE4 = 'mael_hl'
    OFILE5 = 'mael_final'
    OFILE6 = 'mael_terms'
    OFILE7 = 'mael_pcfg' # Parent configurations
    LIM_EV = 45. # If values of the arguments line is given above this value are considered in eV, en cm⁻¹ otherwise

    DESCRIPTION = 'Tool to merge atomic energy levels.'
    EPILOG = '2014-02-11 ThiM'
    parser = ap.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
    parser.add_argument('species', type=str, default=None, help='Symbol of the atomic species (e.g. "Fe" or "FE" or "fe")')
    parser.add_argument('degree', type=str, default=None, help='Ionization degree (e.g. "I" or "i" for neutral species)')
    parser.add_argument('filename', default=None, help='Input atomic energy level filename [e[cm⁻¹] g config_term |\n e[cm⁻¹] g config term |\n e[cm⁻¹] g config term ref]')
    #parser.add_argument('-o', '--output', type=str, default=None, help='name of the ion for naming the output files')
    parser.add_argument('-fl', '--fl_thres', type=float, default=None, help='energy threshold [eV | cm⁻¹] above which fine levels are merged in mean levels (default: 0)')
    parser.add_argument('-sl', '--sl_thres', type=float, default=None, help='energy threshold [eV | cm⁻¹] above which mean levels are merged into super levels (default: ionization)')
    parser.add_argument('-hl', '--hl_thres', type=float, default=None, help='energy threshold [ev | cm⁻¹] above which superlevels are merged into hyperlevels (default: ionization)')
    parser.add_argument('-d', '--delta_e', type=float, default=None, help='energy gap [eV | cm⁻¹] below which levels are merged in superlevels')
    parser.add_argument('-l', '--e_lim', type=float, default=None, help='energy limit [eV | cm⁻¹] below which energy levels are considered')
    parser.add_argument('-i', '--ionization', action='store_true', default=False, help='Only if -l option is used - Add the ionization level from mad.ions')
    parser.add_argument('-p', '--do_not_take_care_of_parity', default=False, action='store_true', help='Do not take care of parity in merging superlevels')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Write detailed structure in output files')
    parser.add_argument('-t', '--temporary', action='store_true', default=False, help='Write temporary files (mael_fl, mael_ml, mael_sl, mael_hl, mael_terms, mael_pcfg)')
    parser.add_argument('-o', '--offset', type=float, default=None, help='Decrease input energy levels with this offset value.')
    parser.add_argument('-e', '--epsilon', type=float, default=0.01, help='Add an espilon to energy if 2 levels have exactely the same values (default: 0.01)')

    arguments = parser.parse_args()

    # Line command arguments
    SPE = arguments.species
    DEG = arguments.degree
    IFILE = arguments.filename
    #ONAME = arguments.output
    FL_THRES = arguments.fl_thres
    SL_THRES = arguments.sl_thres
    HL_THRES = arguments.hl_thres
    DELTA_E = arguments.delta_e
    E_LIM = arguments.e_lim
    I_ADD = arguments.ionization
    TCOP = not arguments.do_not_take_care_of_parity
    VERBOSE = arguments.verbose
    TEMP = arguments.temporary
    OFFSET = arguments.offset
    EPS = arguments.epsilon

    EV_TO_CM = Cst.Q / 100. / Cst.H / Cst.C

    print('__________________')
    print('|    mael.py     |')
    print('TTTTTTTTTTTTTTTTTT')

    # Check the presence of the ion in the dictionnary
    SYM = str(SPE).lower() + str(DEG).lower()
    print("\nIon considered:", SPE, DEG)

    try:
        #ION = ions.get(SYM)
        ION = ions[SYM]
    except KeyError:
        print('Ion not in the dictionnary mad.ions. Add it in mad.py module.')
        quit(1)
    else:
        print('Ionization stage:', ION.get(True)[0], ' cm⁻¹ (', (ION.e/EV_TO_CM).__format__('6.3f'), 'eV)', end=' ')
        print(ION.cfg, ION.term, ION.p, ION.ref)
        #ION.print()

    # Check the existence of the given line command argument's files
    if not os.path.isfile(IFILE):
        print("File", IFILE, "does NOT exist.")
        quit(1)

# Load input file
    #----------------
    # state_list is a list of State instances
    state_list = load_ifile(IFILE)

    # Source of energy levels
    if IFILE.lower().find('nist') != -1:
        for state in state_list:
            state.ref = 'NIST'
    elif IFILE.lower().find('norad') != -1:
        for state in state_list:
            state.ref = 'NORAD'
    elif IFILE.lower().find('topbase') != -1:
        for state in state_list:
            state.ref = 'TOPBASE'

    #Temporary option to manage input file data from MULTI format model atom
    if OFFSET:
        print("OFFSET:", OFFSET)
        for state in state_list:
            state.e = state.e - OFFSET

    nfl = len(state_list) # Number of initial states

    # Remove ionization stage to avoid duplicity
    # The ionization stage will be had to select_list at the end
    if hasattr(ION, 'states'):
        for state in ION.states:
            if state in state_list:
                state_list.remove(state)
            else:
                print("Warning: in ionization Level of", SPE, DEG, "the State component", state.print(), "not found in", SPE.strip())
                quit(1)
    else:
        if ION in state_list:
            state_list.remove(ION)
        else:
            print("Warning: ionization State of", SPE, DEG, "not found in the input file")
            #quit(1)

    # Species degree of the level
    if DEG.strip().upper() == 'I': i = 1
    if DEG.strip().upper() == 'II': i = 2

    for state in state_list:
        state.i = i

    # State with a g = 0 (undefined statistical weight)
    ng0 = 0
    for state in state_list:
        if state.g == 0:
            ng0 += 1

    # Remove states with g = 0 (undefinde statistical weight)
    state_list = [state for state in state_list if state.g != 0]

    #Set of the different multiplicity systems
    mlt_set = set([state.get_S() for state in state_list])
    #List of set of initial single terms ordered by multiplicity
    iterms = []
    for i in mlt_set:
        iterms.append(set([state.remove_label_term(state.term) for state in state_list if state.get_S() == i]))
    #Number of states ordered by multiplicity
    nsobm = []
    for i in mlt_set:
        #List of states by multiplicity
        sbm = [state for state in state_list if state.get_S() == i]
        nsobm.append(len(sbm))

    # Print input data
    print_ifile_data(IFILE, nfl, mlt_set, iterms, nsobm, ng0)
    if TEMP: 
        OFILE1 = OFILE1+'_'+SYM+'.dat'
        write_ofile(OFILE1, state_list, tot=VERBOSE)

    bla = False
    if bla:
        for state in state_list:
            state.print()

# Compute and load mean level
    #----------------------------
    cfg_list = [state.get()[2] for state in state_list]
    term_list = [state.get()[3] for state in state_list]
    cfg_term_list = [cfg+' '+term for cfg, term in pl.array([cfg_list, term_list]).T]
    cfg_term_set = set(cfg_term_list)

    level_list = load_mean_idata(state_list, cfg_term_set)

    #Print mean level and details
    bla = False
    if bla:
        for level in level_list:
            level.print(tot=True)

    # Characteristics of mean input levels
    # List of set of initial mean terms ordered by multiplicity
    mterms = []
    for i in mlt_set:
        mterms.append(set([level.term for level in level_list if level.term and level.get_S() == i]))
    # Number of levels ordered by multiplicity
    nlobm = []
    for i in mlt_set:
        #List of levels by multiplicity
        lbm = [level for level in level_list if level.get_S() == i]
        nlobm.append(len(lbm))

    print_mean_idata(cfg_term_set, mlt_set, nlobm)
    if TEMP: 
        OFILE2 = OFILE2+'_'+SYM+'.dat'
        write_ofile(OFILE2, level_list, tot=VERBOSE)

# Compute and load superlevel
    #----------------------------
    cfg_list = [level.get()[2] for level in level_list]
    p_list = [level.get()[4] for level in level_list]
    cfg_p_list = [cfg+' '+p for cfg, p in pl.array([cfg_list, p_list]).T]
    cfg_p_set = set(cfg_p_list)

    superl_list = load_sl_idata(level_list, cfg_p_set)

    # Print superlevel and details
    bla = False
    if bla:
        for sl in superl_list:
            print(str(sl.nsl).rjust(3), end='')
            sl.print(tot=True)

    # Characteristics of input super levels
    # Number of superlevel (same configuration)
    nsl = len(cfg_p_set)
    # Number of superlevel with clear electronic configuration
    nsl_ec = 0
    for sl in superl_list:
        if sl.term: nsl_ec += 1

    print_sl_idata(nsl, nsl_ec)
    if TEMP: 
        OFILE3 = OFILE3+'_'+SYM+'.dat'
        write_ofile(OFILE3, superl_list, tot=VERBOSE)

# Compute and load hyperlevel
    #----------------------------
    if DELTA_E:
        if DELTA_E < LIM_EV/10.: DELTA_E = DELTA_E * EV_TO_CM
    else:
        DELTA_E = 0.

    hyperl_list = load_hl_idata(superl_list, DELTA_E, TCOP)
    nhl = len(hyperl_list)
    #Number of hyperlevel resulting of merging of superlevels
    nhl_ec = 0

    for hl in hyperl_list:
        if type(hl.cfg) == list:
            nhl_ec += 1

    print_hl_idata(nhl, nhl_ec)

    if TEMP: 
        OFILE4 = OFILE4+'_'+SYM+'.dat'
        write_ofile(OFILE4, hyperl_list, tot=VERBOSE)

# Preparing selection criteria
    #-----------------------------

    # Convert parameter from eV in cm⁻¹
    if E_LIM and E_LIM < LIM_EV: E_LIM = E_LIM * EV_TO_CM
    if FL_THRES and FL_THRES < LIM_EV: FL_THRES = FL_THRES * EV_TO_CM
    if SL_THRES and SL_THRES < LIM_EV: SL_THRES = SL_THRES * EV_TO_CM
    if HL_THRES and HL_THRES < LIM_EV: HL_THRES = HL_THRES * EV_TO_CM
    #if DELTA_E and DELTA_E < LIM_EV/10.: DELTA_E = DELTA_E * EV_TO_CM

    # Print selection parameters if present
    if E_LIM or FL_THRES or SL_THRES or DELTA_E: print("\nSelection parameters:")
    if E_LIM: print(" Levels selected below                      e_lim =", E_LIM.__format__('10.3f'), " cm⁻¹ (", (E_LIM/EV_TO_CM).__format__('6.3f'), ' eV)')
    if FL_THRES: print(" Fine structure taken into account below fl_thres =", FL_THRES.__format__('10.3f'), " cm⁻¹ (", (FL_THRES/EV_TO_CM).__format__('6.3f'), ' eV)')
    if SL_THRES: print(" Super levels taken into account upper   sl_thres =", SL_THRES.__format__('10.3f'), " cm⁻¹ (", (SL_THRES/EV_TO_CM).__format__('6.3f'), ' eV)')
    if HL_THRES: print(" Hyper levels taken into account upper   hl_thres =", HL_THRES.__format__('10.3f'), " cm⁻¹ (", (HL_THRES/EV_TO_CM).__format__('6.3f'), ' eV)')
    if DELTA_E: print(" Merged super levels with                 delta_e <", DELTA_E.__format__('10.3f'), " cm⁻¹ (", (DELTA_E/EV_TO_CM).__format__('6.3f'), ' eV)')
    if E_LIM or FL_THRES or SL_THRES or DELTA_E: print()

# Fine levels consideration

    fl_list = []  # List of selected fine levels  below fl_thres
    ml1_list = [] # List of mean levels corresponding to non-selected fine levels

    nfl = 0      # Number of selected fine levels
    nml1 = 0     # Number of mean levels larger than fl_thres
    nml1_nfl = 0 # Number of mean levels associated to selected fine levels
    nfl_nml1 = 0 # Number of fine levels associated to the non-selected mean levels

    if FL_THRES:
        for level in level_list:
            if level.e < FL_THRES:
                nml1_nfl += 1
                for i in range(level.nfl):
                    fl_list.append(level.states[i])
                    nfl += 1
            else:
                ml1_list.append(level)
                nml1 += 1
                nfl_nml1 += level.nfl
    else:
        ml1_list = level_list

# Superlevels consideration

    sl_list = []  # List of selected superlevels larger than SL_THRES and lower than HL_THRES
    sl_inf_list = [] # List of superlevels lower than SL_THRES
    ml2_list = [] # List of mean levels corresponding to non-selected superlevels

    nsl = 0      # Number of superlevels
    nml2 = 0     # Number of mean levels lower than SL_THRES
    nml2_nsl = 0 # Number of mean levels larger than SL_THRES corresponding to the selected superlevel
    nsl_nml2 = 0 # Number of superlevels associated to the non-selected mean levels

    if SL_THRES:
        for sl in superl_list:
            if sl.e > SL_THRES:
                if HL_THRES:
                    if sl.e < HL_THRES:
                        sl_list.append(sl)
                        nsl += 1
                        nml2_nsl += sl.nml
                else:
                    sl_list.append(sl)
                    nsl += 1
                    nml2_nsl += sl.nml
            else:
                if sl.nml == 0:
                    ml2_list.append(Level(**sl.__dict__))
                    nml2 += 1
                    nsl_nml2 += 1
                else:
                    for i in range(sl.nml):
                        ml2_list.append(sl.levels[i])
                        nml2 += 1
                        nsl_nml2 += sl.nml
                sl_inf_list.append(sl) # will be used to select hyperlevels
    else:
        ml2_list = level_list
        nml2 = len(level_list)

# Mean levels consideration
    ml_list = []
    nfl_nml = 0 # Number of fine levels corresponding to the mean level selected

    for level in level_list:
        if level in ml1_list and level in ml2_list:
            ml_list.append(level)
            nfl_nml += level.nfl

    nml = len(ml_list) # Number of the selected mean levels

# Hyperlevels consideration
    hl_list = []
    nhl = 0 # Number of hyperlevels
    nsl_nhl = 0 # Number of superlevels corresponding the the selected hyperlevels

#    if HL_THRES:
#        for hl in hyperl_list:
#            if hl.e > HL_THRES:
#                hl_list.append(hl)
#                nsl_nhl += hl.nsl
#
    if HL_THRES:
        sl_rest_list = sorted(set(superl_list)-set(sl_list)-set(sl_inf_list), key=op.attrgetter('e'))
        hl_list = load_hl_idata(sl_rest_list, DELTA_E, TCOP) 
        nsl_nhl = [hl.nsl for hl in hl_list]
        nsl_nhl = sum(nsl_nhl)
    
    nhl = len(hl_list)

    #Maximum energy difference for hyperlevel
    de_max = 0

    for hl in hl_list:
        if hl.de > de_max:
            de_max = hl.de

# Print selection criteria
    if FL_THRES:
        print('\nFine structure criterium')
        print(' Number of fine levels selected:', str(nfl).rjust(4), '(corresponding to', str(nml1_nfl).rjust(4), 'mean levels)')

    print('\nMean levels selection')
    print(' Number of mean levels selected:', str(nml).rjust(4), '(corresponding to', str(nfl_nml).rjust(4), 'fine levels)')

    if SL_THRES:
        print('\nSuper level criterium')
        print(' Number of superlevels selected:', str(nsl).rjust(4), '(corresponding to', str(nml2_nsl).rjust(4), 'mean levels)')

    if HL_THRES:
        print('\nHyper level criterium')
        print(' Number of hyperlevels selected:', str(nhl).rjust(4), '(corresponding to', str(nsl_nhl).rjust(4), 'superlevels)')
        print(" with maximum range de_max:     ", de_max, 'cm⁻¹ (', (de_max/EV_TO_CM).__format__('6.3f'), 'eV)')

# Energy threshold criterium
    fl_lim_list = [] # List of selected fine levels below e_lim
    ml_lim_list = [] # List of selected mean levels below e_lim
    sl_lim_list = [] # List of selected superlevels below e_lim
    hl_lim_list = [] # List of selected hyperlevels below e_lim

    nfl_lim = 0 # Number of selected fine levels below e_lim
    nml_lim = 0 # Number of selected mean levels below e_lim
    nsl_lim = 0 # Number of selected superlevels below e_lim
    nhl_lim = 0 # Number of selected hyperlevels below e_lim

    if E_LIM:
        for fl in fl_list:
            if fl.e <= E_LIM:
                fl_lim_list.append(fl)
                nfl_lim += 1
        for ml in ml_list:
            if ml.e <= E_LIM:
                ml_lim_list.append(ml)
                nml_lim += 1
        for sl in sl_list:
            if sl.e <= E_LIM:
                sl_lim_list.append(sl)
                nsl_lim += 1
        for hl in hl_list:
            if hl.e <= E_LIM:
                hl_lim_list.append(hl)
                nhl_lim += 1

    if E_LIM:
        print('\nEnergy limit criterium')
        print(' Number of selected fine levels:', nfl_lim, '/', nfl)
        print(' Number of selected mean levels:', nml_lim, '/', nml)
        print(' Number of selected superlevels:', nsl_lim, '/', nsl)
        print(' Number of selected hyperlevels:', nhl_lim, '/', nhl)

# Combine the different criteria
    if E_LIM:
        select_list = fl_lim_list + ml_lim_list + sl_lim_list + hl_lim_list
        if I_ADD: 
            select_list += [ION]
        nfl = nfl_lim
        nml = nml_lim
        nsl = nsl_lim
        nhl = nhl_lim
    else:
        select_list = fl_list + ml_list + sl_list + hl_list + [ION]

    ntot = nfl + nml + nsl + nhl

    select_list = sorted(select_list, key=op.attrgetter('e'))

    bla = False
    if bla:
        for select in select_list:
            select.print(tot=True)

    print('\n\nTotal number of selected energy levels:', ntot, '(without ionization level)')
    print(' with number of fine levels:', str(nfl).rjust(4))
    print(' with number of mean levels:', str(nml).rjust(4))
    print(' with number of superlevels:', str(nsl).rjust(4))
    print(' with number of hyperlevels:', str(nhl).rjust(4))

    # User output
    OFILE5 = OFILE5+'_'+SYM+'.dat'
    write_ofile(OFILE5, select_list, tot=VERBOSE)

# Determine terms used
    terms_list = []
    for level in select_list[:ntot-1]:
        if type(level.term) == str:
            terms_list.append(level.term)
        elif type(level.term) == list:
            terms_list.extend(level.term)

# Add an espilon if two levels have the same energy value
    for idx, level in enumerate(select_list[1:ntot-1]):
        if abs(level.e-select_list[idx].e) < 1.e-3:
            select_list[idx].print()
            level.print()
            level.e += EPS
            print("Levels for which an espilon of "+str(EPS)+" were added:", end=' '), level.print()

    terms_set = set(terms_list)
    # Sort by spectral term L
    terms_list = sorted(list(terms_set), key=term2num)
    # Sort by multiplicity
    terms_list = sorted(terms_list)

    if TEMP: 
        OFILE6 = OFILE6+'_'+SYM+'.dat'
        write_ofile2(OFILE6, terms_list)

# Determine parent configuration used
    pcfg_list = []
    for level in select_list[:ntot-1]:
        pcfg_list.extend(level.parent_cfg_wo_term())

    pcfg_set = set(pcfg_list)
    pcfg_list = sorted(list(pcfg_set))

    if TEMP: 
        OFILE7 = OFILE7+'_'+SYM+'.dat'
        write_ofile2(OFILE7, pcfg_list)

# Internal output for reading by other python codes
    OFILE0 = OFILE0+'_'+SYM+'.bin'
    ofile = open(OFILE0, 'wb')
    pickle.dump(select_list, ofile, protocol=1)
    ofile.close()
    print('Output binary file:', OFILE0)

