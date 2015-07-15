#!/usr/bin/python
# -*- coding: utf-8 -*-
''' Code to plot Grotrian diagram of a given species in a given degree of ionization.
'''

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
#from __future__ import unicode_literals

import argparse as ap
import pickle
import os

try:
    import pylab as pl
except ImportError:
    raise ImportError('pylab module is not installed.')

from mad import Cst, State, Level, SuperLevel, HyperLevel, L2num, term2num, term2mult, Line, Multiplet, none2nan, ions

if __name__ == '__main__':

    DESCRIPTION = 'Tool to plot Grotrian diagram.'
    EPILOG = '2014-03-24 ThiM'
    IFILE1 = 'ael'
    IFILE2 = 'art_bb'
    IFILE3 = 'act_e'
    OFILE = 'grot'
    IL = True # Ionization Level present in levels
    EV_TO_CM = Cst.Q / 100. / Cst.H / Cst.C
    CT_THEO = False # Plot theoretical collision transitions if True

    parser = ap.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
    parser.add_argument('species', type=str, default=None, help='Symbol of the atomic species (e.g. "Fe" or "FE" or "fe")')
    parser.add_argument('degree', type=str, default=None, help='Ionization degree (e.g. "I" or "i" for neutral species)')
    parser.add_argument('-f1', '--filename1', default='', type=str, help='Filename of energy levels produced by mael.py')
    parser.add_argument('-f2', '--filename2', default='', type=str, help='Filename of radiative b-b transitions produced by mart_bb.py')
    parser.add_argument('-f3', '--filename3', default='', type=str, help='Filename of collision transitions produced by macte_e.py')
    parser.add_argument('-rt', '--rad_transitions', default=False, action='store_true', help='Plot radiative transitions')
    parser.add_argument('-ct', '--col_transitions', default=False, action='store_true', help='Plot collisions transitions')
    parser.add_argument('-pcfg', '--parent_configuration', default=False, action='store_true', help='Plot Grotrian diagram as a function of parent configurations (default: terms)')
    parser.add_argument('-t', '--title', default='', type=str, help='Title of the figure')
    parser.add_argument('-s', '--size', nargs='*', type=float, default=[6, 9], help='Size of the graphic [xsize, ysize]')
    parser.add_argument('-e', '--ext', default='png', type=str, help='Type of extension [png | eps | pdf] (default: png)')
    parser.add_argument('-l', '--level_name', type=float, default=0, help='Energy in eV below which level names are labelled')
    parser.add_argument('-g', '--grid', default=False, action='store_true', help='Plot the grid')
    parser.add_argument('-lbd', '--wavelength', type=float, default=0, help='Energy in eV below which wavelengths are labelled')
    parser.add_argument('-nop', '--no_plot', default=False, action='store_true', help='Do not display interactive plot')
    parser.add_argument('-o', '--output', default=None, help='Output plot name with extension [png | eps | pdf]')


    args = parser.parse_args()

    SPE = args.species
    DEG = args.degree
    IF1 = args.filename1
    IF2 = args.filename2
    IF3 = args.filename3
    RT = args.rad_transitions
    CT = args.col_transitions
    PCFG = args.parent_configuration
    TITLE = args.title
    FIG_SIZE = args.size
    EXT = args.ext
    LN = args.level_name
    GRID = args.grid
    LBD = args.wavelength
    NOP = args.no_plot
    OFN = args.output

    print('__________________')
    print('|    pgrot.py    |')
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

# Check the existence of input level file
    if IF1:
        IFILE1 = IF1
    else:
        IFILE1 = IFILE1+'_'+SYM+'.bin'
    if RT:
        if IF2:
            IFILE2 = IF2
        else:
            IFILE2 = IFILE2+'_'+SYM+'.bin'
    if CT:
        if IF3:
            IFILE3 = IF3
        else:
            IFILE3 = IFILE3+'_'+SYM+'.bin'
   
    if not os.path.isfile(IFILE1): 
        print("File", IFILE1, "does NOT exist.")
        quit(1)
    if RT and not os.path.isfile(IFILE2): 
        print("File", IFILE2, "does NOT exist.")
        quit(1)
    if CT and not os.path.isfile(IFILE3):
        print("File", IFILE3, "does NOT exist.")
        quit(1)

# Read input level data
    ifile = open(IFILE1, 'rb')
    levels = pickle.load(ifile)
    ifile.close()

# Read input radiative b-b transitions data
    if RT:
        ifile = open(IFILE2, 'rb')
        lines = pickle.load(ifile)
        ifile.close()

#Read input collisions transition data
    if CT:
        ifile = open(IFILE3, 'rb')
        cols = pickle.load(ifile)
        ifile.close()

    print('Input energy levels file:            ', IFILE1)
    if RT:
        print('Input radiative b-b transitions file:', IFILE2)
    if CT:
        print('Input collisions b-b transition file:', IFILE3)

    BLA = False
    if BLA:
        for level in levels:
            level.print()

# Remove ionization stage to avoid duplicity
    if hasattr(ION, 'states'):
        for state in ION.states:
            if state in levels:
                try:
                    levels.remove(state)
                except ValueError:
                    levels.remove(ION)
            else:
                print("Warning: in ionization Level of", SPE, DEG, "the State component")
                state.print()
                print("not found in", IFILE1)
                quit(1)
    else:
        if ION in levels:
            levels.remove(ION)
        else:
            print("Warning: ionization State of", SPE, DEG, "not found in", IFILE1)
            IL = False


    terms = [level.remove_label_term(level.term) for level in levels if (level.term or type(level.term) != list)]
    terms = list(set(pl.flatten(terms)))
    terms = sorted(terms, key=term2num)
    terms = sorted(terms, key=term2mult)
    print('Number of terms:', len(terms))
    print(terms)

# Extract levels by class
    fl = []
    ml = []
    sl = []
    hl = []

    for level in levels:
        if isinstance(level, HyperLevel) and hasattr(level, 'superlevels'):
            hl.append(level)
        elif isinstance(level, SuperLevel) and hasattr(level, 'levels'):
            sl.append(level)
        elif isinstance(level, Level) and hasattr(level, 'states'):
            ml.append(level)
        else:
            fl.append(level)

    print('Number of fine levels:', len(fl))
    print('Number of mean levels:', len(ml))
    print('Number of superlevels:', len(sl))
    print('Number of hyperlevels:', len(hl))

# Definition of the plot
    #pl.xticks(pl.arange(len(terms))+1, (terms))
    pl.figure(figsize=FIG_SIZE, dpi=100, facecolor='w', edgecolor='k')
    #ax1 = fig.add_subplot(111)
    #ax2 = ax1.twiny()
    pl.ylabel('Energy [eV]', fontsize=12)
    if len(terms) > 20:
        pl.rc('font', size=10)
        pl.xticks(pl.arange(len(terms))+1, (terms), rotation=60)
    else:
        pl.xticks(pl.arange(len(terms))+1, (terms))
    #ax1.set_xticks(pl.arange(len(terms))[::2]+1)
    #ax1.set_xticklabels(terms[::2])ls
    #ax2.set_xticks(pl.arange(len(terms))[::2]+1.5)
    #ax2.set_xticklabels(terms[1:len(terms)-1:2])
    #ax1.tick_params(axis='x', which='major', labelsize=10)   
    pl.xlim(0.5, len(terms)+0.5)
    if levels[0].e/EV_TO_CM == 0:
        pl.ylim((levels[0].e/EV_TO_CM-0.2), round(ION.e/EV_TO_CM+0.5))
        #pl.ylim((levels[0].e/EV_TO_CM-0.2), round(levels[-1].e/EV_TO_CM+0.5))

    if not IL:
        pl.ylim(-0.2, round(max(levels).e/EV_TO_CM+0.5))

    try:
        pl.ticklabel_format(useOffset=False)
    except AttributeError:
        pass

# Plot collision transitions

    if CT_THEO:
        for i, llev in enumerate(fl+ml+sl+hl):
            for j, ulev in enumerate((fl+ml+sl+hl)[i+1:]):
                term_l = llev.remove_label_term(llev.term)
                term_u = ulev.remove_label_term(ulev.term)
                e_l = llev.e/EV_TO_CM
                e_u = ulev.e/EV_TO_CM                
                #print(terms.index(term_l)+1, terms.index(term_u)+1, [e_l, e_u])
                try:
                    ct_theo_line, = pl.plot([terms.index(term_l)+1, terms.index(term_u)+1], [e_l, e_u], 'k-')
                except ValueError: # To avoid error for collisional ionization
                    pass
        ncol_theo = int(len(levels)*(len(levels)-1)/2)
        print('Number of theoretical collision transitions (without ionization collisions):', ncol_theo)

    if CT:
        for col in cols:
            term_l = col.lower.remove_label_term(col.lower.term)
            term_u = col.upper.remove_label_term(col.upper.term)
            e_l = col.lower.e/EV_TO_CM
            e_u = col.upper.e/EV_TO_CM
            try:
                #pl.plot([terms.index(term_l)+1, terms.index(term_u)+1], [e_l, e_u], 'b-')
                ct_line, = pl.plot([terms.index(term_l)+1, terms.index(term_u)+1], [e_l, e_u], '-', color='0.75')
            except ValueError: # To avoid error for collisional ionization
                pass
        print('Number of collision transitions:', len(cols)-len(levels))


    pl.axhline(y=ION.e/EV_TO_CM, linewidth=0.5, linestyle='--', color='0.75')

    if GRID:
        pl.grid()

    fl_line, ml_line, sl_line, hl_line = None, None, None, None

    if PCFG:
        # To do
        print("'pcfg' not yet a valid option. To do.")
        quit(1)
    else:
        # Plot energy levels
        for lev in fl:
            if not lev.term:
                fl_line = pl.axhline(y=lev.e/EV_TO_CM, linewidth=0.5, linestyle='-', color='c')
            else:
                term = lev.remove_label_term(lev.term)
                fl_line, = pl.plot(terms.index(term)+1, lev.e/EV_TO_CM, 'c_', markersize=12)   
        for lev in ml:
            if type(lev.term) == list or not lev.term:      
                ml_line = pl.axhline(y=lev.e/EV_TO_CM, linewidth=0.5, linestyle='-', color='k')
            else:
                term = lev.remove_label_term(lev.term)
                ml_line, = pl.plot(terms.index(term)+1, lev.e/EV_TO_CM, 'k_', markersize=12)
        for lev in sl:
            if type(lev.term) == list or not lev.term:
                sl_line = pl.axhline(y=lev.e/EV_TO_CM, linewidth=0.5, linestyle='-', color='m')
            else:
                term = lev.remove_label_term(lev.term)
                sl_line, = pl.plot(terms.index(term)+1, lev.e/EV_TO_CM, 'm_', markersize=12)       
        for lev in hl:
            if type(lev.term) == list or not lev.term:
                hl_line = pl.axhline(y=lev.e/EV_TO_CM, linewidth=0.5, linestyle='-', color='y')
            else:
                term = lev.remove_label_term(lev.term)
                hl_line, = pl.plot(terms.index(term)+1, lev.e/EV_TO_CM, 'y_', markersize=12)   

# Plot level names
    if LN:
        size = 6
        for lev in fl+ml+sl+hl:
            if lev.e/EV_TO_CM < LN:
                if len(lev.parent_cfg_wo_term()) > 1:
                    text_cfg = ' '.join(lev.parent_cfg_wo_term())
                    pl.text(0.5, lev.e/EV_TO_CM, text_cfg, fontsize=size, color='grey')
                if type(lev.term) == list or not lev.term:
                    text_cfg = ' '.join(lev.parent_cfg_wo_term())
                    text_term = ' '.join(lev.term)
                    pl.text(0.5, lev.e/EV_TO_CM, text_cfg+' '+text_term, fontsize=size, color='grey')
                else:
                    term = lev.remove_label_term(lev.term)
                    pl.text(terms.index(term)+1, lev.e/EV_TO_CM, lev.parent_cfg_wo_term()[0], fontsize=size, color='grey')

# Plot radiative transitions
    if RT:
        for line in lines:
            if isinstance(line, Multiplet):
                if hasattr(line, 'lines'):
                    for idx, line in enumerate(line.lines):
                        term_l = line.lower.remove_label_term(line.lower.term)
                        term_u = line.upper.remove_label_term(line.upper.term) 
                        e_l = line.lower.e/EV_TO_CM
                        e_u = line.upper.e/EV_TO_CM
                        x_l = terms.index(term_l)+1
                        x_u = terms.index(term_u)+1
                        rt_line, = pl.plot([x_l, x_u], [e_l, e_u], 'r-')  
                         # Display lambda value on the radiative line
                        if e_l < LBD:
                            slope = (e_u-e_l)/(x_u-x_l)
                            angle = pl.arctan(slope)*180./pl.pi
                            offset = 0.15*idx
                            x_c, y_c = (x_l+x_u)/2. + offset, (e_l+e_u)/2. + offset
                            trans_angle = pl.gca().transData.transform_angles(pl.array((angle,)), pl.array((x_c, y_c)).reshape((1, 2)))[0]
                            pl.text(x_c, y_c, format(line.compute_lbd(med='air'), '.1f'), fontsize=6, rotation=trans_angle, horizontalalignment='center', verticalalignment='center')                       
            elif isinstance(line, Line):
                term_l = line.lower.remove_label_term(line.lower.term)
                term_u = line.upper.remove_label_term(line.upper.term) 
                e_l = line.lower.e/EV_TO_CM
                e_u = line.upper.e/EV_TO_CM
                x_l = terms.index(term_l)+1
                x_u = terms.index(term_u)+1
                rt_line, = pl.plot([x_l, x_u], [e_l, e_u], 'r-')
                # Display lambda value on the radiative line
                if e_l < LBD:
                    slope = (e_u-e_l)/(x_u-x_l)
                    angle = pl.arctan(slope)*180./pl.pi
                    x_c, y_c = (x_l+x_u)/2., (e_l+e_u)/2.
                    trans_angle = pl.gca().transData.transform_angles(pl.array((angle,)), pl.array((x_c, y_c)).reshape((1, 2)))[0]
                    pl.text(x_c, y_c, format(line.compute_lbd(med='air'), '.1f'), fontsize=6, rotation=trans_angle, horizontalalignment='center', verticalalignment='center')
                      

        print('Number of radiative transitions:', len(lines))

    #pl.text(len(terms)-1, 1, SPE+' '+DEG, fontsize=25, ha='right')

# Set a title to the plot
    if TITLE:
        pl.title(TITLE)

# Set a legend
    lgd_lines, lgd_labels =  [], []

    if fl_line:
        lgd_lines.append(fl_line)
        lgd_labels.append(' fine levels ('+str(len(fl))+')')
    if ml_line:
        lgd_lines.append(ml_line)
        lgd_labels.append('mean levels ('+str(len(ml))+')')
    if sl_line:
        lgd_lines.append(sl_line)
        lgd_labels.append('super levels ('+str(len(sl))+')')
    if hl_line:
        lgd_lines.append(hl_line)
        lgd_labels.append('hyper levels ('+str(len(hl))+')')
    if RT:
        lgd_lines.append(rt_line)
        lgd_labels.append('radiative tr ('+str(len(lines))+')')
    if CT_THEO:
        lgd_lines.append(ct_theo_line)
        lgd_labels.append('theoretical ('+str(ncol_theo)+')')
    if CT:
        lgd_lines.append(ct_line)
        lgd_labels.append('calculated ('+str(len(cols)-len(levels))+')')       
    leg = pl.legend(lgd_lines, lgd_labels, frameon=False, loc=0, numpoints=1, title=SPE+' '+DEG)
    pl.setp(leg.get_title(), fontsize=25)

# Save Grotrian plot
    OFILE = OFILE + '_' + SYM + '.' + EXT

    if OFN:
        EXT = os.path.basename(OFN).split('.')[1]
        OFILE = OFN
        pl.savefig(OFILE, dpi=100, format=EXT, orientation='landscape', bbox_inches='tight')
    elif EXT in ['png', 'eps', 'pdf']:
        pl.savefig(OFILE, dpi=100, format=EXT, orientation='landscape', bbox_inches='tight')
    else:
        print('Extension not yet supported. Try [eps | png | pdf]')
        quit(1)

    if not NOP:
        pl.show()
 
