#!/usr/bin/python
# -*- coding: utf-8 -*-
''' Code to merge atomic radiative bound-free transitions of a given species in a given ionization degree.
'''

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
#from __future__ import unicode_literals

import argparse as ap
import operator as op
import difflib as dl
import scipy.signal as ss
import struct
import pickle
import sys
import os

try:
    import pylab as pl
except ImportError:
    raise ImportError('pylab module is not installed.')

from mad import ions, eryd, Cst, State, Level, SuperLevel, HyperLevel
from janicki import gaunt_bf
from px import format_px, extract_px

iplot = 0 # For numeroation of PI plots


def gaunt_menzel_and_pickeris(Z, n, Eion):
    ''' Gaunt factor from Menzel & Pickeris (1935)
        n: principal quantum number
        Eion: Ionization energy [cm⁻¹]
        Z: effective ionization degree
    '''
    neff = Z * pl.sqrt(Cst.RYD/(100*Eion))
    #print('neff:', neff)
    return 1.0 - 0.1728*(1/neff)**0.333 * (2*(neff/n)**2-1)

def gaunt_karzas_and_latter(Z, n, l, Eion):
    ''' Gaunt factor from Karzas & Latter (1961)
        n: principal quantum number
        l: secundary quantum number
        Eion: Ionization energy [cm⁻¹]
        Z: effctive ionization degree
    '''
    Eion = 100*Cst.H*Cst.C/Cst.Q * Eion # eV
    #print('Eion:', Eion)
    return gaunt_bf(n, l, Z, Eion)

#def sigma0_bf(Z, n, gaunt):
#    ''' Expression from Mihalas (1978)
#        Return photoionisation cross-section at threshold in m²
#    '''
#    return 7.91e-22 * n * gaunt /Z**2

def sigma0_bf_kramers(n, Eion, gbf):
    ''' Semi-classical relation from Kramers 
        (modified version of Travis & Matsushima (1968) eq. 18 and 19)
        Return photoionisation cross-section at threshold in m²
        Eion: ionization energy of the level [cm⁻¹]
        Z: ionic charge
        n: principal quantum number
        gbf: bound-free Gaunt factor
    '''
    nu0 = Eion * 100 * Cst.C 
    #return 2.815e25 *  Z**4 /(n**5 * nu0**3)
    return 2.815e25 * (Eion*100/Cst.RYD)**2 *1/(n*nu0**3) * gbf

def id_uniq_qm_px(level, qm_levels, idx_table, CUTOFF, DE, OFILE4):
    ''' For a given level [Level|State]
        find the QM level and index of photoionisation table associated and return photoionization table.
    '''
    EV_TO_CM = Cst.Q / 100. / Cst.H / Cst.C

    # Identification by configuration and term
    qm_cfgterm_list = [qm_level.cfg+' '+qm_level.term for qm_level in qm_levels]
    res = dl.get_close_matches(level.cfg+' '+str(level.term), qm_cfgterm_list, n=1, cutoff=CUTOFF)
    try:
        idx = qm_cfgterm_list.index(res[0])
    except IndexError:
        return None, None, None

    # See if theoretical and experimental energy are not too much different
    qm_e = qm_levels[idx].e/EV_TO_CM
    exp_e = level.e/EV_TO_CM
    if abs(qm_e-exp_e) > DE:
        #print('Warning: energy drift >', DE, ' eV:', level.cfg, level.term, ' '.join(res), ' (i.e. ID not reliable)')
        return None, None, None 

    # Find indexes useful for direct access of the binary photoionization file
    try:
        i, nmin, ntot, m, l, p, pos = idx_table[idx] 
    except IndexError:
        return None, None, None

    # Case where level exist in TOPBASE but there is no photoionization table
    if int(ntot) == 0:
        return None, None, None 

    # Extract photoionization from binary file by direct access
    e_ryd0, x_Mb0 = extract_px(OFILE4, i, ntot)
    e_eV0 = e_ryd0 * 13.6057 # energy in eV

    # Exclude null values
    e_eV = pl.compress( x_Mb0 != 0., e_eV0)
    x_Mb = pl.compress( x_Mb0 != 0., x_Mb0)

    return idx, e_eV, x_Mb

def px_smooth(idx, e, x, idx_table, N_HE0, N_US, N_US_HE, WC):
    '''Over sample, smooth and undersample photoionization cross-sections
    '''
    i, nmin, ntot, m, l, p, pos = idx_table[idx]

    try: 
        #case of TOPBASE data
        nmin.index('.')
        nmin = pl.nan
    except ValueError:
        nmin = int(nmin)

    # Keep sampling for high energy values where the variation follow Kramer's law
    if isinstance(int(ntot)-nmin, int):
        N_HE = int(ntot)-nmin
    else:
        N_HE = N_HE0    

    if N_HE >= e.size:
        N_HE = -e.size
        print("Warning: N_HE is larger than photoionization table, select all the table.")

    e_sel = e[:-N_HE]
    e_sel_log = pl.log10(e_sel)
    x_sel = x[:-N_HE]        

    # Interpolate and smooth data
    #e_i = pl.linspace(min(e_sel), max(e_sel), 10000)
    e_i_log = pl.linspace(min(e_sel_log), max(e_sel_log), 10000)
    e_i = 10**e_i_log
    x_i = pl.interp(e_i, e_sel, x_sel)
    x_is = smooth(x_i, WC)      
    e_us = pl.concatenate([e_i[0:10], e_i[::N_US], e[int(ntot)-N_HE::N_US_HE]])
    x_us = pl.concatenate([x_is[0:10], x_is[::N_US], x[int(ntot)-N_HE::N_US_HE]])       

    if x_us.any() == 0.:
        print('x_us = 0')
        quit(1)

    # Conservation of area
    #area = pl.trapz( x_Mb, e_eV)   # total
    #area = pl.trapz( e_sel, x_sel) # selected
    area_i = pl.trapz(x_i, e_i)     # selected interpolated
    area_is = pl.trapz(x_is, e_i)   # selected interpolated and sampled
    #area_us = pl.trapz(x_us, e_us)

    return e_us, x_us, area_i, area_is

def smooth(x, wc):
   ''' Smooth data x by a third order Butterworth low-band filter characterized by a cut frequency wc [rad/s] 
   ''' 
   b,a = ss.butter(3, wc)
   return ss.filtfilt(b, a, x)

def plot_px(ef, xf, a, a_us, level, qm_levels, idx_table, dumb, N_HE0, MPI):
    ''' 
        ef: final smoothed and undersampled energies (in eV)
        xf: final smoothed and undersampled cross-sections (in Mb)
        a: area of the original photoionization
        a_us: area of the smoothed and undersampled photoionization
        dumb: same as previous but as raw data from NORAD/TOPBASE 
    '''
    global iplot 

    iplot += 1
    ofname = format(iplot, '0>4')+'_px.pdf'

    pl.figure(figsize=(12, 6), dpi=100, facecolor='w', edgecolor='k')
    pl.xlabel('E [eV]')
    #pl.xlim([0.1, 200])
    pl.ylim([1e-5, 1e5])
    pl.ylabel('Cross section [Mb]')
    pl.semilogy()
    pl.semilogx()
    title1 = SPE+' '+DEG+': '+str(level.cfg)+' '+str(level.term)+' [g='+str(level.g)+'] \n'
    title2 = QMP.upper()+': '
    #title1 = str(level.cfg)+' '+str(level.term)+' [g='+str(level.g)+'] \n'
    #title2 = "TOPBASE: "

    g_theo = 0
    nb = 0 # Number of photoionization tables on the plot 

    for i, e, x in dumb:

        # Missing photoionization tables
        if i == None:
            continue

        title2 += qm_levels[i].cfg+' '+qm_levels[i].term+' [g='+str(qm_levels[i].g)+'] '
        
        if not (nb+1)%5 and not (nb+1) == len(dumb):
            title2 += '\n' 
        
        g_theo += qm_levels[i].g
        
        idx, nmin, ntot, m, l, p, pos = idx_table[i]

        try: 
            #case of TOPBASE data
            nmin.index('.')
            nmin = pl.nan
        except ValueError:
            nmin = int(nmin)

        # Keep sampling for high energy values where the variation follow Kramer's law
        if isinstance(int(ntot)-nmin, int):
            N_HE = int(ntot)-nmin
        else:
            N_HE = N_HE0

        if N_HE >= e.size:
           N_HE = -e.size
           print("Warning: N_HE is larger than photoionization table, select all the table.")    

        if len(dumb) <= 10:
            scale = nb/10.
        elif len(dumb) <= 20:
            scale = nb/20.
        else:
            scale = nb/50.

        if scale > 1.:
            scale = 1.

        pl.plot(e, x, color=str(scale), label='area ='+format(a[nb], '11.3e')+' ('+format(len(x))+')')
        e_shadow = pl.concatenate([pl.array([e[0]]), e[:-N_HE], pl.array([e[-N_HE-1]])]) 
        x_shadow = pl.concatenate([pl.array([1e-5]), x[:-N_HE], pl.array([1e-5])])
        pl.fill(e_shadow, x_shadow, alpha=0.2)

        nb += 1

    if MPI:
        g = g_theo # if missing photoionization tables
    else:
        g = level.g # If all photoionization tables exist

    if g == g_theo:
        pl.plot(ef, xf, 'r+-', label='area ='+format(sum(a_us), '11.3e')+' ('+format(len(xf))+')')
    else:
        pl.plot(ef, xf, 'r+-', label='area ='+format(level.g)+' / '+format(g_theo)+format(sum(a_us), '11.3e')+' ('+format(len(xf))+')')

    pl.title(title1+title2, fontsize=9)

    pl.legend(frameon=False, fontsize=9)
    pl.savefig(ofname, dpi=100, format='pdf', orientation='landscape', bbox_inches='tight', papertype='a4')
    pl.close()

def px_combine(i_list, e_list, x_list, level, qm_levels):
    ''' Combine differents photoionization cross-sections if necessary and adjust the statistical weight
        i_list: index list of qm_levels for the level considered
        e_list: smoothed and undersampled photoionization energies 
        x_list: smoothed and undersampled photoionization cross-sections
    '''

    if len(i_list) == 1:
        g_theo = qm_levels[i_list[0]].g

        if MPI:
            g = g_theo
        else: 
            g = level.g 

        e = e_list[0]
        x = g / g_theo * x_list[0]

    else:
        
        e_set = set(pl.flatten(e_list))

        eilg = pl.linspace(min(pl.log10(list(e_set))), max(pl.log10(list(e_set))), 50000)
        ei = 10**eilg
        xi_list = []
        for i, x in enumerate(x_list):
            xi = pl.interp(ei, e_list[i], x)
            xi_list.append(xi[::N_US])

        g_theo = pl.sum([qm_levels[i].g for i in i_list])

        if MPI:
            g = g_theo
        else: 
            g = level.g

        e = ei[::N_US]
        x = g / g_theo * pl.sum([xi for xi in xi_list], axis=0)

    return e, x


if __name__ == '__main__':

    IPATH = '/home/tmerle/development/formato2/'
    IFILE1 = 'ael'
    IFILE21 = IPATH+'ad/ael/'
    IFILE22 = IPATH+'ad/rt/bf/'
    OFILE0 = 'ael'
    OFILE1 = 'mart_bf_gaunt'
    OFILE2 = 'mart_bf_qmp_ael'
    OFILE3 = IPATH+'ad/rt/bf/'
    OFILE4 = IPATH+'ad/rt/bf/'
    QM_DATA_PATH = IPATH+'ad/rt/bf/'
    IL = True # Ionization Level present in levels
    MPI = False # Missing PhotoIonization table
    
    DESCRIPTION = 'Tool to merge atomic bound-free radiative transitions.'
    EPILOG = '2014-04-01 ThiM'

    EV_TO_CM = Cst.Q / 100. / Cst.H / Cst.C

    parser = ap.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
    parser.add_argument('species', type=str, default=None, help='Symbol of the atomic species (e.g. "Fe" or "FE" or "fe")')
    parser.add_argument('degree', type=str, default=None, help='Ionization degree (e.g. "I" or "i" for neutral species)')
    parser.add_argument('-i', '--input', type=str, default=None, help='Atomic energy levels file (in output format of mael.py - default ael_xxx.bin)')
    parser.add_argument('-qmp', '--quantum_mechanical_photoionization', type=str, default=None, help='Include TIP/TOPBASE photionisation cross-sections if available [topbase|norad]')
    parser.add_argument('-de', '--delta_e_eV', type=float, default=0.9, help='Useful only if -qmp is used, absolute energy drift tolerance between exp and theoretical values (in eV) for identification')
    parser.add_argument('-c', '--cutoff', type=float, default=0.44, help='Useful only if -qmp is used, cutoff [0, 1] for the identification of configuration and term, default=0.44')
    parser.add_argument('-nus', '--n_under_sampling', type=int, default=150, help='Useful only if -qmp is used, step for undersampling of qm data, default=150')
    parser.add_argument('-n_he', '--n_high_energy', type=int, default=50, help='Useful only if -qmp is used, number of points at high energy keeped for no smoothing, default=50 (used only if not given in NORAD/TOPBASE data).')
    parser.add_argument('-nus_he', '--n_under_sampling_at_high_energy', type=int, default=5, help='Useful only if -qmp is used, step for undersampling of qm data at high energy')
    parser.add_argument('-f', '--frequency_cut', type=float, default=0.01, help='Useful only if -qmp is used, frequency cut of the smoothing function base on a third order Butterworth low pass filter')
    parser.add_argument('-l', '--energy_limit', type=float, default=800., help='Useful only if -qmp is used, wavelenght below which photoionization cross-sections are not considered (in Angstrom), default=911')
    parser.add_argument('-scp', '--semi_classical_photoionization', default=False, action='store_true', help='Use Kramers formula with Gaunt factor')
    parser.add_argument('-p', '--plot', default=False, action='store_true', help='Plot the photoionization cross-sections')

    arguments = parser.parse_args()

    SPE = arguments.species
    DEG = arguments.degree
    INP = arguments.input # Level file
    QMP = arguments.quantum_mechanical_photoionization
    DE = arguments.delta_e_eV
    CUTOFF = arguments.cutoff
    N_US = arguments.n_under_sampling
    N_HE0 = arguments.n_high_energy
    N_US_HE = arguments.n_under_sampling_at_high_energy
    WC = arguments.frequency_cut
    LBD_LIM = arguments.energy_limit
    SCP = arguments.semi_classical_photoionization
    PLOT = arguments.plot

    print('__________________')
    print('|   mart_bf.py   |')
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

    if QMP:
        if QMP == 'norad':
            IFILE21 = IFILE21+'norad_'+SYM+'.dat'
            IFILE22 = IFILE22+'norad_'+SYM+'_px.dat'
        elif QMP == 'topbase':
            IFILE21 = IFILE21+'topbase_'+SYM+'.dat'
            IFILE22 = IFILE22+'topbase_'+SYM+'_px.dat'            

    # Check the existence of the given line command argument's files
    if not os.path.isfile(IFILE1): 
        print("\nFile", IFILE1, "does NOT exist.")
        quit(1)

    if QMP:
        if not os.path.isfile(IFILE21) or not os.path.isfile(IFILE22):
            print("\nWARNING:")
            print("File", IFILE21,"and/or", IFILE22, "does NOT exist.")
            print("\nContinue without quantum mechanical bf data.\n")
            QMP = False

    # Read input levels data
    ifile = open(IFILE1, 'rb')
    levels = pickle.load(ifile)
    ifile.close()

    # Remove ionization level
    if ION in levels:
            levels.remove(ION)
    else:
        print("Warning: ionization State of", SPE, DEG, "not found in", SPE.strip())
        IL = False # Not ionization level present in levels

    print("Number of input levels for which bf transitions are requested:", len(levels))

# Semi Classical Photoionization using Gaunt factor
    if SCP:

        print('\nSemi-Classical Photoinization cross-section using Gaunt factor...')

        gaunt_mp_list = []
        gaunt_kl_list = []
        sigma0_mp_list = []
        sigma0_kl_list = []
        sigma0_kramers_list = []
        tot = []    

        for level in levels:
            #level.print()
            #print('n:', level.get_pqn(), 'l:', level.get_sqn())
            #print('Eion:', ION.e-level.e,'cm⁻¹')
            Z = level.i
            Eion = ION.e-level.e
            pqn = level.get_pqn()
            if pqn > 1000:
                pqn = Z*pl.sqrt(Cst.RYD*100/level.e)
                pqn = int(pqn)
            sqn = level.get_sqn()
            if pqn:
                gaunt_mp = gaunt_menzel_and_pickeris(Z, pqn, Eion)
                sigma0_mp = sigma0_bf_kramers(pqn, Eion, gaunt_mp)*1e22
                sigma0_kramers = sigma0_bf_kramers(pqn, Eion, 1.0)*1e22
                if sqn != None:
                    gaunt_kl = gaunt_karzas_and_latter(Z, pqn, sqn, Eion)
                    sigma0_kl = sigma0_bf_kramers(pqn, Eion, gaunt_kl)*1e22
                    level.sbf = sigma0_kl*1e-18
                else:
                    gaunt_kl = pl.nan
                    sigma0_kl = pl.nan
                    level.sbf = sigma0_mp*1e-18
            else:
                gaunt_mp = pl.nan
                sigma0_mp = pl.nan
                gaunt_kl = pl.nan
                sigma0_kl = pl.nan                 
                sigma0_kramers = pl.nan

            gaunt_mp_list.append(gaunt_mp)
            gaunt_kl_list.append(gaunt_kl)
            sigma0_mp_list.append(sigma0_mp)
            sigma0_kl_list.append(sigma0_kl)
            sigma0_kramers_list.append(sigma0_kramers)
            
            level.__format__()
            tot.append((level.e_p, level.cfg, Eion.__format__('12.3f'), Z, pqn, sqn, gaunt_mp.__format__('8.2e'), gaunt_kl.__format__('8.2e'), sigma0_mp.__format__('8.2e'), sigma0_kl.__format__('8.2e'), sigma0_kramers.__format__('8.2e')))
            #print('gaunt_mp:', gaunt_mp.__format__('10.3e'), 'sigma0_mp:', sigma0_mp.__format__('10.3e'), 'Mb')
            #print('gaunt_kl:', gaunt_kl.__format__('10.3e'), 'sigma0_kl:', sigma0_kl.__format__('10.3e'), 'Mb')    

        OFILE1 += '_'+SYM+'.dat'
        ofile = open(OFILE1, 'w')
        print('#     E[cm⁻¹]                       cfg    Eion[cm⁻¹]  Z    n    l  gbf_mp   gbf_kl sbf_mp[Mb] sbf_kl[Mb] sbf_kramers[Mb]', file=ofile)
        pl.savetxt(ofile, tot, fmt='%13s %25s %13s %2s %4s %4s %8s %8s %8s %8s %8s')
        ofile.close()

        print('Output ascii file:', OFILE1)

# Quantum Mechanical Photoionization cross section from TOPBASE/NORAD  
    if QMP:
        print('\nQuantum Mechanical Photoionization cross-section using '+QMP.upper()+' data...')

        # naming intermediate files
        if QMP == 'norad':
            OFILE3 = OFILE3+'norad_'+SYM+'_idx.dat'
            OFILE4 = OFILE4+'norad_'+SYM+'_px.bin'
        elif QMP == 'topbase':
            OFILE3 = OFILE3+'topbase_'+SYM+'_idx.dat'
            OFILE4 = OFILE4+'topbase_'+SYM+'_px.bin'    

        # format photoionization tables of all qm levels from ascii to binary with generation of an index table
        if not os.path.isfile(OFILE3) or not os.path.isfile(OFILE4):
            print(" Format photoionization table in binary format for saving time...")
            format_px(QMP, IFILE22, OFILE3, OFILE4)
        else:
            print(" Index and binary photoionization tables already computed. Skip this part.")


        # Read NORAD/TOPBASE ael data for mean levels (i.e. without fine structure)
        try:
            dumb = pl.loadtxt(IFILE21, dtype='str')
        except ValueError:
            print("Problem reading input file:", IFILE21)
            quit(1)       
            
        if dumb.shape[1] == 4:
            qm_levels = [Level(e=e, g=g, cfg=cfg, term=term, ref=QMP) for e, g, cfg, term in dumb]
        elif dumb.shape[1] == 6:
            qm_levels = [Level(e=e, g=g, cfg=cfg, term=term, p=p, ref=ref) for e, g, cfg, term, p, ref in dumb]
        else:
            print("Input format problem.")
            quit(1)     

        # open index table
        try:
            idx_table = pl.loadtxt(OFILE3, dtype='str')
        except ValueError:
            raise ValueError('Problem reading file '+OFILE3)

        OFILE2 += '_'+SYM+'.dat'
        ofile = open(OFILE2, 'w')
        print('# cfg                     term     g qm_cfg                    qm_term  qm_g e[eV] qm_e[eV]', file=ofile)

        nqmp = 0 # Number of levels for which QM data exists

        for level in levels:
            MPI = False
            dumb = []

            # Identification of index, photoionization energy and cross-sections in qm_levels and idx_table
            if isinstance(level, HyperLevel):
                hl = level
                if hasattr(hl, 'superlevels'):
                    for sl in hl.superlevels:
                        if hasattr(sl, 'levels'):
                            for lev in sl.levels:
                                i, e, x = id_uniq_qm_px(lev, qm_levels, idx_table, CUTOFF, DE, OFILE4)
                                if not i:
                                    MPI = True
                                dumb.append((i, e, x))
                        else:
                            i, e, x = id_uniq_qm_px(level, qm_levels, idx_table, CUTOFF, DE, OFILE4)
                            dumb.append((i, e, x))
                else:
                    i, e, x = id_uniq_qm_px(level, qm_levels, idx_table, CUTOFF, DE, OFILE4)
                    dumb.append((i, e, x))
                    
            elif isinstance(level, SuperLevel):
                sl = level
                if hasattr(sl, 'levels'):
                    for lev in sl.levels:
                        i, e, x = id_uniq_qm_px(lev, qm_levels, idx_table, CUTOFF, DE, OFILE4)
                        if not i:
                            MPI = True # Missing Photoionization table
                        dumb.append((i, e, x))                       
                else:
                    i, e, x = id_uniq_qm_px(level, qm_levels, idx_table, CUTOFF, DE, OFILE4)
                    dumb.append((i, e, x))
                    
            elif isinstance(level, Level) or isinstance(level, State):
                i, e, x = id_uniq_qm_px(level, qm_levels, idx_table, CUTOFF, DE, OFILE4)
                dumb.append((i, e, x))
                
            else:
                print('level is not an instance of classes [HyperLevel|SuperLevel|Level|State].')
                quit(1)

            i_list, e_list, x_list = [], [], []
            area_list, area_is_list = [], []
            qm_cfg_list, qm_term_list, qm_e_list, qm_g_list = [], [], [], []

            # Smooth data (possible several tables) for the selected level
            for i, e, x in dumb:
                if i != None:
                    try:
                        e_us, x_us, area_i, area_is = px_smooth(i, e, x, idx_table, N_HE0, N_US, N_US_HE, WC)
                    except:
                        print(qm_levels[i].print())
                        print(i, e, x)
                        quit(1)
                    i_list.append(i)
                    e_list.append(e_us)
                    x_list.append(x_us)
                    area_list.append(area_i)
                    area_is_list.append(area_is)
                    qm_cfg_list.append(qm_levels[i].cfg)
                    qm_term_list.append(qm_levels[i].term)
                    qm_e_list.append(qm_levels[i].e)
                    qm_g_list.append(qm_levels[i].g)
                else:
                    print("Photoionization not found for one of the (sub)level:", end=' ')
                    level.print()

            if i_list:
                nqmp += 1
                # Possible combination of photoionization for HyperLevel|SuperLevel or fraction for State
                ef, xf = px_combine(i_list, e_list, x_list, level, qm_levels)
                # Remove negative cross-sections
                ef = pl.compress( xf > 0., ef)
                xf = pl.compress( xf > 0., xf)
                
                # Compute the wavelength at the threshold and compare with limit LBD_LIM fixed by the user 
                e_thres = Cst.H * Cst.C * 1.e10 / Cst.Q / ef[0]
                x_thres = xf[0] * 1.e-18 
                if e_thres < LBD_LIM:
                    print("Warning: ionization energy corresponds to wavelenght lower than the limit:", LBD_LIM, ' Å')   

                # Remove cross-section for energies lower than LBD_LIM A
                ef = pl.compress( Cst.H * Cst.C * 1.e10 / Cst.Q / ef >= LBD_LIM, ef)
                xf = pl.compress( Cst.H * Cst.C * 1.e10 / Cst.Q / ef >= LBD_LIM, xf)
                
                # Plot photoionization cross-section for the selected level
                if PLOT:
                    plot_px(ef, xf, area_list, area_is_list, level, qm_levels, idx_table, dumb, N_HE0, MPI)

                qm_e = pl.vdot(qm_g_list, qm_e_list)/sum(qm_g_list)

                print(format(level.cfg, '25'),
                      format(level.term, '5'),
                      format(level.g, '4'),
                      format(' '.join(qm_cfg_list), '25'),
                      format(' '.join(qm_term_list), '5'),
                      format(sum(qm_g_list), '4'), 
                      format(level.e/EV_TO_CM, '5.3f'), 
                      format(qm_e/EV_TO_CM, '5.3f'), file=ofile, end='')

                de = abs(level.e-qm_e)/EV_TO_CM
                if de > DE:
                    print(format(de, '6.3f'), file=ofile)
                else:
                    print(file=ofile)
            else:
                # No photoionisation found
                print(format(level.cfg, '25'),
                      format(level.term, '5'),
                      format(level.g, '4'), file=ofile)

            # Add the (combined,) smoothed and undersampled photoionisation table to the level instance
            if i_list:
                if e_thres > LBD_LIM:
                    pe = Cst.H * Cst.C * 1.e10 / Cst.Q / ef  # eV => Angstrom
                    px = xf * 1.e-18                         # Mb => cm²
                    pe_px_list = sorted(set(zip(pe, px)), reverse=True)
                    pe = [e for e, x in pe_px_list]
                    px = [x for e, x in pe_px_list]
                    level.pe = pe
                    level.px = px
                    level.sbf = level.px[0]
                else:
                    level.sbf = x_thres

        print('\n Number of input levels for which QM data are available:', nqmp)
                #print(' Among them, number of levels with energy drift >', DE, ' eV:', list(crit).count(True), ' (i.e. identification not reliable)')   
        
        ofile.close()    
        print(' Output ascii file:', OFILE2)    

#####

#        #qm_levels, idx_qm_levels = qm_id_ael(levels, IFILE21, OFILE2)#

#        # Select levels for which one wants photoionization table
#        idx_qm_levels = [i for i in idx_qm_levels if i < len(idata)]
#        idata = idata[idx_qm_levels] #

#        print( ' Number of levels with truly photoionisation tables available:', len(idata))#
#

#        for idx, (i, j, k, l, m, n, o) in enumerate(idata):#

#            try: 
#                #case of TOPBASE data
#                j.index('.')
#                j = pl.nan
#            except ValueError:
#                j = int(j)#

#            print(format(i,'>8'), format(j, '>5'), format(k, '>5'), l, m, n, format(o, '>2'), end=' ')#

#            # Extract photoionization from binary file by direct access
#            e_ryd0, x_Mb0 = extract_px(OFILE4, i, k)
#            e_eV0 = e_ryd0 * 13.6 # energy in eV#

#            # Exclude null values
#            e_eV = pl.compress( x_Mb0 != 0., e_eV0)
#            x_Mb = pl.compress( x_Mb0 != 0., x_Mb0)#

#            # Keep sampling for high energy values where the variation follow Kramer's law
#            if isinstance(int(k)-j, int):
#                N_HE = int(k)-j#

#            e_sel = e_eV[:-N_HE]
#            x_sel = x_Mb[:-N_HE]#

#            # Interpolate and smooth data
#            e_i = pl.linspace(min(e_sel), max(e_sel), 10000)
#            x_i = pl.interp(e_i, e_sel, x_sel)
#            x_is = smooth(x_i, WC)#
#

#            e_us = pl.concatenate([e_i[0:10], e_i[::N_US], e_eV[int(k)-N_HE::N_US_HE]])
#            x_us = pl.concatenate([x_is[0:10], x_is[::N_US], x_Mb[int(k)-N_HE::N_US_HE]])#

#            if x_us.any() == 0.:
#                print('x_us = 0')
#                quit(1)#

#            # Conservation of area#

#            #area = pl.trapz( x_Mb, e_eV)   # total
#            #area = pl.trapz( e_sel, x_sel) # selected
#            area_i = pl.trapz(x_i, e_i)     # selected interpolated
#            area_is = pl.trapz(x_is, e_i)   # selected interpolated and sampled
#            #area_us = pl.trapz(x_us, e_us)#

#            #n_area = len(x_Mb)
#            #n_area_us = len(x_us)#

#            #area_er = abs(area - area_us)/ area
#            area_er = abs(area_i - area_is)/ area_i#

#            # For fine level take a fraction of the photoionization table
#            print(levels[idx].g, qm_levels[idx_qm_levels[idx]].g,  levels[idx].cfg, qm_levels[idx_qm_levels[idx]].cfg, end=' ')
#            x_us = levels[idx].g / qm_levels[idx_qm_levels[idx]].g * x_us#

#            if area_er < 0.05:
#                #print(format(area, '10.3e'), format(area_i, '10.3e'), format(area_is, '10.3e'), format(area_us, '10.3e'),  format(len(x_us)), format(area_er*100,'3.1f'), end=' %')
#                print(format(area_i, '10.3e'), format(area_is, '10.3e'),  format(len(x_us)), format(area_er*100,'3.1f'), end=' %')
#            else:
#                #print(format(area, '10.3e'), format(area_i, '10.3e'), format(area_is, '10.3e'), format(area_us, '10.3e'),  format(len(x_us)), format(area_er*100,'3.1f'), end=' % WARNING: bad undersampling...')
#                print(format(area_i, '10.3e'), format(area_is, '10.3e'),  format(len(x_us)), format(area_er*100,'3.1f'), end=' % WARNING: bad undersampling...')#

#            if len(x_Mb0) != len(x_Mb):
#                print('  ', len(x_Mb0)-len(x_Mb), " null values.")
#            else:
#                print('')#

#            #plot_px_old(e_eV, x_Mb, e_us, x_us, area, area_us, (i, j, k, l, m, n, o))
#            plot_px_old(e_eV, x_Mb, e_us, x_us, area_i, area_is, (i, j, k, l, m, n, o))

#####


# Internal output for reading by other python code
    if IL:
        levels.append(ION)

    OFILE0 = OFILE0+'_'+SYM+'.bin'
    ofile = open(OFILE0, 'wb')
    pickle.dump(levels, ofile, protocol=1)
    ofile.close()

    print('\nOutput binary file:', OFILE0)


