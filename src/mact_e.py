#!/usr/bin/python
# -*- coding: utf8 -*-
''' Code to merge atomic electron collisions of a given species in a given ionization degree.
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
import scipy.special as ss 
import scipy.optimize as so
import scipy.interpolate as si

try: 
    import pylab as pl 
except ImportError: 
    raise ImportError('pylab module not installed')

from mad import Cst, ions, State, Collision, orbital_radius

iplot = 0 # For numeroation of cross-section plots

cst = Cst.H * Cst.C / (Cst.Q * 1e-10)

def ups_sce_fisher(gi, f, Eo, dn):
    ''' Compute the effective electron collision strength
        with semi-classical formula from Fisher et al. 1996, Phys. Rev. A
        Eo [eV] transition energy
        8pi/sqrt(3)* 13.6057 = 197.42
    '''
    cst = 197.42/Eo * gi * f
    def func(T):
        return cst * average_gaunt_fisher(dn, Eo, T)
    return func

def average_gaunt_fisher(dn, Eo, T):
    ''' Compute the maxwellian average of the Gaunt factor for electron collision
        Only for neutral species
    '''
    x_vec = pl.linspace(0., 8., 50) # sqrt(E/Eo)
    x0 = Eo * Cst.Q / (Cst.K * T)    # Eo/kT
    g_func = make_gaunt_fisher(dn)
    g_vec = [g_func(x) * pl.exp(-x**2*x0) for x in x_vec]
    return pl.trapz(g_vec, x_vec**2*x0)

def make_gaunt_fisher(dn):
    ''' Implementation of the semi-classical formula of Fisher et al. 1996, Phys. Rev.A
        Only for neutral species
        dn is the difference of the quantum principal number in the transition
    ''' 
    if dn == 0:
        def f(x): # x=sqrt(E/Eo)
            x = x + 1 # change incident energy to energy after collision
            return (0.33-0.3/x**2 +0.08/x**4) * pl.log(x**2)
    else:
        def f(x): # x=sqrt(E/Eo)
            x = x + 1 # change incident energy to energy after collision
            return (0.276-0.18/x**2) * pl.log(x**2)
    return f 

def ups_sce_vanregemorter(gi, f, Eo, deg):
    ''' Compute the effectice electron collision strength
        with semi-classical formula from Van Regemorter 1962, ApJ
        Eo [eV] transition energy
    '''    
    cst = 197.42/Eo * gi * f
    def func(T):
        #return cst * average_gaunt_vanregemorter(deg, Eo, T)
        return cst * average_gaunt_vanregemorter2(deg, Eo, T)
    return func

def average_gaunt_seaton(ups, line):
    ''' Compute the maxwellian average of the Gaunt factor for electron collision
        from the Upsilon (only for comparison purpose)
    '''
    global cst
    cst_ups = 197.42/(cst/line.lbd) * line.lower.g * line.f
    return ups/cst_ups

def average_gaunt_vanregemorter2(deg, Eo, T):
    ''' Compute the maxwellian average of the Gaunt factor for electron collision
    '''
    x_vec = pl.linspace(0.2, 8., 50) # sqrt(E/Eo)
    x0 = Eo * Cst.Q / (Cst.K * T)    # Eo/kT
    g_func = make_gaunt_vanregemorter(deg)
    g_vec = [g_func(x) * pl.exp(-x**2*x0) for x in x_vec]

    return pl.trapz(g_vec, x_vec**2*x0)

def make_gaunt_vanregemorter(deg):
    ''' according to Van Regemorter 1962, ApJ, Table 1. and Fig. 1
    '''
    x = [0.0, 0.2, 0.4, 0.6, 0.8, 1., 2., 3., 4., 5., 6., 7., 8.]# sqrt(E/Eo)
    #x = [0.04, 0.16, 0.36, 0.64, 1., 4., 9., 16., 25., 36.] # E/Eo
    gi = [0.0, 0.015, 0.034, 0.057, 0.084, 0.124, 0.328, 0.561, 0.775, 0.922, 1.040, 1.12, 1.20]
    gii = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.328, 0.561, 0.775, 0.922, 1.040, 1.12, 1.20]

    if deg == 'I':
        return si.interp1d(x, gi, kind='cubic')
    elif deg == 'II':
        return si.interp1d(x, gii, kind='cubic')
    else:
        print("Z must be 'I' or 'II'")
        quit(1)

def average_gaunt_vanregemorter(deg, Eo, T):
    ''' Interpolate the maxwellian average of the Gaunt factor for electron collisions
        according to Van Regemorter 1962, ApJ
        deg: ionization degree
        Eo: transition energy [eV]
        T: temperature [K]
    '''
    x = [0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1., 2., 4., 10.]
    av_gi = [1.160, 0.956, 0.758, 0.493, 0.331, 0.209, 0.100, 0.063, 0.040, 0.023]
    av_gii = [1.160, 0.977, 0.788, 0.554, 0.403, 0.290, 0.214, 0.201, 0.200, 0.200]

    x0 = Eo * Cst.Q / (Cst.K * T)  

    if x0 <= 0.005:
        return 0.276 * (-0.57722 - pl.log(x0)) 

    if deg.upper().strip() == 'I':
        if x0 > 10.:
            return 0.066 * x0**-0.5
        else:
            return pl.interp(x0, x, av_gi)
    elif deg.upper().strip() == 'II':
        if x0 > 10.:
            return 0.200
        else:
            return pl.interp(x0, x, av_gii)
    else:
        print("Z must be 'I' or 'II'")
        quit(1)

def phi(x):
    ''' Function for the IPM (Seaton, 1962)
    '''
    return x * ss.k0(x) * ss.k1(x)

def zeta(x):
    ''' Function for the IPM (Seaton, 1962)
    '''
    return x**2 * (ss.k0(x)**2 + ss.k1(x)**2)

def beta0(x, x0, r, T):
    ''' argument of phi function for weak coupling
        line: Line of Multiplet instance
        x: energy of colliding electron before collision [without dimension]
        x0: energy of the collision transition [without dimension]
        r: hydrogenoid orbital radius
        T: temperature [K]
    '''
    #return r * pl.sqrt(Cst.K * T * (x+x0)/13.6057/Cst.Q) * x0 / (2*x + x0)
    return r * pl.sqrt(Cst.K * T * (x)/13.6057/Cst.Q) * x0 / (2*x - x0)

def func_seaton(b1, c1):
    '''Implicit equation to solve for Seaton IPM method
    '''
    return ss.k0(b1)**2 + ss.k1(b1)**2 - c1

def beta1(x, x0, f, T):
    ''' argument of the implicit equation for the strong coupling
        x: energy of colliding electron before collision [without dimension]
	    x0: energy of the collision transition [without dimension]
	    f: oscillator strength
	    T: temperature [K]
    '''
    #c1 = Cst.K * T * x0 / (8. * Cst.RYD * Cst.H * Cst.C * f) * ((2*x+x0)/x0)**2
    c1 = Cst.K * T * x0 / (8. * Cst.RYD * Cst.H * Cst.C * f) * ((2*x-x0)/x0)**2
    b1 = so.brentq(func_seaton, 0, 1000, args=c1)
    return b1

def sigma_weak_seaton(x, ion, line, T):
    ''' Compute collision cross-section for weak coupling in the IPM formula (Seaton 1962)
        x: energy of colliding electron before collision [without dimension]
        ion: Ionization instance of State of Level class
	    line: radiative transtion instance of Line class
	    T: temperature [K]
    '''
    global cst
    
    r_i = orbital_radius(ion, line.lower)
    r_j = orbital_radius(ion, line.upper)
    #print(r_i, r_j)
    
    if r_i <= r_j:
        r = r_i
    else:
        r = r_j   
    
    de = cst / line.lbd
    x0  = de * Cst.Q / Cst.K / T
    #print(r, de, x0)
    b0 = beta0(x, x0, r, T)
    if x < x0:
        return 0.
    else:
        #return 8 * (13.6057*Cst.Q/Cst.K/T)**2 * line.f * (1/x0) * (1/(x+x0)) * phi(b0)
        return 8 * (13.6057*Cst.Q/Cst.K/T)**2 * line.f * (1/x0) * (1/(x)) * phi(b0)

def sigma_strong_seaton(x, ion, line, T):
    ''' Compute collision cross-section for strong coupling in the IPM formula (Seaton 1962)
        x: energy of colliding electron before collision [without dimension]
	    ion: Ionization instance of State of Level class
	    line: radiative transtion instance of Line class
	    T: temperature [K]
    '''
    global cst
    de = cst / line.lbd
    x0  = de * Cst.Q / Cst.K / T    
    b1 = beta1(x, x0, line.f, T)
    if x < x0:
        return 0.
    else:
        #return 8 * (13.6057*Cst.Q/Cst.K/T)**2 * line.f * (1/x0) * (1/(x+x0)) * (0.5*zeta(b1) + phi(b1))
        return 8 * (13.6057*Cst.Q/Cst.K/T)**2 * line.f * (1/x0) * (1/(x)) * (0.5*zeta(b1) + phi(b1))

def sigma_seaton(x, ion, line, T):
    ''' Select the lowest cross-section from the weak and strong coupling
    '''
    q0 = sigma_weak_seaton(x, ion, line, T)
    q1 = sigma_strong_seaton(x, ion, line, T)

    if q0 < q1:
        return q0
    else:
        return q1

def ups_sce_seaton(ion, line):
    ''' Integrate from the IPM cross-sections
    '''
    global cst
    de = cst / line.lbd
    def func(T):
        x0  = de * Cst.Q / Cst.K / T
        #x_vec = pl.linspace(0.0001, 100, 300)
        x_vec_log = pl.linspace(-5, 4, 100)
        x_vec = 10**x_vec_log
        q_vec = [sigma_seaton(i+x0, ion, line, T) for i in x_vec]
        ups_cst = line.lower.g * Cst.K * T / Cst.RYD / Cst.H / Cst.C 
        return ups_cst * pl.trapz(q_vec*(x_vec+x0)*pl.exp(-(x_vec)), x_vec)
        #return ups_cst * pl.trapz(q_vec*(x_vec)*pl.exp(-(x_vec-x0)), x_vec)
    return func

def plot_ex(ion, line):
    ''' Plot the collisions cross-section with electrons
    '''
    global cst, iplot
    iplot += 1
    ofname = format(iplot, '0>4')+'_ex.pdf'
    T = 5000.
    x0  = cst / line.lbd * Cst.Q / Cst.K / T
    x_vec_log = pl.linspace(-3, 4, 100)
    x_vec = 10**x_vec_log
    q0 = [sigma_weak_seaton(i, ion, line, T) for i in x_vec]
    q1 = [sigma_strong_seaton(i, ion, line, T) for i in x_vec]
    q = [sigma_seaton(i, ion, line, T) for i in x_vec]
    ll = line.lower
    ul = line.upper
    #title='Na I: '+str(ll.cfg)+' '+str(ll.term)+' [g='+str(ll.g)+'] <=> '+str(ul.cfg)+' '+str(ul.term)+' [g='+str(ul.g)+'], Ro = '+format(orbital_radius(ion, line.lower), '4.2f')+' a$_0$, f = '+str(line.f)
    title = SPE+' '+DEG+': '+str(ll.cfg)+' '+str(ll.term)+' [g='+str(ll.g)+'] '+str(ul.cfg)+' '+str(ul.term)+' [g='+str(ul.g)+'] Ro = '+str(orbital_radius(ion, line.lower))+' f = '+str(line.f)
    pl.figure(figsize=(12, 6), dpi=100, facecolor='w', edgecolor='k')
    pl.plot(x_vec*Cst.K*T/Cst.Q, q0, 'k', lw=1, label='Weak coupling')
    pl.plot(x_vec*Cst.K*T/Cst.Q, q1, 'k', lw=2, label='Strong coupling')
    pl.plot(x_vec*Cst.K*T/Cst.Q, q, 'r+', label='IPM (Seaton 1962)')
    pl.semilogx()
    #pl.semilogy()
    pl.xlim([1e-2, 1e4])
    pl.xlabel('E [eV]')
    pl.ylabel('Cross-section [$\pi a_0^2$]')
    pl.legend(frameon=False, loc=0)
    pl.title(title)
    #bbox_inches='tight'
    pl.savefig(ofname, dpi=100, format='pdf', orientation='landscape', papertype='a4')
    pl.close()

def plot_ups(col):
    ''' Plot the electronic collision strength as a function of temperature
    '''
    global iplot
    t = col.t_list
    ups = col.ups_list
    ref = col.ref
    #Name of the output plot file
    ofname = format(iplot, '0>4')+'_ups_e.pdf'    
    #Figure
    pl.figure(dpi=100, facecolor='w', edgecolor='k')
    #Title of the figure
    ll = col.lower
    ul = col.upper
    title = SPE+' '+DEG+': '+str(ll.cfg)+' '+str(ll.term)+' [g='+str(ll.g)+'] '+str(ul.cfg)+' '+str(ul.term)+' [g='+str(ul.g)+'], $\lambda =$'+format(col.lbd, '8.1f')+' $\AA$'
    #Axes
    pl.xlabel('T [K]')
    pl.ylabel('$\Upsilon_e$')
    pl.semilogy()
    pl.ylim([1e-6, 1e5])
    pl.title(title)
    #Plot
    pl.plot(t, ups, 'k', label=col.ref+' ('+col.type.lower()+')')
    pl.legend()
    pl.savefig(ofname, dpi=100, format='pdf', orientation='landscape', papertype='a4')
    pl.close()

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

def load_act(sym, source):
    ''' Load atomic collision transition file from mechanical quantum computations
        3 ascii input files:
        ael.dat: atomic energy levels in the mael.py format
        temp.dat: list of temperatures
        ups.dat: effective collision strengths
    '''
    path = '/home/tmerle/development/formato2/ad/ct/e/qm/'+SYM+'/'+source+'/'

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

    for i, (l, u) in enumerate(zip(l_list, u_list)):
        #print(l, u, ups_list[i])
        col = Collision(type='QM', lower=ael_list[l-1], upper=ael_list[u-1], t_list=t_list, ups_list=ups_list[i], ref=source, kw='UPS_E')
        col_list.append(col)

    return col_list


if __name__ == '__main__':

    #t_list = [1000, 2000, 5000, 10000]
    t_list = [1000, 2000, 4000, 6000, 8000, 10000, 15000, 20000]

    EV_TO_CM = Cst.Q / 100. / Cst.H / Cst.C

    IFILE1 = 'ael'
    IFILE2 = 'art_bb'
    OFILE0 = 'act_e'
    OFILE1 = 'mact_e' 
    QM_DATA = ['zatsarinny2009', 'pelan1997', 'zhang1995']
    RT = True

    DESCRIPTION = 'Tool to merge atomic collision transition with electrons'
    EPILOG = '2014-03-20 ThiM'
    parser = ap.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
    parser.add_argument('species', type=str, default=None, help='SYMbol of the atomic species (e.g. "Fe" or "FE" or "fe")')
    parser.add_argument('degree', type=str, default=None, help='Ionization degree (e.g. "I" or "i" for neutral species)')
    #parser.add_argument('-i', '--input', default=False, action='store_true', help='Atomic energy levels file (in output format of oomael.py - default levels.bin)')
    parser.add_argument('-f1', '--filename1', default='', type=str, help='Filename of energy levels produced by mael.py')
    parser.add_argument('-f2', '--filename2', default='', type=str, help='Filename of radiative b-b transitions produced by mart_bb.py')
#parser.add_argument('-o', '--output', type=str, default=None, help='name of the ion for naming the output files')
    parser.add_argument('-lte', '--lte', default=False, action='store_true', help='Enforce collisions to ensure LTE')
    parser.add_argument('-sce', '--semi_classical_e', default=None, type=str, help='Choose semi-classical formula for electron collisions [Seaton | VanRegemorter | Fisher]')
    parser.add_argument('-flow', '--forbidden_line', default=1e-5, type=float, help='Useful only if -sce is used, effective collision strength of 1.0 for semi-forbidden lines with f lower than -flow')
    parser.add_argument('-qme', '--quantum_mechanical_e', default=None, type=str, help='Choose quantum mechanical data if available for electron collisions')
    parser.add_argument('-t', '--temperature', nargs='*', default=t_list, help='Temperature list')
    parser.add_argument('-de', '--delta_e', type=float, default=50, help='allowed energy gap [cm⁻¹] for identification with quantum mechanical data (useful when dealing with fine structure')
    parser.add_argument('-p', '--plot', default=False, action='store_true', help='Plot the effective collisions strengths (and cross-section in case of Seaton)')

    arguments = parser.parse_args()

    SPE = arguments.species
    DEG = arguments.degree
    IF1 = arguments.filename1
    IF2 = arguments.filename2
    #INP = arguments.input
    #ONAME = arguments.output
    lte = arguments.lte
    SCE = arguments.semi_classical_e
    FLOW = arguments.forbidden_line
    QME = arguments.quantum_mechanical_e
    t_list = map(float, arguments.temperature)
    DE = arguments.delta_e
    PLOT = arguments.plot

    print('__________________')
    print('|    mact_e.py   |')
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
    try:
        ifile = open(IFILE2, 'rb')
        lines_list = pickle.load(ifile)
        ifile.close()
    except:
        print("File", IFILE2, "does NOT exist.")
        RT=False

    nlev = len(levels_list)
    ncol = int((nlev-1) * ((nlev-1)-1)/2)
    if RT:
        nlin = len(lines_list)

    print('\nKind of e⁻ bb collisions:', QME,' ', SCE)
    print('Kind of e⁻ bf collisions:', 'Seaton')

    print('\nNumber of energy levels:           ', str(nlev).rjust(6))
    print('\nNumber of e⁻ collision transitions:', str(ncol).rjust(6))
    if RT:
        print(' with oscillator strengths:', str(nlin).rjust(6))
        print(' forbidden transitions:    ', str(ncol-nlin).rjust(6))

# Electron semi-classical collisions
    qm_col_list = []
    lu_qm_col_list = [] # list of the level indexes
    allow_col_list = []
    forbd_col_list = []

    if QME:
        # Be careful
        # Treatment implemented only for qm data relative to mean levels
        # Manage distribution over fine levels included in list of energy levels selected
        # Not yet implemented for the reverse i.e. for qm fine level data (excepted for pelan1997 which is given between fine levels)

        nqmcs = 0 # Number of qm collisions selected

        print("\nQuantum mechanical electron collisions:", QME)

        if QME in QM_DATA:
            init_qm_col_list = load_act(SYM, QME)

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
                            ups = [levels_list[i].g*levels_list[j].g/(gl*gu) * float(upse) for upse in init_qm_col.ups_list]
                            col = Collision(type=init_qm_col.type, lower=levels_list[i], upper=levels_list[j], t_list=init_qm_col.t_list, ups_list=ups, kw=init_qm_col.kw, ref=QME)
                            col.l = i+1
                            col.u = j+1
                            qm_col_list.append(col)
                            lu_qm_col_list.append((i+1, j+1))
                            iplot += 1
                            if PLOT:
                                plot_ups(col)
        else:
            print("QME:", QME, " not yet implemented.")
            print()

        qm_col_list = sorted(qm_col_list, key=op.attrgetter('l', 'u', 'lbd'))   

        print(" Number of initial collisions:            ", len (init_qm_col_list))
        print(" Number of initial unselected collisions: ", len(init_qm_col_list)-nqmcs)
        print(" Number of initial selected collisions:   ", nqmcs)
        print(" Number of initial splitted selected collisions:   ", len(qm_col_list))

    if SCE and RT:
        nsccs = 0 # Number of semi-classic collision selection
        nrfcs = 0  # Number of radiatively forbidden collision selection 
        print("\nSemi-classical electron collisions:", SCE)
        print(' Temperatures selected:', t_list)

        # Allowed transitions
        if SCE == 'Fisher':
            if DEG != 'I':
                print("Fisher formula only for neutral species.")
                quit(1)
            for line in lines_list:
                # Use SCE if not QM 
                if (line.l, line.u) not in lu_qm_col_list: 
                    Eo =  cst / line.lbd
                    gi = line.lower.g
                    f = line.f
                    dn = line.dn()
                    if f < FLOW:
                        ups = list(pl.ones(len(t_list)))
                        type_str = 'F LOWER THAN '+str(FLOW)
                    else:             
                        upse = ups_sce_fisher(gi, f, Eo, dn)
                        ups = [upse(float(t)) for t in t_list]
                        type_str= 'ALLOWED'
                    col = Collision(type=type_str, lower=line.lower, upper=line.upper, t_list=t_list, ups_list=ups, lbd=line.lbd, kw='UPS_E', ref=SCE)
                    col.l = line.l
                    col.u = line.u
                    allow_col_list.append(col)
                    nsccs += 1
                    iplot += 1
                    if PLOT:
                        plot_ups(col)

        elif SCE == 'VanRegemorter':
            for line in lines_list:
                # Use SCE if not QM 
                if (line.l, line.u) not in lu_qm_col_list: 
                    Eo =  cst / line.lbd
                    gi = line.lower.g
                    f = line.f
                    if f < FLOW:
                        ups = list(pl.ones(len(t_list)))
                        type_str = 'F LOWER THAN '+str(FLOW)
                    else:             
                        upse = ups_sce_vanregemorter(gi, f, Eo, DEG)
                        ups = [upse(float(t)) for t in t_list]
                        type_str= 'ALLOWED'
                    col = Collision(type=type_str, lower=line.lower, upper=line.upper, t_list=t_list, ups_list=ups, lbd=line.lbd, kw='UPS_E', ref=SCE)
                    col.l = line.l
                    col.u = line.u
                    allow_col_list.append(col)
                    nsccs += 1
                    iplot += 1
                    if PLOT:
                        plot_ups(col)

        elif SCE == 'Seaton':
            for line in lines_list:
                # Use SCE if not QM 
                if (line.l, line.u) not in lu_qm_col_list:                 
                    Eo =  cst / line.lbd
                    gi = line.lower.g
                    f = line.f
                    #print(line.lower.cfg, line.upper.cfg)
                    if PLOT:
                        plot_ex(ION, line)
                    if f < FLOW:
                        ups = list(pl.ones(len(t_list)))
                        type_str = 'F LOWER THAN '+str(FLOW)
                    else:            
                        upse = ups_sce_seaton(ION, line) 
                        ups = [upse(float(t)) for t in t_list]
                        type_str= 'ALLOWED'
                    col = Collision(type=type_str, lower=line.lower, upper=line.upper, t_list=t_list, ups_list=ups, lbd=line.lbd, kw='UPS_E', ref=SCE)
                    col.l = line.l
                    col.u = line.u
                    allow_col_list.append(col)
                    nsccs += 1
                    iplot += 1
                    if PLOT:
                        plot_ups(col)
        else:
            print("Semi-classical formulæ [Fisher|VanRegemorter|Seaton]")
            quit(1)
        # Forbidden transitions
        for i in range(1, nlev):
            # List of index (low, up) energy levels
            idx = [(line.l, line.u) for line in lines_list if line.l == i]
            for j in range(i+1, nlev):
                # Forbidden transition are not in SCE (idx) nor in QM (lu_qm_col_list)
                if (i, j) not in idx+lu_qm_col_list:
                    ups =  list(pl.ones(len(t_list)))
                    col = Collision(type='FORBIDDEN', lower=levels_list[i-1], upper=levels_list[j-1], t_list=t_list, ups_list=ups, kw='UPS_E', ref='DEFAULT')
                    col.l = i
                    col.u = j
                    forbd_col_list.append(col)
                    nrfcs += 1
                    iplot += 1
                    if PLOT:
                        plot_ups(col)

        print(' Number of collisions a la ', SCE, ':', nsccs)
        print(' Number of radiatively forbidden collisions:', nrfcs)


    # Merge QME and SCE
    # TODO


#    for allow_col in allow_col_list:
#        print(allow_col)

#    for forbd_col in forbd_col_list:
#        print(forbd_col)

# Electron bound-free collisions

    bf_col_list = []

    for i in range(1, nlev):
        # For the moment just copy the bf cross-section at the threshold in unit of cm⁻²
        # To improve 
        ups = levels_list[i-1].sbf
        col = Collision(type='BF', lower=levels_list[i-1], upper=levels_list[nlev-1], ups_list=[ups], kw='CI', ref='SEATON 1962')
        col.l = i
        col.u = nlev
        bf_col_list.append(col)

# Merge bb allowed, forbidden and bf electron collisions
    col_list = allow_col_list+forbd_col_list+qm_col_list
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
