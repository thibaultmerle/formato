#!/usr/bin/python
# -*- coding: utf-8 -*-
''' This module defines the constants, the classes and some input atomic data.
'''

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
#from __future__ import unicode_literals

import re

try:
    import pylab as pl
except ImportError:
    raise ImportError('pylab module is not installed.')

nan = pl.nan

class Cst:
    ''' C: celerity of light [m/s]
        Q: electric charge [C]
        H: Planck constant [Js]
        K: Boltzmann constant [J/K]
        RYD: Rydberg constant [m⁻1]
        EPS0: Vacuum permittivity [F/m]
        ME: Electron mass [kg]
        A0: Bohr radius [m]
        AMU: Atomic mass unit [kg]
    '''
    C = 2.99792458e8    # m/s
    Q = 1.602176565e-19 # C
    H = 6.62606957e-34  # Js
    K = 1.3807e-23      # J/K
    RYD = 10973732      # m⁻¹
    EPS0 = 8.854187817e-12 # F/m
    ME = 9.10938215e-31    # kg
    A0 = 5.2917720859e-11  # m
    AMU = 1.660538921e-27 # kg

    # Dictionnary of atomic masses exressed in unit of AMU
    MASS = {'H': 1.008, 'O': 15.994, 'NA': 22.990, 'MG': 24.305, 'AL':26.982, 'K': 39.098, 'CA': 40.078, 'FE': 55.845}
 
#-------------------------------------------------------------------------------
#class Species(object):
#    ''' atomic species with a given ionization stage
#    '''
#    def __init__(self, sym, deg):
#        self.sym = sym
#        self.deg = deg
#        self.mass = Cst.MASS[sym.upper()] 

class State():
    ''' atomic energy state
         (defined by a configuration, a term, a parity and a statistical weight)
    '''
    def __init__(self, **kwargs):
        ''' float e:  energy level wrt the ground state [cm⁻¹]
            int g:    statistical weight
            str cfg:  electronic configuration
            str term: spectral term in LS, jj or j[K] coupling
            str p:    parity (even 'e', odd 'o')
            str ref:  energy level reference
            int i: ionization degree
            float sbf: photoionization cross-section at the threshold [cm⁻²]
        '''
        self.e = nan
        self.g = 0
        self.cfg = ''
        self.term = ''
        self.p = ''
        self.ref = ''
        self.e_p = None
        self.g_p = None
        self.cfg_p = None
        self.term_p = None
        self.p_p = None
        self.ref_p = None
        self.i = 1
        self.sbf = nan
        self.__dict__.update(kwargs)

        if self.cfg and not self.term:
            self.cfg, self.term = self.divide_cfg_term(fmt_term=False)

        if not self.p and self.term or self.cfg:
            self.p = self.extract_parity()

        if self.term:
            if type(self.term) == list:
                for i in range(len(self.term)):
                    self.term[i] = self.remove_label_term(self.term[i])
                self.term.sort(key=term2num)
            #else:
            #    self.term = self.remove_label_term(self.term)

        setattr(self, 'term', str(self.term))
        setattr(self, 'cfg', str(self.cfg))
        self.__format__()


    def __format__(self):
        ''' Format attribute for print in standard output and write into file
        '''
        setattr(self, 'e', float(self.e))
        setattr(self, 'g', int(float(self.g)))
        setattr(self, 'sbf', float(self.sbf))

        self.e_p = self.e.__format__('12.4f')
        self.g_p = self.g.__format__('5.0f')
        self.sbf_p = self.sbf.__format__('8.2e')

        if self.cfg:
            self.cfg_p = self.cfg.__format__('35s')
        if self.term:
            self.term_p = self.term.__format__('10s')
        if self.p:
            self.p_p = self.p.__format__('1s')
        if self.ref:
            self.ref_p = str(self.ref)

    def print(self, tot=False, ret=True):
        ''' Print all the data formatted of the instance of State
        '''
        self.__format__()
        if ret:
            print('   ', self.e_p, self.g_p, self.cfg_p, self.term_p, self.p_p, self.ref_p)
        else:
            # For printing in Line class
            print(self.e_p, self.g_p, self.cfg_p, self.term_p, self.p_p, self.ref_p, end='')

    def get(self, fmt=False):
        ''' Get the attributes of state instance
            If fmt=False, get attributes with defined types (float, int, str, ...)
            If fmt=True, get attributes formatted for printing (str only)
        '''
        if fmt:
            return self.e_p, self.g_p, self.cfg_p, self.term_p, self.p_p, self.ref_p
        else:
            return self.e, self.g, self.cfg, self.term, self.p, self.ref

    def get_S(self):
        ''' Get the multiplicity of the state
            Return a int if LS coupling
            Return a str in jK and jj coupling
            Return None if empty
        '''
        #Delete lower beggining letter in the term name
        if self.term:
            self.term = self.term.replace(' ', '')
            if self.term.find('[') != -1:
                return 'jK'
            if self.term[:1].isalpha():
                if self.term[:1].islower():
                    if self.term[1:2].isdigit():  #e.g. 'a5D'
                        return int(self.term[1:2])
                elif self.term[1:2].isdigit():    #e.g. 'A5D'
                        return int(self.term[1:2])                        
                else:                            #e.g. 'H'
                    return
            else:
                if self.term[:1].isdigit():      #e.g. '3F'
                    return int(self.term[0:1])
                elif self.term[:1] == '(':       #e.g. '(1/2,1/2)'
                    return 'jj'
        else:
            return ''

    def get_L(self):
        ''' Get the orbital momentum
             Return a capital letter (LS), or '[' (jK)
        '''
        #self.term = self.term.replace(' ', '')
        if self.term:
            for i in self.term:
                if i.isalpha():
                    if i.isupper():
                        return str(i)

    def get_pqn(self):
        ''' Get the principal quantum number from the electronic configuration
            get_pqn('3p6.4s2') returns 4
            get_pqn('3p6.4s.10p') returns 10
        '''
        if self.cfg:
            if type(self.cfg) == str:
                cfg = self.parent_cfg_wo_term()[0]
                orbs = cfg.split('.')
                norb = len(orbs)
                orb = orbs[norb-1]
                m = re.search('^[0-9]*', orb)
                if m:
                    return int(m.group()) 
            else:
                print("get_pqn: Warning self.cfg is not a string.")

    def get_sqn(self):
        ''' Get the secondary quantum number from the electronic configuration
            get_sqn('3p6.4s2') returns 0
            get_sqn('3p6.4s.5p') returns 1
        '''
        sqn_list = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm']

        if self.cfg:
            if type(self.cfg) == str:
                cfg = self.parent_cfg_wo_term()[0]
                orbs = cfg.split('.')
                norb = len(orbs)
                orb = orbs[norb-1]
                for c in orb:
                    if c.isalpha() and c.islower():
                        try:
                            return sqn_list.index(c)
                        except ValueError:
                            return None
            elif type(self.cfg) == list:
                sqn = []
                for cfg in self.cfg:
                    cfg = self.parent_cfg_wo_term()[0]
                    orbs = cfg.split('.')
                    norb = len(orbs)
                    orb = orbs[norb-1]
                    for c in orb:
                        if c.isalpha() and c.islower():
                            sqn.append(sqn_list.index(c))  
                    return sqn
            else:
                print("get_sqn: Warning self.cfg is neither a string nor a list.")

    def remove_label_term(self, term):
        '''Remove the convention label in terms of complex atomic spectra (if exist)
            i.e. seniority index
            a, b, c, ... for even-parity LS term
            z, y, x, ... for odd-parity LS term
            e.g. 'a5D' => '5D'
        '''
        if type(term) == str:
            if term[:1].isalpha():
                # Case where seniority index is a lower letter
                if term[:1].islower():
                    term = term[1:]
                # Case where seniority index is a capital letter
                elif term[1:2].isdigit():
                    term = term[1:]
                #else:
                #    term = ' '+term
            return term
        elif type(term) == list:
            term_list = []
            for t in term:
                term_list.append(self.remove_label_term(t))
            return term_list
        else:
            raise TypeError('term must be str or list of str type.')

    def divide_cfg_term(self, fmt_term=True):
        '''Divide cfg_term in cfg and term
            If fmt_term is True, remove the seniority index (lower letter before multiplicity)
        '''

        dumb = self.cfg.split('_')
        cfg = dumb[0]
        term = ''
        try:
            term = dumb[1]
        except IndexError:
            pass
        if fmt_term:
            term = self.remove_label_term(term)

        return cfg, term

    def extract_parity(self):
        '''Extract parity of term
            Return 'e' for even and 'o' for odd parities
        '''
        if self.term:
            try:
                self.term.index('*')
            except ValueError:
                return 'e'
            else:
                return 'o'
        else:
            try:
                self.cfg.index('e')
            except ValueError:
                return 'o'
            else:
                return 'e'

    def delta_e(self):
        ''' Compute the maximum difference in energy among all the energies present.
        '''
        e_list = []
        if hasattr(self, 'superlevels'):
            for i in range(self.nsl):
                if hasattr(self.superlevels[i], 'levels'):
                    for j in range(self.superlevels[i].nml):
                        if hasattr(self.superlevels[i].levels[j], 'states'):
                            for k in range(self.superlevels[i].levels[j].nfl):
                                e_list.append(self.superlevels[i].levels[j].states[k].e)
                        else:
                            e_list.append(self.superlevels[i].levels[j].e)
                else:
                    e_list.append(self.superlevels[i].e)
        else:
            e_list.append(self.e)

        return max(e_list) - min(e_list)

    def parent_cfg_wo_term(self):
        ''' Remove terms in parenthesis inside the configuration
             3d7.(4F).4p          => 3d7.4p
             3d6.(5D).4s.4p.(1P*) => 3d6.4s.4p
             3d6.4s.(6D<7/2>).4f  => 3d6.4s.4f
        '''
        pcfg = []
        if self.cfg:
            if type(self.cfg) == str:
                pcfg.append(re.sub('\.\([a-zA-Z0-9\*<>/]*[\)\.|\)$]', '', self.cfg))
            elif type(self.cfg) == list:
                for cfg in self.cfg:
                    pcfg.append(re.sub('\.\([a-zA-Z0-9\*<>/]*[\)\.|\)$]', '', cfg))
            else:
                print('Type for cfg unsupported:', type(self.cfg), self.cfg)
                quit(1)
        return pcfg


    def __eq__(self, other):
        ''' Define equality between 2 levels of class State, Level, SuperLevel or HyperLevel
            self and other are not necessarily of the same class
        '''
        de = 5.0 # cm⁻¹
        crit1, crit2, crit3, crit4 = False, False, False, False
        e_list = [other.e]
        g_list = [other.g]
        cfg_list = [other.cfg]
        term_list = [other.term]
        if hasattr(other, 'superlevels'):
            e_list.extend([sl.e for sl in other.superlevels])
            g_list.extend([sl.g for sl in other.superlevels])
            cfg_list.extend([sl.cfg for sl in other.superlevels])
            term_list.extend([sl.term for sl in other.superlevels])
            for sl in other.superlevels:
                if hasattr(sl, 'levels'):
                    e_list.extend([l.e for l in sl.levels])
                    g_list.extend([l.g for l in sl.levels])
                    term_list.extend([l.term for l in sl.levels])
                    for level in sl.levels:
                        if hasattr(level, 'states'):
                            e_list.extend([s.e for s in level.states])
                            g_list.extend([s.g for s in level.states])
                            cfg_list.extend([s.cfg for s in level.states])
                            term_list.extend([s.term for s in level.states])
        if hasattr(other, 'levels'):
            e_list.extend([l.e for l in other.levels])
            g_list.extend([l.g for l in other.levels])
            term_list.extend([l.term for l in other.levels])
            for level in other.levels:
                if hasattr(level, 'states'):
                    e_list.extend([s.e for s in level.states])
                    g_list.extend([s.g for s in level.states])
                    cfg_list.extend([s.cfg for s in level.states])
                    term_list.extend([s.term for s in level.states])
        if hasattr(other, 'states'):
            e_list.extend([s.e for s in other.states])
            g_list.extend([s.g for s in other.states])
            cfg_list.extend([s.cfg for s in other.states])
            term_list.extend([s.term for s in other.states])

        e_list = list(set(e_list))
        g_list = list(set(g_list))
        cfg_list = [str(cfg) for cfg in pl.flatten(cfg_list)]
        cfg_list = list(set(cfg_list))
        term_list = [str(term) for term in pl.flatten(term_list)]
        term_list = list(set(term_list))

        #print(e_list)
        #print(g_list)
        #print(cfg_list)
        #print(term_list)

        for e in e_list:
            if abs(self.e - e) <= de:
                crit1 = True
        if g_list.count(self.g):
            crit2 = True
        for cfg in cfg_list:
            #print(type(self.cfg), self.cfg, type(cfg), cfg)
            if type(self.cfg) == list:
                for c in self.cfg:
                    if c.strip() == cfg.strip():
                        crit3 == True
            elif type(self.cfg) == str:
                if self.cfg.strip() == cfg.strip():
                    crit3 = True
        for term in term_list:
            if type(self.term) == list:
                for t in self.term:
                    if t.strip() == term.strip():
                        crit4 = True
            elif type(self.term) == str:
                if self.term.strip() == term.strip():
                    crit4 = True
                #else:
                #   print(self.term.strip())
                #   print(term.strip())
            else:
                #print("Pb with term.", self.term, other.term, type(self.term), type(other.term), len(self.term), len(other.term))
                print("Pb with term.", self.term, other.term, type(self.term), type(other.term))
                quit(1)

        #print(crit1, crit2, crit3, crit4)
        if crit1 and crit2 and crit3 and crit4:
            return True
        else:
            return False

    def __lt__(self, other):
        return self.e < other.e

    def __gt__(self, other):
        return self.e > other.e

    def __hash__(self):
        ''' Hash must return an integer
            useful for set for example
        '''
        return int(self.e*1000)
#-------------------------------------------------------------------------------
class Level(State):
    ''' atomic energy level
        (weighted mean of states of the same configuration, term and parity)
    '''
    def __init__(self, *args, **kwargs):
        ''' args must be of the State class
            kwargs must be attributes e, g, cfg, term, p, ref, i and sbf
        '''
        State.__init__(self, **kwargs)

        if args:
            self.states = []
            for state in args:
                if isinstance(state, State):
                    self.states.append(state)
                else:
                    raise TypeError
            self.g = int(self.mean_g())
            self.e = float(self.mean_e())
            self.cfg = self.check_cfg()
            self.term = self.check_term()
            self.p = self.check_parity()
            self.ref = self.check_ref()
            self.nfl = len(self.states)
            self.i = self.states[0].i
            self.sbf = self.mean_sbf()

        self.de = self.delta_e()
        self.__format__()

    def print(self, tot=False, ret=True):
        ''' Print all the data formatted of the instance of Level
        '''
        self.__format__()

        if ret:
            print('  ', self.e_p, self.g_p, self.cfg_p, self.term_p, self.p_p, self.ref_p)
        else:
            print(self.e_p, self.g_p, self.cfg_p, self.term_p, self.p_p, self.ref_p, end='')

        if tot:
            if hasattr(self, 'states'):
                for i in range(self.nfl):
                    print('    ', end='')
                    self.states[i].print()

    def get_states(self):
        ''' Get attributes of individual states
        '''
        states = []
        if self.states:
            for i in range(self.nfl):
                states.append(self.states[i].get())
        return states

    def add_state(self, state=State()):
        ''' Add state to an existing instance of Level
        '''
        print("to do.")
        quit(1)

    def mean_g(self):
        ''' Compute mean statistical weight of the mean energy level
        '''
        if self.states:
            # Case where we have the same lower level in a multiplet
            if len(set(self.states)) == 1:
                return self.states[0].g
            else:
                g = []
                for state in self.states:
                    g.append(state.g)
                return sum(g)
        else:
            return 0

    def mean_e(self):
        ''' Compute mean energy level weighted by statistical weights
        '''
        if self.states:
            e, g = [], []
            for state in self.states:
                e.append(state.e)
                g.append(state.g)
            if set(g) == set([0]):
                return sum(e) / len(e)
            else:
                return sum(pl.array(g)*pl.array(e)) / sum(g)
        else:
            return nan

    def mean_sbf(self):
        ''' Compute mean photoionization cross-section at the threshold
        '''
        if self.states:
            sbf, g  = [], []
            for state in self.states:
                state.sbf = nan2zero(state.sbf)
                sbf.append(state.sbf)
                g.append(state.g)
            if set(g) == set([0]):
                return sum(sbf)/len(sbf)
            else:                
                return sum(pl.array(g)*pl.array(sbf)) / sum(g)
        else:
            return nan


    def check_cfg(self):
        ''' Check if the configurations of states are identical
            Return the configuration if so
        '''
        if self.states:
            cfg = []
            for state in self.states:
                if state.cfg:
                    cfg.append(state.cfg)
            if len(set(cfg)) == 1:
                return cfg[0]
            #else:
            #   print('check_cfg: ', cfg)

    def check_term(self):
        ''' Check if the terms of states are identical
            Return the term if so
        '''
        if self.states:
            term = []
            for state in self.states:
                if state.term:
                    term.append(state.term)
            if len(set(term)) == 0:
                return term
            elif len(set(term)) == 1:
                return term[0]
            #else:
            #   print('check_term:', term)

    def check_parity(self):
        ''' Check if the parity of states are identical
            Return the parity if so
        '''
        if self.states:
            parity = []
            for state in self.states:
                if state.p:
                    parity.append(state.p)
            if len(set(parity)) == 0:
                return parity
            elif len(set(parity)) == 1:
                return parity[0]
            else:
                print('check_parity', parity)

    def check_ref(self):
        ''' Check if the references are identical
             Return all the reference if None
        '''
        if self.states:
            ref = []
            for state in self.states:
                if state.ref:
                    ref.append(state.ref)
            ref = [i for i in pl.flatten(ref)]
            if len(set(ref)) == 0:
                return ref
            elif len(set(ref)) == 1:
                return ref[0]
            elif len(set(ref)) > 1:
                return ref
#-------------------------------------------------------------------------------
class SuperLevel(Level):
    ''' atomic energy superlevel
         (mean of levels of the same configuration and parity)
    '''
    def __init__(self, *args, **kwargs):
        ''' args must be of the Level class
             kwargs must be attributes e, g, cfg, term, p, ref, i and sbf
        '''
        Level.__init__(self, **kwargs)

        if args:
            self.levels = []
            for level in args:
                if isinstance(level, Level):
                    self.levels.append(level)
                else:
                    raise TypeError
            self.g = int(self.mean_g())
            self.e = float(self.mean_e())
            self.cfg = self.check_cfg()
            self.term = self.check_term()
            if type(self.term) == list:
                self.term.sort(key=term2num)
            self.p = self.check_parity()
            self.ref = self.check_ref()
            self.nml = len(self.levels)
            self.i = self.levels[0].i
            self.sbf = self.mean_sbf()

        self.de = self.delta_e()
        self.__format__()

    def print(self, tot=False, ret=True):
        ''' Print all the data formatted of the instance of SuperLevel
        '''
        self.__format__()

        if ret:
            print(' ', self.e_p, self.g_p, self.cfg_p, self.term_p, self.p_p, self.ref_p)
        else:
            print(self.e_p, self.g_p, self.cfg_p, self.term_p, self.p_p, self.ref_p, end='')

        if tot:
            if hasattr(self, 'levels'):
                for i in range(self.nml):
                    print('    ', end='')
                    self.levels[i].print(tot=True)

    def get_levels(self):
        ''' Get attributes of individual levels
        '''
        levels = []
        if self.levels:
            for i in range(self.nml):
                levels.append(self.levels[i].get())
        return levels

    def mean_g(self):
        ''' Compute statistical weight of the superlevel
        '''
        if self.levels:
            g = []
            for level in self.levels:
                g.append(level.g)
            return sum(g)
        else:
            return 0

    def mean_e(self):
        ''' Compute energy of the superlevel weighted by statistical weights
        '''
        if self.levels:
            e, g = [], []
            for level in self.levels:
                e.append(level.e)
                g.append(level.g)
            if set(g) == set([0]):
                return sum(e) / len(e)
            else:
                return sum(pl.array(g)*pl.array(e)) / sum(g)
        else:
            return nan

    def mean_sbf(self):
        ''' Compute mean photoionization cross-section at the threshold
        '''
        if self.levels:
            sbf, g  = [], []
            for level in self.levels:
                level.sbf = nan2zero(level.sbf)
                sbf.append(level.sbf)
                g.append(level.g)
            if set(g) == set([0]):
                return sum(sbf)/len(sbf)
            else:
                return sum(pl.array(g)*pl.array(sbf)) / sum(g)
        else:
            return nan

    def check_cfg(self):
        ''' Check if the configurations of levels are identical
            Return the configuration if so
        '''
        if self.levels:
            cfg = []
            for level in self.levels:
                if level.cfg:
                    cfg.append(level.cfg)
            if len(set(cfg)) == 1:
                return cfg[0]
            else:
                print('SuperLevel: check_cfg:', cfg)
                quit(1)

    def check_term(self):
        ''' Check if the terms of levels are identical
            Return the term if so
            or assembling of terms
        '''
        if self.levels:
            term = []
            for level in self.levels:
                if level.term:
                    term.append(level.term)
            if len(set(term)) == 0:
                return term
            elif len(set(term)) == 1:
                return term[0]
            else:
                return term

    def check_parity(self):
        ''' Check if the parity of levels are identical
            Return the parity if so
        '''
        if self.levels:
            parity = []
            for level in self.levels:
                if level.p:
                    parity.append(level.p)
            if len(set(parity)) == 0:
                return parity
            elif len(set(parity)) == 1:
                return parity[0]
            else:
                print('SuperLevel: check_parity', parity)

    def check_ref(self):
        ''' Check if the references are identical
            Return all the references if None
        '''
        if self.levels:
            ref = []
            for level in self.levels:
                if level.ref:
                    ref.append(level.ref)
            ref = [i for i in pl.flatten(ref)]
            if len(set(ref)) == 0:
                return ref
            elif len(set(ref)) == 1:
                return ref[0]
            elif len(set(ref)) > 1:
                return ref
#-------------------------------------------------------------------------------
class HyperLevel(SuperLevel):
    ''' atomic energy hyperlevel
        (mean of levels of the same parity with energy difference lower than delta_e)
    '''
    def __init__(self, tcop, *args, **kwargs):
        ''' de: energy range for merge of superlevels
            tcop: take care of parity if True
            args must be of the SuperLevel class
            kwargs must be attributes e, g, cfg, term, p and ref
        '''
        self.tcop = tcop
        SuperLevel.__init__(self, **kwargs)

        if args:
            self.superlevels = []
            for sl in args:
                if isinstance(sl, SuperLevel):
                    self.superlevels.append(sl)
                else:
                    raise TypeError
            self.g = int(self.mean_g())
            self.e = float(self.mean_e())
            self.cfg = self.check_cfg()
            self.term = self.check_term()
            self.nsl = len(self.superlevels)
            self.i = self.superlevels[0].i
            self.sbf = self.mean_sbf()

        if self.tcop and hasattr(self, 'superlevels'):
            self.p = self.check_parity()
        else:
            self.p = None

        self.de = self.delta_e()

        self.__format__()

    def print(self, tot=False, ret=True):
        ''' Print all the data formatted of the instance of SuperLevel
        '''
        self.__format__()
        if tot:
            if hasattr(self, 'nsl'):
                cfg = "'hyperlevel:" + str(self.nsl) + " superlevel(s)'"
            else:
                cfg = "'hyperlevel: 0 superlevel(s)'"
            cfg = cfg.__format__('35s')
            term = "'" + str(len(self.term)) + " term(s)'"
            term = term.__format__('10s')
            print('', self.e_p, self.g_p, cfg, term, self.p_p, self.ref_p)
            if hasattr(self, 'superlevels'):
                for i in range(self.nsl):
                    print('    ', end='')
                    self.superlevels[i].print(tot=True)
        else:
            if ret:
                print('', self.e_p, self.g_p, self.cfg_p, self.term_p, self.p_p, self.ref_p)
            else:
                print(self.e_p, self.g_p, self.cfg_p, self.term_p, self.p_p, self.ref_p, end='')

    def get_superlevels(self):
        ''' Get attributes of individual levels
        '''
        superlevels = []
        if self.superlevels:
            for i in range(self.nsl):
                superlevels.append(self.superlevels[i].get())
        return superlevels

    def mean_g(self):
        ''' Compute statistical weight of the hyperlevel
        '''
        if self.superlevels:
            g = []
            for sl in self.superlevels:
                g.append(sl.g)
            return sum(g)
        else:
            return 0

    def mean_e(self):
        ''' Compute energy of the hyperlevel weighted by statistical weights
        '''
        if self.superlevels:
            e, g = [], []
            for sl in self.superlevels:
                e.append(sl.e)
                g.append(sl.g)
            if set(g) == set([0]):
                return sum(e)/len(e)
            else:
                return sum(pl.array(g)*pl.array(e)) / sum(g)
        else:
            return nan

    def mean_sbf(self):
        ''' Compute mean photoionization cross-section at the threshold
        '''
        if self.superlevels:
            sbf, g  = [], []
            for sl in self.superlevels:
                sl.sbf = nan2zero(sl.sbf)
                sbf.append(sl.sbf)
                g.append(sl.g)
            if set(g) == set([0]):
                return sum(sbf)/len(sbf)
            else:
                return sum(pl.array(g)*pl.array(sbf)) / sum(g)
        else:
            return nan            

    def check_cfg(self):
        ''' Check if the configurations of superlevels are identical
            Return all the configuration if so
        '''
        if self.superlevels:
            cfg = []
            for sl in self.superlevels:
                if sl.cfg:
                    cfg.append(sl.cfg)
            if len(set(cfg)) == 1:
                return cfg[0]
            else:
                return cfg

    def check_term(self):
        ''' Check if the terms of superlevels are identical
            Return the term if so
            or assembling of terms
        '''
        if self.superlevels:
            term = []
            for sl in self.superlevels:
                if sl.term:
                    if type(sl.term) == str:
                        term.append(sl.term)
                    else:
                        term.extend(sl.term)

        return term

    def check_parity(self):
        ''' Check if the parity of superlevels are identical
            Return the parity if so
        '''
        if self.superlevels:
            parity = []
            for sl in self.superlevels:
                if sl.p:
                    parity.append(sl.p)
            if len(set(parity)) == 0:
                return parity
            elif len(set(parity)) == 1:
                return parity[0]
            else:
                print('HyperLevel: check_parity', parity)
                quit(1)
#-------------------------------------------------------------------------------
class Line:
    ''' line between two atomic energy states
    '''
    def __init__(self, **kwargs):
        ''' State lower: instance of class (State, Level, SuperLevel, HyperLevel)
            State upper: instance of class (State, Level, SuperLevel, HyperLevel)
            f:           oscillator strength
            f_acc:       accuracy (NIST notation)
            f_ref:       oscillator strength reference
            lbd:         wavelength in Å
            gr:          radiative broadening in s⁻¹
            gv:          H collision broadening (Van der Waals)
            gs:          e⁻ collision broadening (Stark)
            rtype:       type of transition (E1, E2, M1, ...)
            lev_ref:     energy level reference
            lin_ref:     energy level reference
        '''
        self.lower = None
        self.upper = None
        self.f = nan
        self.f_acc = ''
        self.lbd = nan
        self.gr = nan
        self.gv = nan
        self.gs = nan
        self.f_ref = ''
        self.rtype = ''
        self.lev_ref = ''
        self.__dict__.update(kwargs)

        setattr(self, 'f', float(self.f))
        setattr(self, 'lbd', float(self.lbd))
        setattr(self, 'gr', float(self.gr))
        setattr(self, 'gv', float(self.gv))
        setattr(self, 'gs', float(self.gs))

        if pl.isnan(self.lbd) and self.lower and self.upper:
            self.lbd = self.compute_lbd(med='air')
        #elif not pl.isnan(self.lbd):
        #    self.lbd = vac2air(self.lbd)
        #else:
        #    pass

        self.__format__()

    def __format__(self):
        self.f_p = self.f.__format__('10.3e')
        self.lbd_p = self.lbd.__format__('14.4f')
        self.gr_p = self.gr.__format__('9.2e')
        self.gv_p = self.gv.__format__('9.3f')
        self.gs_p = self.gs.__format__('9.2e')
        if self.f_acc:
            self.f_acc_p = self.f_acc.__format__('3s')
        else:
            self.f_acc_p = None
        if self.f_ref:
            self.f_ref_p = str(self.f_ref)
        else:
            self.ref_p = None
        if self.rtype:
            self.rtype_p = str(self.rtype).__format__('3s')
        else:
            self.rtype_p = None
        if self.lev_ref:
            self.lev_ref_p = self.lev_ref.__format__('15s')
        else:
            self.lev_ref_p = None

    def print(self, tot=False, lev=False, ret=False):
        ''' Print all the data formatted of the instance of Line
        '''
        self.__format__()
        if lev:
            self.lower.print(ret=ret)
            self.upper.print(ret=ret)

        print(self.f_p, self.gr_p, self.gv_p, self.gs_p, self.lbd_p, self.f_acc_p)

        if tot:
            if hasattr(self, 'lines'):
                for line in self.lines:
                    print(' ', end='')
                    line.print()

    def compute_lbd(self, med='vac'):
        ''' Compute wavelength in Angstroem
        '''
        e1 = self.lower.e
        e2 = self.upper.e
        lbd = 1.e8/abs(e1 - e2) # Angstrom in vacuum

        if med == 'air':
            lbd = vac2air(lbd)
        elif med == 'vac':
            pass
        else:
            print('compute_lbd: med not yet implemented')
            quit(1)

        return lbd

    def dn(self):
        ''' Compute the difference of the principal quantum number of the transition
        '''
        n1 = self.lower.get_pqn()
        n2 = self.upper.get_pqn()
        return abs(n2-n1)

    def f2A(self):
        ''' Convert oscillator strength into deexcitation probability
            (from Biemont - Spectroscopie atomique) 
            or exact formulation Eq. (3.84) Merle PhD report
            lbd in [A]
        '''        
        return 6.669e15*self.lower.g*self.f/(self.upper.g*self.lbd**2)
#-------------------------------------------------------------------------------
class Multiplet(Line):
    ''' transition array between two levels
    '''
    def __init__(self, *args, **kwargs):
        ''' transitions between 2 means levels
        '''
        Line.__init__(self, **kwargs)

        if args:
            self.lines = []
            lower_states_list = []
            upper_states_list = []
            for line in args:
                if isinstance(line, Line):
                    self.lines.append(line)
                    lower_states_list.append(line.lower)
                    upper_states_list.append(line.upper)
                else:
                    raise TypeError

            self.nl = len(self.lines)

            self.lower = Level(*lower_states_list)
            self.upper = Level(*upper_states_list)
            self.f = self.mean_f()
            self.lbd = self.mean_lbd()
            self.gr = self.mean_gr()
            self.gv = self.mean_gv()
            self.gs = self.mean_gs()

        self.__format__()

    def mean_f(self):
        ''' Compute the multiplet oscillator strength
        '''
        gf = 0.
        if hasattr(self, 'lines'):
            for line in self.lines:
                #print(line.lower.g)
                f_temp = nan2zero(line.f)
                gf += line.lower.g * f_temp
            #print(self.lower.g)
            return gf/self.lower.g

    def mean_gr(self):
        ''' Compute the multiplet radiative broadening
        '''
        gr = 0.
        if hasattr(self, 'lines'):
            for i in range(self.nl):
                gr_temp = nan2zero(self.lines[i].gr)
                gr += gr_temp
            return gr

    def mean_gv(self):
        ''' Compute the mean H broadening
        '''
        gv = 0.
        if hasattr(self, 'lines'):
            for i in range(self.nl):
                gv_temp = nan2zero(self.lines[i].gv)
                gv += gv_temp
            return gv/self.nl

    def mean_gs(self):
        ''' Compute the mean e⁻ broadening
        '''
        gs = 0.
        if hasattr(self, 'lines'):
            for i in range(self.nl):
                gs_temp = nan2zero(self.lines[i].gs)
                gs += gs_temp
            return gs/self.nl

    def mean_lbd(self):
        ''' Compute the mean wavelength
        '''
        lbd = 0.
        if hasattr(self, 'lines'):
            for line in self.lines:
                lbd += line.lbd
            lbd = lbd/float(self.nl)
            return lbd
#-------------------------------------------------------------------------------
class Collision:
    ''' Collision transition between two energy states
    '''
    def __init__(self, **kwargs):
        self.type = ''
        self.lower = None
        self.upper = None
        self.t_list = []
        self.ups_list = []
        self.lbd = nan
        self.ref = ''
        self.kw = ''
        self.__dict__.update(kwargs)

        self.nt = len(self.t_list)
        self.nups = len(self.ups_list)

        setattr(self, 'lbd', float(self.lbd))

        if pl.isnan(self.lbd) and self.lower and self.upper:
            self.lbd = self.compute_lbd(med='air')
        #elif not pl.isnan(self.lbd):
        #    self.lbd = vac2air(self.lbd)
        #else:
        #    pass

        self.__format__()

    def __format__(self):
        self.lbd_p = self.lbd.__format__('10.3f')
        self.t_list_p = []
        self.ups_list_p = []
        for i in range(self.nt):
            self.t_list_p.append(self.t_list[i].__format__('9.2e'))
        for i in range(self.nups):
            self.ups_list_p.append(self.ups_list[i].__format__('9.2e'))


    def print(self, tot=False, lev=False, ret=False):
        ''' Print all the data formatted of the instance of Collision
        '''
        self.__format__()
        
        if lev:
            self.lower.print(ret=ret)
            self.upper.print(ret=ret)

        if tot:
            print(format(self.kw, '<12'), self.ref, self.lbd_p, self.type)
            if hasattr(self, 'l') and hasattr(self, 'u'):
                ups = ''
                for i in range(self.nups):
                    ups += ' '+str(self.ups_list_p[i])
                print(str(self.l).rjust(4), str(self.u).rjust(4), ups)
            else:
                print('self.l and self.u required to print with tot = True')
                quit(1)
        else:
            print(self.ups_list_p)

    def compute_lbd(self, med='vac'):
        ''' Compute wavelength in Angstroem
        '''
        e1 = self.lower.e
        e2 = self.upper.e
        if e1 == e2:
            print("Collision: compute_lbd: e1 == e2:", e1, ", do e2=e1+0.001 cm⁻1")
            e2 = e1 + 0.001
            #quit(1)
        lbd = 1.e8/abs(e1 - e2) # Angstrom in vacuum

        if med == 'air':
            lbd = vac2air(lbd)
        elif med == 'vac':
            pass
        else:
            print('compute_lbd: med not yet implemented')
            quit(1)

        return lbd
#-------------------------------------------------------------------------------
def term2num(term):
    ''' Convert spectral term into corresponding number
    '''
    L = ['S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'Q', 'R', 'T', 'U', '[', '(']
    if term:
        if len(term) == 1 or term[1] == '*' or term[0] == '(':
            t = term[0]
        else:
            t = term[1]
        if t in L: 
            return L.index(t)+1 
#-------------------------------------------------------------------------------
def L2num(L):
    ''' Convert spectral term L into corresponding number
    '''
    L_list = ['S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'Q', 'R', 'T', 'U', '[', '(']
    if L in L_list: 
        return L_list.index(L)+1 
    else:
        return None
#-------------------------------------------------------------------------------
def term2mult(term):
    ''' Convert spectral term into corresponding multiplicity
    '''
    mult_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9',\
                 'S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'Q', 'R', 'T', 'U',\
                 '1[','2[', '3[', '4[', '5[', '6[', '7[', '8[', '9[', '(', '']
    if term[0:2] in mult_list:
        return mult_list.index(term[0:2]) + 1
    elif term[0] in mult_list:
        return mult_list.index(term[0]) + 1 
#-------------------------------------------------------------------------------
def none2nan(chain):
    ''' Convert 'None' string in pl.nan
        return the same string otherwise
    '''
    if str(chain).strip() == 'None':
        return nan
    else:
        return chain
#-------------------------------------------------------------------------------
def nan2zero(val):
    ''' Convert pl.nan to 0.0 float
        return the same value otherwise
    '''
    if val == val:
        return val
    else:
        return 0.
#-------------------------------------------------------------------------------
def vac2air(x):
    ''' Conversion of wavelengths from vacuum to the air in Angström in the standard condition
        of pressure and temperature (IAU Morton (1991, ApJS, 77, 119))
        http://www.sdss.org/DR6/products/spectra/vacwavelength.html
    '''
    return x / (1.0 + 2.735182E-4 + 131.4182E0 / x**2 + 2.76249E8 / x**4)
#-------------------------------------------------------------------------------
def air2vac(x):
    ''' Conversion of wavelengths from vacuum to the air in Angström in the standard condition
        of pressure and temperature (IAU Morton (1991, ApJS, 77, 119))
        http://www.sdss.org/DR6/products/spectra/vacwavelength.html
    '''
    return x * (1.0 + 2.735182E-4 + 131.4182E0 / x**2 + 2.76249E8 / x**4)   
#-------------------------------------------------------------------------------
def neff(ion, level):
    ''' Compute effective principal quantum number
    '''
    i = level.i
    e = 100 * Cst.H * Cst.C / Cst.Q * (ion.e - level.e) # Ionization energy of the level in eV
    return i * pl.sqrt(13.6057/e) 
#-------------------------------------------------------------------------------
def orbital_radius(ion, level):
    ''' Normally only valid for hydrogenoid ions
    '''
    ns = neff(ion, level)
    l = level.get_sqn()
    if type(l) == int:
        return (3*ns**2-l*(l+1)) / (2*level.i) # in unit of Cst.A0 (Bohr radius)
    elif type(l) == list:
        or_list = []
        for ll in l:
            or_list.append((3*ns**2-ll*(ll+1)) / (2*level.i))
        return pl.mean(or_list)

#-------------------------------------------------------------------------------

# First energy ionization of each species
ions = {  'oi': State(e=109837.026, g=4, cfg='2s2.2p3', termp='4S*', p='o', ref='NIST', i=2),
          'oii': Level(State(e=283270.95, g=1, cfg='2s2.2p2', term='3P', p='e', ref='NIST', i=3),\
                       State(e=283384.1, g=3, cfg='2s2.2p2', term='3P', p='e', ref='NIST', i=3),\
                       State(e=283577.1, g=5, cfg='2s2.2p2', term='3P', p='e', ref='NIST', i=3)), 
          'nai':  State(e=41449.451, g=1, cfg='2p6', term='1S', p='e', ref='NIST', i=2),
          'naii': Level(State(e=381390.2, g=4, cfg='2s2.2p5', term='2P*', p='o', ref='NIST', i=3), \
                        State(e=382756.5, g=3, cfg='2s2.2p5', term='2P*', p='o', ref='NIST', i=3)),
          'mgi':  State(e=61671.05, g=2, cfg='3s', term='2S', p='e', ref='NIST', i=2),
          'mgii': State(e=121267.61, g=1, cfg='2p6', term='1S', p='e', ref='NIST', i=3),
          'ali': State(e=48278.48, g=2, cfg='3s2', term='1S*', p='o', ref='NIST', i=2),
	      'ki':   State(e=35009.814, g=1, cfg='3p6', term='1S', p='e', ref='NIST', i=2),
          'kii':  State(e=255072.8, g=4, cfg='3p5', term='2P*', p='o', ref='NIST', i=3),
          'cai':  State(e=49305.9240, g=2, cfg='4s', term='2S', p='e', ref='NIST', i=2),
          'caii': State(e=95751.87, g=1, cfg='3p6', term='1S', p='e', ref='NIST', i=3),
          'fei':  State(e=63737.704,g=10, cfg='3d6.(5D).4s', term='6D', p='e', ref='NIST', i=2),
          'feii': State(e=130655.4, g=9, cfg='3d6', term='5D', p='e', ref='NIST', i=3),
          'bai': State(e=42034.91, g=2, cfg='6s', term='2S', p='e', ref='NIST', i=2),
          'baii': State(e=80686.30, g=1, cfg='5p6', term='1S', p='e', ref='NIST',i=3)}

# Enhancement factor for H broadening
fh = {'nai': 2.0, 'naii': 2.5, \
        'mgi': 2.5, 'mgii': 2.5, \
	'ali': 2.5,\
        'sii': 1.3, \
        'ki':  2.5, 'kii':  2.5, \
        'cai': 1.8, 'caii': 1.4, \
        'fei': 1.4, 'feii': 2.5}

#Theoretical ionization energy from TIPTOP base
eryd = {'fei': 0.56999, 'feii': 1.1782}



if __name__ == '__main__':

    #Test the dictionary of ionization stage
    print('\nTest dictionary ions.\n')
    for key, val in ions.items():
        print(key.__format__('5s'),': ', end='')
        val.print()

    mgi_ion = ions.get('mgi')

    #Test state class
    print('\nTest State class.\n')

    state1 = State(e=0.0000, g=9, cfg='3d6.4s2', term='5D', p='e', ref='NIST')
    state2 = State(e=415.9330, g=7, cfg='3d6.4s2', term='a5D', p='e', ref='NIST')
    state3 = State(e=704.0070, g=5, cfg='3d6.4s2', term='5D', ref='NIST')
    state4 = State(e=888.1320, g=3, cfg='3d6.4s2_5D', p='e', ref='NIST')
    state5 = State(e=978.0740, g=1, cfg='3d6.4s2_a5D', p='e', ref='NIST')

    state6 = State(e=19350.8910, g=11, cfg='3d6.(5D).4s.4p.(3P*)', term='7D*')
    state7 = State(e=19562.4390, g=9, cfg='3d6.(5D).4s.4p.(3P*)', term='7D*')
    state8 = State(e=19757.0320, g=7, cfg='3d6.(5D).4s.4p.(3P*)', term='7D*')
    state9 = State(e=19912.4950, g=5, cfg='3d6.(5D).4s.4p.(3P*)', term='7D*')
    state10 = State(e=20019.6350, g=3, cfg='3d6.(5D).4s.4p.(3P*)', term='7D*')

    for state in [state1, state2, state3, state4, state5, state6, state7, state8, state9, state10]:
        state.print()

    state1 = State(e=21850.405, g=1, cfg='3s.3p', term='3P*', p='o', ref='NIST')
    state2 = State(e=58661.9780, g=9, cfg='58661e')
    state3 = State(e=21853.405, g=1, cfg='3s.3p', term='3P*')

    state = state1

    print('Test method print():              ',), state.print()
    print('Test method get():                   ', state.get())
    print('Test method get(True):               ', state.get(True))
    print('Test method get_S():                 ', state.get_S())
    print('Test method get_L():                 ', state.get_L())
    print('Test method remove_label_term("a5D"):', state.remove_label_term('a5D'))
    print('Test method divide_cfg_term():       ', state.divide_cfg_term())
    print('Test method divide_cfg_term(True):   ', state.divide_cfg_term(fmt_term=True))
    print('Test method extract_parity():        ', state.extract_parity())
    print('Test method __lt__() and __gt__():   ', state < state1, state > state1)
    print('Test method __eq__():                ', state == state2, state == state3)

    #Test Level class
    print('\nTest Level class.\n')

    level1 = Level(state1, state2, state3)
    # Print the mean level
    level1.print()
    # Print the fine structure
    level1.print(tot=True)

    print('Test method get_states():  ', level1.get_states())
    print('Test method mean_g():      ', level1.mean_g())
    print('Test method mean_e():      ', level1.mean_e())
    print('Test method check_cfg():   ', level1.check_cfg())
    print('Test method check_term():  ', level1.check_term())
    print('Test method check_parity():', level1.check_parity())
    print('Test method check_ref():   ', level1.check_ref())

#   ground = Level(ground)
#   level_list = [ground, level1]
#   level_test = level_list[0]
#   level_test.print()
#   print(level_test.get_states())
#   print(level_test.get_S())

    # Test SuperLevel class
    print('\nTest SuperLevel class.\n')

    level2 = Level(state2)

    sl1 = SuperLevel(level1, level2)
    sl2 = SuperLevel(Level(e=3000, g=12, cfg='toto', term='3H*'), Level(e=3001, g=8, cfg='toto', term='5P*'))

    # Print the superlevel
    sl2.print(tot=True)

    print('Test method get_levels():  ', sl1.get_levels())
    print('Test method mean_g():      ', sl1.mean_g())
    print('Test method mean_e():      ', sl1.mean_e())
    print('Test method check_cfg():   ', sl1.check_cfg())
    print('Test method check_term():  ', sl1.check_term())
    print('Test method check_parity():', sl1.check_parity())
    print('Test method check_ref():   ', sl1.check_ref())

    # Test HyperLevel class
    print('\nTest HyperLevel class.\n')

    state4 = State(e=62079.3220, g=13, cfg='62079e', term=None, p=None, ref=None)
    state5 = State(e=62192.7350, g=11, cfg='62192e', term=None, p='e', ref=None)

    level4 = Level(state4)
    level5 = Level(state5)

    sl2 = SuperLevel(level4)
    sl3 = SuperLevel(level5)
    sl_list = [sl1, sl2, sl3]
    de = 200.

    hyperlevel = HyperLevel(False, *sl_list)

    hyperlevel.print(tot=True)
    print(hyperlevel.delta_e())

    HyperLevel(True, *[sl1]).print(tot=True)

    # Test of list of instance
    print('\nTest of list of State instance.\n')
    state_list = [state1, state2]
    print('state_list[0].e', state_list[0].e)

    # Test Line class
    print('\nTest Line class.\n')

    ground = State(e=0.000, g=1, cfg='2p6.3s2', term='1S', p='e', ref='NIST')
    state1 = State(e=21850.405, g=1, cfg='3s.3p', term='3P*', p='o', ref='NIST')
    state2 = State(e=21870.464, g=3, cfg='3s.3p', term='3P*', p='o', ref='NIST')
    state3 = State(e=21911.178, g=5, cfg='3s.3p', term='3P*', p='o', ref='NIST')
    ion = State(e=61671.05, g=2, cfg='2p6.3s', term='2S', p='e', ref='NIST')

    l1 = Line(lower=ground, upper=ion, f=1.0, f_acc='A', lbd=1e10/(float(ion.e)-float(ground.e))/100., gr=1e8, gv=1.50, ref='test')
    l2 = Line(lower=ground, upper=state1, f=0.1, lbd=5000., gr=1e8, gv=1.50, ref='test')

    l1 = Line(lower=ground, upper=state1, f=10e-11, f_acc='B+', lbd=4563.881, ref='test')
    l2 = Line(lower=ground, upper=state2, f=2.382e-6, f_acc='D', lbd=4572.3767, gr=2.54e+02, ref='NIST')
    l3 = Line(lower=ground, upper=state3, f=1.820e-12, f_acc='D', lbd=4576.574, gr=5.8e-04, ref='NIST')

    l1.print(tot=True, lev=True)

    # Test Multiplet class
    print('\nTest Multiplet class.\n')

    m1 = Multiplet(lower=level1, upper=level2, f=1.00, gr=1e8)
    m1.print(tot=True, lev=True, ret=True)

    print()

    m2 = Multiplet(l1, l2, l3)
    m2.print(tot=True)

    # Test Collision class
    print('\n Test Collision class.\n')

    c1 = Collision(lower=state1, upper=state2, type='ALLOWED', \
                    t_list=[2000, 4000, 6000, 8000], ups_list=[0.1, 0.2, 0.3, 0.5], lbd=5000, ref='Fisher', kw='TEP')

    c1.l = 2
    c1.u = 10 
    c1.print(tot=True)
