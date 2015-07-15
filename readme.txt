FoRMaTo2.0 beta (work in progress)

To get information about mael.py, mart.py etc. invoked them with the -h option.


mad.py    Merge Atomic Data module
------ 

Definition of constants

  Cst        class of the phyisycal constants and atomic masses

Definition of classes for levels and lines

	State      (fine level) 
	Level      (mean level - weighted mean of fine level with same configuration and same spectral term)
	SuperLevel (mean of Levels with same configuration but different spectral terms)
	HyperLevel (mean of SuperLevels with difference energy lower than -d value and 
              with same parity (except if -p option is set))
  Line       (line as defined in atomic line database (between 2 states))
  Multiplet  (combination of lines between Levels/SuperLevel/HyperLevel)

Definition of data dictionnaries

  ions     Energy level of ionization stage
  fh       Default enhancement factor for broadening of lines by H collisions in the Unsold formulation

Definition of useful fonctions

  none2nan(chain)   convert 'None' into nan 
  nan2zero(val)     convert nan to 0  
  vac2air(x)        convert x (in A) from vaccum to air


mael.py  Tool to merge atomic energy levels
-------

usage: mael.py [-h] [-fl FL_THRES] [-sl SL_THRES] [-hl HL_THRES] [-d DELTA_E]
               [-l E_LIM] [-i] [-p] [-v] [-t] [-o OFFSET] [-e EPSILON]
               species degree filename

Tool to merge atomic energy levels.

positional arguments:
  species               Symbol of the atomic species (e.g. "Fe" or "FE" or
                        "fe")
  degree                Ionization degree (e.g. "I" or "i" for neutral
                        species)
  filename              Input atomic energy level filename [e[cm⁻¹] g
                        config_term | e[cm⁻¹] g config term | e[cm⁻¹] g
                        config term ref]

optional arguments:
  -h, --help            show this help message and exit
  -fl FL_THRES, --fl_thres FL_THRES
                        energy threshold [eV | cm⁻¹] above which fine
                        levels are merged in mean levels (default: 0)
  -sl SL_THRES, --sl_thres SL_THRES
                        energy threshold [eV | cm⁻¹] above which mean
                        levels are merged into super levels (default:
                        ionization)
  -hl HL_THRES, --hl_thres HL_THRES
                        energy threshold [ev | cm⁻¹] above which
                        superlevels are merged into hyperlevels (default:
                        ionization)
  -d DELTA_E, --delta_e DELTA_E
                        energy gap [eV | cm⁻¹] below which levels are
                        merged in hyperlevels
  -l E_LIM, --e_lim E_LIM
                        energy limit [eV | cm⁻¹] below which energy levels
                        are considered
  -i, --ionization      Only if -l option is used - Add the ionization level
                        from mad.ions
  -p, --do_not_take_care_of_parity
                        Do not take care of parity in merging superlevels
  -v, --verbose         Write detailed structure in output files
  -t, --temporary       Write temporary files (mael_fl, mael_ml, mael_sl,
                        mael_hl, mael_terms, mael_pcfg)
  -o OFFSET, --offset OFFSET
                        Decrease input energy levels with this offset value.
  -e EPSILON, --epsilon EPSILON
                        Add an espilon to energy if 2 levels have exactely the
                        same values (default: 0.01)


mart_bb.py  Tool to merge atomic bound-bound radiative transitions
----------

usage: mart_bb.py [-h] [-i INPUT] [-s SOURCE] [-f] [-abo]
                  [-fh ENHANCEMENT_FACTOR] [-v] [-t]
                  species degree filename

Tool to merge atomic bound-bound radiative transitions.

positional arguments:
  species               Symbol of the atomic species (e.g. "Fe" or "FE" or
                        "fe")
  degree                Ionization degree (e.g. "I" or "i" for neutral
                        species)
  filename              Bound-bound radiative transitions file

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Atomic energy levels file (in output format of mael.py
                        - default ael_xxx.bin)
  -s SOURCE, --source SOURCE
                        Source of atomic line database [nist | kald | vald |
                        multi] (default nist)
  -f, --f_nan           Select also lines with no f values
  -abo, --abo_theory    Apply ABO theory where possible.
  -fh ENHANCEMENT_FACTOR, --enhancement_factor ENHANCEMENT_FACTOR
                        enhancement factor for Van der Waals broadening
                        according to Unsold
  -v, --verbose         Write detailed structure in output files
  -t, --temporary       Write temporary files (mart_bb_init, mart_bb_ll)

mart_bb.py uses barklem.f90 routine which need to be compiled and encapsulated using f2py:


f2py -h barklem.pyf -m barklem barklem.f90
f2py -c --fcompiler=gfortran -m barklem barklem.pyf barklem.f90


mart_bf.py  Tool to merge atomic bound-free transitions
----------
usage: mart_bf.py [-h] [-i INPUT] [-qmp QUANTUM_MECHANICAL_PHOTOIONIZATION]
                  [-de DELTA_E_EV] [-c CUTOFF] [-nus N_UNDER_SAMPLING]
                  [-n_he N_HIGH_ENERGY]
                  [-nus_he N_UNDER_SAMPLING_AT_HIGH_ENERGY] [-f FREQUENCY_CUT]
                  [-l ENERGY_LIMIT] [-scp]
                  species degree

Tool to merge atomic bound-free radiative transitions.

positional arguments:
  species               Symbol of the atomic species (e.g. "Fe" or "FE" or
                        "fe")
  degree                Ionization degree (e.g. "I" or "i" for neutral
                        species)

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Atomic energy levels file (in output format of mael.py
                        - default ael_xxx.bin)
  -qmp QUANTUM_MECHANICAL_PHOTOIONIZATION, --quantum_mechanical_photoionization QUANTUM_MECHANICAL_PHOTOIONIZATION
                        Include TIP/TOPBASE photionisation cross-sections if
                        available [topbase|norad]
  -de DELTA_E_EV, --delta_e_eV DELTA_E_EV
                        Useful only if -qmp is used, absolute energy drift
                        tolerance between exp and theoretical values (in eV)
                        for identification
  -c CUTOFF, --cutoff CUTOFF
                        Useful only if -qmp is used, cutoff [0, 1] for the
                        identification of configuration and term, default=0.44
  -nus N_UNDER_SAMPLING, --n_under_sampling N_UNDER_SAMPLING
                        Useful only if -qmp is used, step for undersampling of
                        qm data, default=150
  -n_he N_HIGH_ENERGY, --n_high_energy N_HIGH_ENERGY
                        Useful only if -qmp is used, number of points at high
                        energy keeped for no smoothing, default=50 (used only
                        if not given in NORAD/TOPBASE data).
  -nus_he N_UNDER_SAMPLING_AT_HIGH_ENERGY, --n_under_sampling_at_high_energy N_UNDER_SAMPLING_AT_HIGH_ENERGY
                        Useful only if -qmp is used, step for undersampling of
                        qm data at high energy
  -f FREQUENCY_CUT, --frequency_cut FREQUENCY_CUT
                        Useful only if -qmp is used, frequency cut of the
                        smoothing function base on a third order Butterworth
                        low pass filter
  -l ENERGY_LIMIT, --energy_limit ENERGY_LIMIT
                        Useful only if -qmp is used, wavelenght below which
                        photoionization cross-sections are not considered (in
                        Angstrom), default=911
  -scp, --semi_classical_photoionization
                        Use Kramers formula with Gaunt factor


mart_bf.py use janicki.f90 and px.f90 which need to be encapsulated with f2py:

f2py -h janicki.pyf -m janicki janicki.f90
f2py -c -m janicki janicki.pyf janicki.f90

f2py -m px -c px.f90 --fcompiler=gfortran


mact_e.py  Tool to merge atomic collision transition with electrons
----------

usage: mact_e.py [-h] [-f1 FILENAME1] [-f2 FILENAME2] [-lte]
                 [-sce SEMI_CLASSICAL_E] [-flow FORBIDDEN_LINE]
                 [-qme QUANTUM_MECHANICAL_E]
                 [-t [TEMPERATURE [TEMPERATURE ...]]] [-de DELTA_E] [-p]
                 species degree

Tool to merge atomic collision transition with electrons

positional arguments:
  species               SYMbol of the atomic species (e.g. "Fe" or "FE" or
                        "fe")
  degree                Ionization degree (e.g. "I" or "i" for neutral
                        species)

optional arguments:
  -h, --help            show this help message and exit
  -f1 FILENAME1, --filename1 FILENAME1
                        Filename of energy levels produced by mael.py
  -f2 FILENAME2, --filename2 FILENAME2
                        Filename of radiative b-b transitions produced by
                        mart_bb.py
  -lte, --lte           Enforce collisions to ensure LTE
  -sce SEMI_CLASSICAL_E, --semi_classical_e SEMI_CLASSICAL_E
                        Choose semi-classical formula for electron collisions
                        [Seaton | VanRegemorter | Fisher]
  -flow FORBIDDEN_LINE, --forbidden_line FORBIDDEN_LINE
                        Useful only if -sce is used, effective collision
                        strength of 1.0 for semi-forbidden lines with f lower
                        than -flow
  -qme QUANTUM_MECHANICAL_E, --quantum_mechanical_e QUANTUM_MECHANICAL_E
                        Choose quantum mechanical data if available for
                        electron collisions
  -t [TEMPERATURE [TEMPERATURE ...]], --temperature [TEMPERATURE [TEMPERATURE ...]]
                        Temperature list
  -de DELTA_E, --delta_e DELTA_E
                        allowed energy gap [cm⁻¹] for identification with
                        quantum mechanical data (useful when dealing with fine
                        structure
  -p, --plot            Plot the effective collisions strengths


mact_h.py  Tool to merge atomic collision transition with H atoms
---------

usage: mact_h.py [-h] [-f1 FILENAME1] [-f2 FILENAME2] [-lte]
                 [-sch SEMI_CLASSICAL_H] [-sh SCALING_FACTOR] [-fl]
                 [-qmh QUANTUM_MECHANICAL_H]
                 [-t [TEMPERATURE [TEMPERATURE ...]]] [-de DELTA_E] [-p]
                 species degree

Tool to merge atomic collision transition with neutral hydrogen atoms

positional arguments:
  species               SYMbol of the atomic species (e.g. "Fe" or "FE" or
                        "fe")
  degree                Ionization degree (e.g. "I" or "i" for neutral
                        species)

optional arguments:
  -h, --help            show this help message and exit
  -f1 FILENAME1, --filename1 FILENAME1
                        Filename of energy levels produced by oomael.py
  -f2 FILENAME2, --filename2 FILENAME2
                        Filename of radiative b-b transitions produced by
                        oomart_bb.py
  -lte, --lte           Enforce collisions to ensure LTE
  -sch SEMI_CLASSICAL_H, --semi_classical_h SEMI_CLASSICAL_H
                        Choose semi-classical formula for electron collisions
                        [Drawin]
  -sh SCALING_FACTOR, --scaling_factor SCALING_FACTOR
                        Scaling factor for semi-classic Drawin formula
  -fl, --forbidden_lines
                        Include forbidden lines using Drawin formula with f=1
  -qmh QUANTUM_MECHANICAL_H, --quantum_mechanical_h QUANTUM_MECHANICAL_H
                        Choose quantum mechanical data if available for
                        hydrogen collisions
  -t [TEMPERATURE [TEMPERATURE ...]], --temperature [TEMPERATURE [TEMPERATURE ...]]
                        Temperature list
  -de DELTA_E, --delta_e DELTA_E
                        allowed energy gap [cm⁻¹] for identification with
                        quantum mechanical data (useful when dealing with fine
                        structure
  -p, --plot            Plot the effective collisions strengths



formato.py  Tool to create model atom from previous output binary files
----------

usage: formato.py [-h] species [degree [degree ...]]

FoRMATo Tool to create model atom for NLTE radiative transfer in late-type
stars.

positional arguments:
  species     Symbol of the atomic species (e.g. "Fe" or "FE" or "fe")
  degree      Ionization degrees (e.g. "I", "II" or "I II" for neutral,
              ionized or combined species)

optional arguments:
  -h, --help  show this help message and exit



pgrot.py  Tool to plot Grotrian diagram
--------

usage: pgrot.py [-h] [-f1 FILENAME1] [-f2 FILENAME2] [-f3 FILENAME3] [-rt]
                [-ct] [-pcfg] [-t TITLE] [-s [SIZE [SIZE ...]]] [-e EXT]
                [-l LEVEL_NAME] [-g] [-lbd WAVELENGTH] [-nop]
                species degree

Tool to plot Grotrian diagram.

positional arguments:
  species               Symbol of the atomic species (e.g. "Fe" or "FE" or
                        "fe")
  degree                Ionization degree (e.g. "I" or "i" for neutral
                        species)

optional arguments:
  -h, --help            show this help message and exit
  -f1 FILENAME1, --filename1 FILENAME1
                        Filename of energy levels produced by mael.py
  -f2 FILENAME2, --filename2 FILENAME2
                        Filename of radiative b-b transitions produced by
                        mart_bb.py
  -f3 FILENAME3, --filename3 FILENAME3
                        Filename of collision transitions produced by
                        macte_e.py
  -rt, --rad_transitions
                        Plot radiative transitions
  -ct, --col_transitions
                        Plot collisions transitions
  -pcfg, --parent_configuration
                        Plot Grotrian diagram as a function of parent
                        configurations (default: terms)
  -t TITLE, --title TITLE
                        Title of the figure
  -s [SIZE [SIZE ...]], --size [SIZE [SIZE ...]]
                        Size of the graphic [xsize, ysize]
  -e EXT, --ext EXT     Type of extension [png | eps | pdf] (default: png)
  -l LEVEL_NAME, --level_name LEVEL_NAME
                        Energy in eV below which level names are labelled
  -g, --grid            Plot the grid
  -lbd WAVELENGTH, --wavelength WAVELENGTH
                        Energy in eV below which wavelengths are labelled
  -nop, --no_plot       Do not display interactive plot


Examples of use
---------------
fl = fine levels
ml = mean levels
sl = superlevels
hl = hyperlevels


WARNING: 
There is an increasing energy order to follow when using -fl, -sl and -hl options.
Basically -fl < -sl < -hl.
If this is not followed mael.py can produce strange output. Always check the standard outputs!!!

The -hl option need the -d option to work properly and cannot be use alone (without -fl or -sl).



- Model not taking into account fine structure (default model)

mael.py Fe I ad/ael/nist_fei.dat        => 315 ml (corresponding to 846 fl) + ionization stage
mart_bb ad/rt/nist_fei.dat              => 777 multiplets (corresponding to 2402 transitions)



- Model taking into account all fine structure

mael.py Fe I ad/ael/nist_fei.dat -fl 15 => 846 fl + ionization stage
mart_bb ad/rt/nist_fei.dat              => 2402 transtions



- Model taking into account fine levels below 1 eV 
  (basically only fine levels of the ground stage and the first excited level of Fe I)

mael.py Fe I ad/ael/nist_fei.dat -fl 1  => 324 = 10 fl + 313 ml + ionization level
mart_bb.py ad/rt/nist_fei.dat           => 981 multiplets (corresponding to 2402 transitions)



- Model taking into account fine levels below 1 eV and superlevels above 5 eV

mael.py Fe I ad/ael/nist_fei.dat -fl 1 -sl 5  => 118 = 10 fl + 40 ml + 67 sl + ionization stage
mart_bb.py ad/rt/nist_fei.dat                 => 508 multiplets (corresponding to 2402 transtions)



- Model taking into account fine levels below 1 eV, superlevels above 5 eV and hyperlevels above 6.5 eV (with sl energy diffence lower than 0.05 eV )

mael.py Fe I ad/ael/nist_fei.dat -fl 1 -sl 5 -hl 6.5 -d 0.05  =>  76 = 10 fl + 40 ml + 13 sl + 12 hl + ionization stage
mart_bb.py ad/rt/nist_fei.dat                                 =>  507 multiplets (corresponding to 2402 transitions)



To compile and encapsulate the fortran routines
-----------------------------------------------

f2py -h barklem.pyf -m barklem barklem.f90
f2py -c --fcompiler=gfortran -m barklem barklem.pyf barklem.f90

f2py -h janicki.pyf -m janicki janicki.f90
f2py -c -m janicki janicki.pyf janicki.f90

f2py -m px -c px.f90 --fcompiler=gfortran
