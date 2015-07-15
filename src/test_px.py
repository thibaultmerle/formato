#!/usr/bin/python

# Plot the photoionization for all the data

from __future__ import print_function

import pylab as pl
import scipy.signal as ss
from px import format_px, extract_px
from mart_bf import id_uniq_qm_px, px_smooth, plot_px
from mad import *

flag = 'topbase'

#px_org_fn = '../ad/rt/bf/norad_fei_px.dat'
#px_bin_fn = '../ad/rt/bf/norad_fei_px.bin'
#idx_fn = '../ad/rt/bf/norad_fei_idx.dat'
#opath='../run/fei/plot/bf/'

#px_org_fn = '../ad/rt/bf/norad_feii_px.dat'
#px_bin_fn = '../ad/rt/bf/norad_feii_px.bin'
#idx_fn = '../ad/rt/bf/norad_feii_idx.dat'
#opath='../run/feii/plot/bf/'

#px_org_fn = '../ad/rt/bf/topbase_mgi_px.dat'
#px_bin_fn = '../ad/rt/bf/topbase_mgi_px.bin'
#idx_fn = '../ad/rt/bf/topbase_mgi_idx.dat'
#opath = './'

ael_fn = '../ad/ael/topbase_nai.dat'
px_org_fn = '../ad/rt/bf/topbase_nai_px.dat'
px_bin_fn = '../ad/rt/bf/topbase_nai_px.bin'
idx_fn = '../ad/rt/bf/topbase_nai_idx.dat'
opath = './'

CUTOFF = 1
DE = 0.5
N_HE0 = 20
N_US = 500
N_US_HE = 10
WC = 0.01
MPI = False

# Test for Na I

state1 = State(e=30266.99, g=2, cfg='2p6.4p', term='2P*')   
state2 = State(e=30272.58, g=4, cfg='2p6.4p', ter='2P*')

level = Level(state1, state2)

idx_table = pl.loadtxt(idx_fn, dtype='str')

# Read NORAD/TOPBASE ael data for mean levels (i.e. without fine structure)
try:
    dumb = pl.loadtxt(ael_fn, dtype='str')
except ValueError:
    print("Problem reading input file:", ael_fn)
    quit(1)       
    
if dumb.shape[1] == 4:
    qm_levels = [Level(e=e, g=g, cfg=cfg, term=term, ref=QMP) for e, g, cfg, term in dumb]
elif dumb.shape[1] == 6:
    qm_levels = [Level(e=e, g=g, cfg=cfg, term=term, p=p, ref=ref) for e, g, cfg, term, p, ref in dumb]
else:
    print("Input format problem.")
    quit(1)  


i, e, x = id_uniq_qm_px(level, qm_levels, idx_table, CUTOFF, DE, px_bin_fn)

e_us, x_us, area_i, area_is = px_smooth(i, e, x, idx_table, N_HE0, N_US, N_US_HE, WC)

# Remove negative cross-sections
ef = pl.compress( x_us > 0., e_us)
xf = pl.compress( x_us > 0., x_us)

plot_px(ef, xf, [area_i,], [area_is,], level, qm_levels, idx_table, ((i,e,x),), N_HE0, MPI)




# def smooth(x,wc):
#    ''' Smooth data x by a third order Butterworth low-band filter characterized by a cut frequency wc [rad/s] 
#    ''' 
#    b,a = ss.butter(3, wc)
#    return ss.filtfilt(b, a, x)
# 
# # Convert photoionization tables of all the levels in binary format
# #format_px(flag, px_org_fn, idx_fn, px_bin_fn)
# 
# # Read index of levels
# idx = pl.loadtxt(idx_fn, dtype='str')
# 
# for i, j, k, l, m, n, o in idx:
# 
#     try: 
#         #case of TOPBASE data
#         j.index('.')
#         j = pl.nan
#     except ValueError:
#         j = int(j)
# 
#     print(i, j, k, l, m, n, o, end=' ')
#     # Extract photoionization from binary file by direct access
#     e_ryd0, x_Mb0 = extract_px(px_bin_fn, i, k)
#     e_eV0 = e_ryd0 * 13.6 # energy in eV
# 
#     # Exclude null values
#     e_eV = pl.compress( x_Mb0 != 0., e_eV0)
#     x_Mb = pl.compress( x_Mb0 != 0., x_Mb0)
# 
#     # Keep sampling for high energy values where the variation follow Kramer's law
#     if isinstance(int(k)-j, int):
#         nkeep = int(k)-j
#     else:
#         nkeep = 10
# 
#     e_sel = e_eV[:-nkeep]
#     x_sel = x_Mb[:-nkeep]
# 
#     # Interpolate and smooth data
#     e_i = pl.linspace(min(e_sel), max(e_sel), 10000)
#     x_i = pl.interp(e_i, e_sel, x_sel)
#     x_is = smooth(x_i, 0.01)
# 
#     # Undersampling number
#     n_us = 150
#     n_us_keep = 13
# 
#     #e_us = e_i[::n_us]
#     #x_us = x_is[::n_us]
# 
#     e_us = pl.concatenate([e_i[::n_us], e_eV[int(k)-nkeep::n_us_keep]])
#     x_us = pl.concatenate([x_is[::n_us], x_Mb[int(k)-nkeep::n_us_keep]])
# 
#     if x_us.any() == 0.:
#         print('x_us = 0')
#         quit(1)
# 
#     # Conservation of area
# 
#     area = pl.trapz( x_Mb, e_eV)
#     #area_i = pl.trapz(x_i, e_i)
#     #area_is = pl.trapz(x_is, e_i)
#     area_us = pl.trapz(x_us, e_us)
# 
#     n_area = len(x_Mb)
#     n_area_us = len(x_us)
# 
# 
#     #print(format(area, '10.3e'), format(area_i, '10.3e'), format(area_is, '10.3e'), format(area_us, '10.3e'), end=' ')
#     print(format(area, '10.3e'), format(area_us, '10.3e'), format(n_area), format(n_area_us), end=' ')
# 
#     if len(x_Mb0) != len(x_Mb):
#         print('  ', len(x_Mb0)-len(x_Mb), " null values.")
#     else:
#         print('')
# 
# 
#     # Plot section
# 
#     pl.figure(figsize=(12, 6), dpi=100, facecolor='w', edgecolor='k')
#     pl.xlabel('E [eV]')
#     pl.xlim([0, 16])
#     pl.ylabel('Cross section [Mb]')
#     #pl.semilogx()
#     pl.semilogy()
#     pl.title(i+' '+str(j)+' '+k+' '+l+' '+m+' '+n+' '+o)
#     ax = pl.plot(e_eV, x_Mb, 'k', label='area = '+format(area, '11.3e')+' ('+format(n_area)+')')
#     pl.plot(e_us, x_us, 'r+-', label='area = '+format(area_us, '11.3e')+'('+format(n_area_us)+')')
#     pl.legend(frameon=False)
#     ofname=opath+l+m+n+'_'+o+'.png'
#     pl.savefig(ofname, dpi=100, format='png', orientation='landscape', bbox_inches='tight')
#     pl.close()
# 
# 
# #pl.show()

#######################################""