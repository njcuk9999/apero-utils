import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as ius
import etienne_tools as et
from bisector import bisector


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# input file for which we want an absolute RV
file = '2585295o_pp_s1d_w_tcorr_AB.fits'
#file = '2585275o_pp_s1d_w_tcorr_AB.fits'
#file = 'Template_s1d_GL699_sc1d_v_file_AB.fits'

doplot = True

teff = 3500
logg = 5.50

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
# logg = 4.0, 4.5, 5.0, 5.5
#teff = 2700 to 6000 in 100 K steps
# dowload all masks here :
# www.astro.umontreal.ca/~artigau/all_model_masks.tar.gz

import os
os.system('wget www.astro.umontreal.ca/~artigau/all_model_masks.tar.gz')
"""

# get the mask corresponding to the temperature of your object
mask = Table.read('HiResFITS/lte0{0}-{1:.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.csv'.format(teff,logg))

# read s1d or template as a table
tbl = Table.read(file)
hdr = fits.getheader(file,ext=1)

if 'Template' not in file:
    # read header to get BERV
    hdr = fits.getheader(file,ext=1)
    berv = hdr['BERV']*1000
else:
    print('This is a template, no berv adjustment')
    berv = 0

# wavelength and flux from table. Converting to numpy array to speed-up
wave = np.array(tbl['wavelength'])
sp = np.array((tbl['flux']))

# save for the mask
wave_mask = np.array(mask['wavelength'])
weight_mask = np.array(mask['weight'])

# speedier than keep all the mask lines
good_mask = (wave_mask>np.min(wave))*(wave_mask<np.max(wave))
weight_mask,wave_mask = weight_mask[good_mask], wave_mask[good_mask]

# valid
valid_sp = np.isfinite(sp)



spl_valid = ius(wave,valid_sp,k=1,ext=1)

good_line = np.ones_like(wave_mask,dtype = bool)

# coarse scan to reject lines that are affected at some point by a gap in the spectrum
dv = np.arange(-200,200,3,dtype = float)*1000

for i in range(len(dv)):
    good_line*=(spl_valid(et.doppler(wave_mask,-dv[i]+berv))>0.99)

# rejecting lines bad at least one
wave_mask = wave_mask[good_line]
weight_mask = weight_mask[good_line]

spl = ius(wave[valid_sp],sp[valid_sp],k=1,ext=1)

# finer scan
dv = np.arange(-200,200,.5,dtype = float)*1000
ccf = np.zeros_like(dv)

for i in range(len(dv)):
    ccf[i] = np.sum(weight_mask*spl(et.doppler(wave_mask,-dv[i]+berv)))
ccf/=np.nanmedian(ccf)

#cen, ew, amp,expo, zp, slope
p0 =dv[np.argmin(ccf)],5000,np.min(ccf)-np.nanmedian(ccf),1.4,np.nanmedian(ccf),0
fit = et.fit_super_gauss(dv,ccf,p0)

if doplot:
    plt.plot(dv/1000,ccf)
    plt.xlabel( 'velocity [km/s]')
    plt.ylabel('contrast')
    plt.plot(dv/1000,et.super_gauss(dv,*fit))
    plt.title(hdr['OBJECT'])
    plt.show()

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Object : {}'.format(hdr['OBJECT']))
print('~~~~ Super-gaussian fit properties ~~~~')
print(' velocity : {:.2f} km/s'.format(fit[0]/1000))


fwhm = 2*fit[1]*(-2*np.log(0.5))**(1/fit[3])
print(' FWHM : {:.2f} km/s'.format(fwhm/1000))


# get bisector properties
depth, bisector_position, width_ccf = bisector(dv/1000, ccf,doplot=doplot,low_high_cut=0.2 )

v0 = ius(depth[np.argsort(depth)],bisector_position[np.argsort(depth)])(0.5)
fwhm = ius(depth[np.argsort(depth)],width_ccf[np.argsort(depth)])(0.5)

print('~~~~ bisector properties ~~~~')
print('Radial velocity : {0:.3f} km/s'.format(v0))
print('FWHM : {0:.3f} km/s'.format(fwhm))
