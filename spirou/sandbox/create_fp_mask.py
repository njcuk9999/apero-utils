import numpy as np
from astropy.table import Table
import os

"""
Creation of a FP CCF mask from the cavity lenght properties
The are no arguments to pass, you need to edit the first few
lines in this code. Variables should be pretty self-explanatory.

Logic behind the code:
    We have a cavity length and we know that the FP peaks will
    range from ~10 000 to ~25 000. We take the cavity length and
    derive the corresponding wavelengths to be used as a CCF 
    mask. The process has to be done in an iterative way as 
    the exact wavelength of a line affects the cavity length. The
    computation coverges in a few steps. The CCF mask is saved 
    with a blue/red edge that are -0.5*dv and 0.5*dv from the central
    line position. Peak weights are all set to 1.0. 
"""

cavity_file = 'apero-drs/apero/data/spirou/calib/cavity_length_ll_fit.dat'
outname = 'smart_fp_mask.mas'
dv = 1.0 # width of the CCF
c = 299792.458 # velocity of light

# define cuts in wavelenght
cut = [970,2500]

# read the cavity length as used by the DRS
print('Read {0}'.format(cavity_file))
f = open(cavity_file, 'r')
tbl = f.read()
f.close()

# transform into a numpy array
fit = np.array(tbl.split('\n'))
fit = fit[fit != '']
fit = np.array(fit,dtype = float)

# start with a broader range of FP N values and clip later
# down
n_fp_fpeak = np.arange(9000,27000)
# placeholder for wavelength, needs to be iterated-on with
# the cavity length polynomial
# the starting wavelength is the midpoint between
# the extremities of the domain
wave_fp_peak = np.ones(len(n_fp_fpeak))*np.mean(cut)

# we perform the following loop
#
# Take the wavelength and derive a cavity length
# Take the Nth peak and derive a line wavelength from the
# cavity length
# take the new wavelength and dereive a new cavity length
# find an update wavelength for the Nth peak
#  .... converge down to a 1e-9 error in wavelength

dwave = 1e9
while abs(dwave) > 1e-9:
    # keep track of the central line to check convergence
    prev = wave_fp_peak[len(wave_fp_peak)//2]

    # derive a new wavelength for each fp peak
    wave_fp_peak = np.polyval(fit, wave_fp_peak) * 2 / n_fp_fpeak

    # check convergnce
    dwave = prev - wave_fp_peak[len(wave_fp_peak)//2]

# keep lines within the domain
keep = (wave_fp_peak > cut[0]) & (wave_fp_peak < cut[1])

wave_fp_peak = wave_fp_peak[keep]
n_fp_fpeak = n_fp_fpeak[keep]

print('Write {0}'.format(outname))

# export everything into the mas file
f = open(outname, 'w')
for i in range(len(wave_fp_peak)):
    f.write('{0:20.10f}{1:20.10f}{2:20.10f}\n'.format(wave_fp_peak[i]*(1-0.5*dv/c),wave_fp_peak[i]*(1+0.5*dv/c),1.0))
f.close()

os.system('say "ça marche les amis!"&')