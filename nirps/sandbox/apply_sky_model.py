from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from etienne_tools import lowpassfilter
from tqdm import tqdm
from etienne_tools import sigma
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.ndimage import binary_erosion

# retrieve A and B spectra with corresponding wavelengths
wave_A = fits.getdata('NIRPS_2022-11-25T12_00_56_233_pp_e2dsff_A_wavesol_ref_A.fits')
sp_A = np.array(fits.getdata('NIRPS_2022-12-06T01_30_31_491_pp_e2dsff_A.fits'),dtype = float)

wave_B = fits.getdata('NIRPS_2022-11-25T12_00_56_233_pp_e2dsff_B_wavesol_ref_B.fits')
sp_B = np.array(fits.getdata('NIRPS_2022-12-06T01_30_31_491_pp_e2dsff_B.fits'),dtype = float)

# retrived model
model_A = fits.getdata('sky_model_A.fits')
model_B = fits.getdata('sky_model_B.fits')
# retrieve region IDs
reg_id = fits.getdata('sky_reg_ID.fits')

weights = np.array(reg_id != 0,dtype = float).ravel()
weights2 = np.array(reg_id != 0,dtype = float).ravel()

for i in range(5):
    weights2 = binary_erosion(weights2,structure = np.ones(3))
    print(np.mean(weights2))
    weights+=np.array(weights2,dtype = float)

weights/=np.max(weights)
weights = weights.reshape(model_A.shape)

# get wavelength grid matching the models
wave1 = fits.getdata('sky_wavemap.fits')

nord = sp_B.shape[0]

sp_B_original = np.array(sp_B)

sp_B[~np.isfinite(sp_B)] = 0.0
for iord in range(nord):
    sp_B[iord] -= lowpassfilter(sp_B[iord], 51)
    for ite in range(2):
        nsig = sp_B[iord] / sigma(sp_B[iord])
        g = np.zeros_like(nsig)
        g[nsig > 3] = np.nan
        sp_B[iord] -= lowpassfilter(g + sp_B[iord], 101)

    # map sp_B onto reference grid
    sp_B[iord] = ius(wave_B[iord],sp_B[iord],ext=1,k=3)(wave1[iord])

sky_A = np.zeros_like(sp_A)
sky_B = np.zeros_like(sp_B)

for i in tqdm(np.unique(reg_id)[1:],leave = False):
    g = reg_id == i
    #fit amplitude in B fiber and apply to A
    amp = np.sum(sp_B[g]*model_B[g])/np.sum(model_B[g]**2)
    sky_A[g] = model_A[g]*amp*weights[g]
    sky_B[g] = model_B[g]*amp*weights[g]

# spline back to the file's wavelength grid
for iord in range(nord):
    sky_A[iord] = ius(wave1[iord],sky_A[iord],ext=1,k=3)(wave_A[iord])
    sky_B[iord] = ius(wave1[iord],sky_B[iord],ext=1,k=3)(wave_B[iord])

# some plots for fun
fig, ax = plt.subplots(nrows = 2, ncols = 1, sharex = True)
ax[0].plot(wave_A.ravel(),sp_A.ravel(),alpha = 0.5)
ax[0].plot(wave_A.ravel(),sp_A.ravel()-sky_A.ravel())
ax[1].plot(wave_B.ravel(),sp_B_original.ravel(),alpha = 0.5)
ax[1].plot(wave_B.ravel(),sp_B_original.ravel()-sky_B.ravel())
ax[1].plot(wave_B.ravel(),sky_B.ravel())
plt.show()

fits.writeto('sky_subtract_A.fits',sp_A - sky_A,overwrite = True)
fits.writeto('sky_subtract_B.fits',sp_B_original - sky_B,overwrite = True)
