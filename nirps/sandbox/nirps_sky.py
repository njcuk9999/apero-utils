from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from etienne_tools import lowpassfilter
import glob
from tqdm import tqdm
from etienne_tools import sigma
from scipy.ndimage import binary_erosion, binary_dilation
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from astropy.table import Table

# This code loads and processes sky data from NIRPS to identify and mask out absorption lines.
# This is useful for removing sky lines that could affect the accuracy of spectral analysis.

# get wavelength grid for A spectra
waveref_A = fits.getdata('NIRPS_2022-11-25T12_00_56_233_pp_e2dsff_A_wavesol_ref_A.fits')
# get reference A files
all_sky_A = np.array(glob.glob('sky/*_pp_e2dsff_A.fits'))

# get wavelength grid for B spectra
waveref_B = fits.getdata('NIRPS_2022-11-25T12_00_56_233_pp_e2dsff_B_wavesol_ref_B.fits')
# get corresponding reference B files
all_sky_B = np.array([f.replace('_A.','_B.') for f in all_sky_A])

# find number of orders
nord = fits.getdata(all_sky_A[0]).shape[0]

# place holder for the A and B fiber sky 'cubes'
all_sky_map_A = np.zeros([len(all_sky_A),len(waveref_A)*4088])
all_sky_map_B = np.zeros_like(all_sky_map_A)

# loading all sky data
for i in tqdm(range(len(all_sky_A))):
    sky_A = np.array(fits.getdata(all_sky_A[i]),dtype = float)
    # avoid NaNs for the spline later
    sky_A[~np.isfinite(sky_A)] = 0.0

    sky_B = np.array(fits.getdata(all_sky_B[i]),dtype = float)
    # avoid NaNs for the spline later
    sky_B[~np.isfinite(sky_B)] = 0.0

    if i == 0:
        # this should be the master grid
        wave1 = waveref_A

    for iord in range(nord):
        # high-pass the sky
        sky_A[iord] -= lowpassfilter(sky_A[iord],101)
        sky_B[iord] -= lowpassfilter(sky_B[iord],101)

        # spline onto the master grid
        sky_A[iord] = ius(waveref_A[iord],sky_A[iord],ext=1,k=3)(wave1[iord])
        sky_B[iord] = ius(waveref_B[iord],sky_B[iord],ext=1,k=3)(wave1[iord])

    # expressed as a ravelled e2ds
    all_sky_map_A[i] = sky_A.ravel()
    all_sky_map_B[i] = sky_B.ravel()

# get ravel-e2ds of the wavelegnth grid for dilation stuff later
waveref2 = wave1.ravel()

# get median sky spectrum, only used for identification of lines
v = np.nanmedian(all_sky_map_A,axis=0)
v[~np.isfinite(v)] = 0

# find positive excursions in sky signal
nsig = v/sigma(v)
# lines are >5 sigma positive excursion
line = np.array(nsig>5,dtype = int)
# erode features that are too narrow
line = binary_erosion(line,structure = np.ones(5))
# dilate to get wings of lines
line = binary_dilation(line,structure = np.ones(27))


reg_id = (np.cumsum(line!=np.roll(line,1)))
reg_id[(reg_id % 2) == 0] = 0
reg_id[reg_id != 0] = (reg_id[reg_id != 0]+1)//2

# put the line mask onto a magic grid to avoid errors at order overlaps
magic_grid = np.array(Table.read('NIRPS_2022-12-02T00_59_50_241_pp_s1d_v_A.fits')['wavelength'])
magic_mask = np.zeros_like(magic_grid, dtype = bool)


# put everything in the magic grid referential
for i in np.unique(reg_id)[1:]:
    g = (i == reg_id)
    magic_mask[ (magic_grid>np.min(waveref2[g]))*(magic_grid<np.max(waveref2[g]))] = True

# now in the space of the magic grid
reg_id_magic = (np.cumsum(magic_mask!=np.roll(magic_mask,1)))
reg_id_magic[(reg_id_magic % 2) == 0] = 0
reg_id_magic[reg_id_magic != 0] = (reg_id_magic[reg_id_magic != 0]+1)//2

# fill the map with unique values and common ID for overlapping orders
reg_id = np.zeros_like(reg_id, dtype = int)
for ureg in np.unique(reg_id_magic)[1:]:
    wave_min = np.min(magic_grid[reg_id_magic == ureg])
    wave_max = np.max(magic_grid[reg_id_magic == ureg])

    g = (waveref2>wave_min)*(waveref2<wave_max)
    reg_id[g] = ureg


# plots that you can remove, for debug
plt.plot(waveref2, v)
plt.plot(waveref2[reg_id !=0], v[reg_id !=0],'g.',alpha =  0.3)

for reg in np.unique(reg_id):
    if reg == 0:
        continue

    gg = reg_id == reg
    plt.text(np.mean(waveref2[gg]),0,'{}'.format(reg))
plt.show()

# construct a model with all lines normalized to a median of 1 in fiber A
model_A = np.zeros_like(v)
model_B = np.zeros_like(v)

for ii in tqdm(np.unique(reg_id)[1:]):
    if ii == 145:
        doplot = True
    else:
        doplot = False
    if doplot:
        fig, ax = plt.subplots(nrows = 2,ncols=1,sharex = True)
    gg = reg_id.ravel() == ii

    all_A = []
    all_B = []
    for i in range(len(all_sky_A)):
        tmp_A = all_sky_map_A[i][gg]
        tmp_B = all_sky_map_B[i][gg]

        amp = np.nansum(tmp_A)
        all_A = np.append(all_A, tmp_A/amp)
        all_B = np.append(all_B, tmp_B/amp)

        if doplot:
            ax[0].plot(waveref2[gg],tmp_A/amp,alpha = 0.5,color = 'orange')
            ax[1].plot(waveref2[gg],tmp_B/amp,alpha = 0.5,color = 'orange')

    med_A = np.nanmedian(all_A.reshape(len(all_sky_A),len(tmp_A)),axis=0)
    med_B = np.nanmedian(all_B.reshape(len(all_sky_B),len(tmp_B)),axis=0)

    if doplot:
        ax[0].plot(waveref2[gg], med_A, color='blue')
        ax[1].plot(waveref2[gg], med_B, color='blue')

    model_A[gg] = med_A
    model_B[gg] = med_B
    if doplot:
        plt.show()

model_A = model_A.reshape(sky_A.shape)
model_B = model_B.reshape(sky_B.shape)
reg_id = reg_id.reshape(sky_A.shape)

fits.writeto('sky_model_A.fits',model_A,overwrite = True)
fits.writeto('sky_model_B.fits',model_B,overwrite = True)
fits.writeto('sky_reg_ID.fits',reg_id,overwrite = True)
fits.writeto('sky_wavemap.fits',wave1,overwrite = True)

