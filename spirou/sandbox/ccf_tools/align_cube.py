import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from bisector import *
from astropy.time import Time
from ccf2rv import *
from scipy.interpolate import InterpolatedUnivariateSpline as ius


wavestart = 965
waveend = 2500
speed_of_light_kms = 299792.458
binvelo = 1.00

# work out number of wavelength points
flambda = np.log(waveend / wavestart)
nlambda = np.round((speed_of_light_kms / binvelo) * flambda)
# updating end wavelength slightly to have exactly 'step' km/s
waveend = np.exp(nlambda * (binvelo / speed_of_light_kms)) * wavestart
# get the wavegrid
index = np.arange(nlambda) / nlambda
wavegrid = wavestart * np.exp(index * np.log(waveend / wavestart))



# take a BigCube_s1d_v file from the DRS and cross-correlated all
# spectra to get a Template that is NOT RV shifted
cube_file = 'BigCube_s1d_TOI-1278_sc1d_v_file_AB.fits'

cube = fits.getdata(cube_file)

#cube = cube[:,159900:161000]
tbl_cube  = fits.getdata(cube_file, ext=1)

files_in_cube = np.array(tbl_cube['Filename'])
for i in range(len(files_in_cube)):
    files_in_cube[i] = files_in_cube[i].split('o')[0]


object = 'TOI-1278'
mask =  'gl846_neg'
exclude_orders = [-1]
# number of median-absolute deviations within an epoch to consider a point discrepant
tbl = get_object_rv(object,mask =mask,
                    method = 'template',force = True,
                    exclude_orders = exclude_orders,
                    snr_min = 20.0, sanitize = False,
                    dvmax_per_order = 3.0, bandpass = 'YJHK',
                    doplot = True, do_blacklist = True)



fig,ax = plt.subplots(nrows = 1, ncols =2,sharex = True, sharey = True)
ax[0].imshow(cube,vmin = 0.5,vmax = 1.5,aspect = 'auto')
ax[0].set(title = 'Before alignment', ylabel = 'Nth frame', xlabel = 'pixels')
ax[1].set(title = 'After alignment', ylabel = 'Nth frame', xlabel = 'pixels')

xpix = np.arange(cube.shape[1])
for i in range(cube.shape[0]):

    if  np.max(files_in_cube[i] == tbl['ODOMETER']  )   == True:
        rv =  tbl['RV'][np.where(files_in_cube[i] == tbl['ODOMETER'])[0][0]]


        #plt.plot(cube[i,:],color = 'red',alpha = 0.2)

        g = np.isfinite(cube[i,:])
        smask = ius(xpix,g,k=1)
        spline = ius(xpix[g],cube[i,:][g],k=3)

        mask = smask(xpix+rv)
        tmp = spline(xpix+rv)

        tmp[mask < 0.9] = np.nan
        cube[i, :] = tmp

        #plt.plot(cube[i, :],color = 'green',alpha = 0.2)
        print(rv)
    else:
        print('no odometer {0} in list'.format(files_in_cube[i]) )
        cube[i,:] = np.nan

ax[1].imshow(cube,vmin = 0.5,vmax = 1.5,aspect = 'auto')
ax[0].imshow(cube,vmin = 0.5,vmax = 1.5,aspect = 'auto')

plt.show()

sp = np.nanmedian(cube,axis=0)

fig, ax = plt.subplots(nrows = 2, ncols=1,sharex = True)
ax[0].plot(wavegrid, sp)
ax[1].plot(wavegrid, np.gradient(np.gradient(sp)))
plt.show()