import etienne_tools as et
from astropy.io import fits
import glob
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

files = glob.glob('tellurics/TRAPPIST-1/*s1d*.fits')

wave = Table.read(files[0])['wavelength']

cube = np.zeros([len(files),len(wave)])

for i in range(len(files)):
    print(i,len(files))
    flux = Table.read(files[i])['flux']
    hdr = fits.getheader(files[i],ext=1)
    keep = np.isfinite(flux)
    wave2 = et.doppler(wave, -1000 * hdr['BERV'])
    tmp = ius(wave2[keep], flux[keep], k=1, ext=3)(wave)
    mask = ius(wave2, np.isfinite(flux), k=1, ext=0)(wave)
    tmp[mask < 0.99] = np.nan
    cube[i] = tmp

cube_hf = np.array(cube)

for i in range(len(files)):
    print(i, len(files))
    cube_hf[i] -= et.lowpassfilter(cube_hf[i], width=201)

g = (wave > 1500) * (wave < 1700)


ref = np.nanmedian(cube_hf, axis=0)
for ite in range(3):
    for i in range(len(files)):
        gg = g * np.isfinite(cube_hf[i]) * np.isfinite(ref)
        amp = np.nansum(cube_hf[i][gg] * ref[gg]) / np.nansum(ref[gg] ** 2)
        cube_hf[i] /= amp
        cube[i] /= amp
        print(ite, i, len(files), amp)
ref = np.nanmedian(cube_hf, axis=0)

amp = np.nanmedian(cube[:,g],axis=1)
lowf = et.lowpassfilter(cube[np.argsort(amp)[len(files)//4]])

plt.plot(lowf)
plt.plot(ref)
plt.plot(ref+lowf)
plt.show()

template = ref+lowf
template/=np.nanmedian(template[g])
tbl = Table.read(files[i])
tbl['flux'] = template

sp, hdr = fits.getdata(files[0], header=True)
hdu1 = fits.PrimaryHDU()

hdr0 = hdu1.header

keys = np.array([key for key in hdr.keys()])
for key in keys[10:]:
    if key not in hdr0.keys():
        hdr0[key] = hdr[key]
hdu1.header = hdr0
# convert back from dictionnary to table and save
hdu2 = fits.BinTableHDU(tbl)

new_hdul = fits.HDUList([hdu1, hdu2])
new_hdul.writeto('Template_s1d_TRAPPIST-1_sc1d_v_file_AB.fits', overwrite=True)


#fits.writeto('test1.fits', cube, overwrite=True)
#fits.writeto('test2.fits', cube_hf, overwrite=True)