import numpy as np
from astropy.io import fits
from astropy.table import Table
import glob
# import matplotlib.pyplot as plt
import etienne_tools as et
from astropy import constants
from scipy.interpolate import InterpolatedUnivariateSpline as ius
import os
from tqdm import tqdm


def mk_harps_template(files, outname, obj, teff):
    # files -> output from a glob.glob with all s1d files
    # outname -> full name of output template file
    # obj -> for proper header handling, name of object
    # teff -> for proper header handling, temperature in K

    grid = et.get_magic_grid(wave0=378, wave1=692, dv_grid=500)

    cube = np.zeros([len(grid), len(files)])

    tbl2 = Table()
    tbl2['FILES'] = files
    print('load all files')
    for i in tqdm(range(len(files))):
        sp, hdr = fits.getdata(files[i], header=True)

        berv = hdr['HIERARCH ESO DRS BERV']
        wave1 = (np.arange(len(sp)) * hdr['CDELT1'] + hdr['CRVAL1']) / 10.0
        sp /= np.nanmedian(sp[(wave1 > 600) * (wave1 < 630)])
        sp[sp == 0] = np.nan

        g = np.isfinite(sp)
        tmp1 = ius(wave1[g], sp[g], k=3, ext=1)(grid)
        tmp2 = ius(wave1, (np.isfinite(sp) * 1.0), k=1, ext=1)(grid)
        tmp1[tmp2 < 0.99] = np.nan
        cube[:, i] = tmp1

    cube[cube == 0] = np.nan
    nfinite = np.mean(np.isfinite(cube), axis=1)
    med = np.nanmedian(cube, axis=1)
    med[nfinite < 0.5] = np.nan

    tbl = Table()
    tbl['wavelength'] = grid
    tbl['flux'] = med
    tbl['eflux'] = np.zeros_like(med)
    tbl['rms'] = np.nan + np.zeros_like(med)

    hdr['OBJECT'] = obj  # 'PROXIMA'
    hdr['OBJTEMP'] = teff  # 3042

    hdu0 = fits.PrimaryHDU()
    hdu1 = fits.BinTableHDU(tbl)
    hdu2 = fits.BinTableHDU(tbl2)

    hdu0.header = et.harps2spirou(hdr)
    hdu1.header['OBJECT'] = obj
    hdu1.header['OBJTEMP'] = hdr['OBJTEMP']

    # convert back from dictionnary to table and save
    new_hdul = fits.HDUList([hdu0, hdu1, hdu2])
    new_hdul.writeto(outname, overwrite=True)

