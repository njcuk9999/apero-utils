from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
from etienne_tools import snail, sigma
from scipy.constants import c
import os

def gp_project(x,y,yerr,ww = 1e-1,xmin = 980,xmax = 1850, npts = 100000):

    xv = np.linspace(xmin,xmax,npts)
    weights = np.zeros_like(xv)
    yv = np.zeros_like(xv)

    for i in snail(range(len(x)), leave = False):
        dd =(xv-x[i])/ww

        w2 = np.exp(-0.5*dd**2)/yerr[i]**2
        weights+=w2
        yv+=w2*y[i]

    yv/=weights

    return xv,yv

path_to_calibs = '/Volumes/courlan/lbl_NIRPS_HE/calib/'

all_calibs = glob.glob(f'{path_to_calibs}/*_wave_fplines_*.fits')


for ifile,file in (enumerate(all_calibs)):
    print(ifile, len(all_calibs))
    h = fits.getheader(file)
    outname = path_to_calibs+h['WAVEFILE'].replace('.fits','_sliky.csv')
    if os.path.isfile(outname):
        print('already done, file {} of {}'.format(ifile,len(all_calibs)))
        continue

    tbl = Table.read(file)
    wave = tbl['WAVE_MEAS'].data.data #/ 1000
    nsig = tbl['NSIG'].data.data

    off = (tbl['WAVE_MEAS'] / tbl['WAVE_REF'] - 1).data.data * c
    sig1 = sigma(off * nsig)
    sig = sig1 / nsig

    ord= np.argsort(wave)
    wave = wave[ord]
    off = off[ord]
    sig = sig[ord]

    x = wave
    y = off
    yerr = sig
    nsig = y/yerr

    valid = np.isfinite(x + y + yerr) & (np.abs(nsig) < 8)
    x = x[valid]
    y = y[valid]
    yerr = yerr[valid]
    nsig = nsig[valid]

    x_pred, pred = gp_project(x, y, yerr, ww=1e-1, xmin=960, xmax=1850, npts=100000)

    plt.errorbar(x, y, yerr=yerr, fmt='.', alpha=0.5)
    plt.plot(x_pred, pred)
    plt.show()


    tbl = Table()
    tbl['wavelength'] = x_pred
    tbl['corr'] = pred

    tbl.write(outname, format = 'ascii.csv', overwrite = True)