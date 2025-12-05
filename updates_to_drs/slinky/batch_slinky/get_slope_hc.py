from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
from etienne_tools import snail, sigma
from scipy.constants import c

def odd_ratio_linfit(x, y, yerr):
    g = np.isfinite(y + yerr + x)
    x = x[g]
    y = y[g]
    yerr = yerr[g]

    w = np.ones(len(x))

    sum = 1.0
    sum0 = 0.0
    while np.abs(sum0 - sum) > 1e-6:
        sum0 = np.sum(w)
        fit, sig = np.polyfit(x, y, 1, w=w / yerr, cov=True)
        errfit = np.sqrt(np.diag(sig))

        res = (y - np.polyval(fit, x)) / yerr
        p1 = np.exp(-0.5 * res ** 2)
        p2 = 1e-6

        w = p1 / (p1 + p2)
        # print(fit, sum,sum0)
        sum = np.sum(w)

    return fit, errfit

path_to_calibs = '/Volumes/courlan/lbl_NIRPS_HE/calib/'

all_calibs = glob.glob(f'{path_to_calibs}/*_wave_hclines_*.fits')

tbl_params = Table()
tbl_params['WAVEFILE'] = np.zeros(len(all_calibs),dtype = 'U999')
tbl_params['mjd'] = np.zeros(len(all_calibs),dtype = float)
tbl_params['slope'] = np.zeros_like(all_calibs, dtype=float)
tbl_params['offset_1600'] = np.zeros_like(all_calibs, dtype=float)
tbl_params['slope_err'] = np.zeros_like(all_calibs, dtype=float)
tbl_params['offset_1600_err'] = np.zeros_like(all_calibs, dtype=float)
tbl_params['file'] = all_calibs

for ifile in snail(range(len(all_calibs))):
    h = fits.getheader(all_calibs[ifile])
    tbl_params['WAVEFILE'][ifile] = h['WAVEFILE']
    tbl = Table.read(all_calibs[ifile])

    wave = np.array(tbl['WAVE_MEAS'].data.data)# / 1000
    wave -= 1600
    nsig = tbl['NSIG'].data.data

    off = (tbl['WAVE_MEAS'] / tbl['WAVE_REF'] - 1).data.data * c
    sig1 = sigma(off * nsig)
    sig = sig1 / nsig

    x = wave
    y = off
    yerr = sig

    keep = np.isfinite(x + y + yerr)
    x = x[keep]
    y = y[keep]
    yerr = yerr[keep]
    nsig = y/yerr
    keep = (np.abs(nsig) < 8)
    x = x[keep]
    y = y[keep]
    yerr = yerr[keep]

    fit, errfit = odd_ratio_linfit(x, y, yerr)

    tbl_params['slope'][ifile] = fit[0]
    tbl_params['offset_1600'][ifile] = fit[1]
    tbl_params['slope_err'][ifile] = errfit[0]
    tbl_params['offset_1600_err'][ifile] = errfit[1]
    tbl_params['mjd'][ifile] = h['MJDMID']

tbl_params = tbl_params[np.argsort(tbl_params['mjd'])]

tbl_params.write(path_to_calibs+'cavity_slope.csv')