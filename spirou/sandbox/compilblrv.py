import glob
from astropy.io import fits
import os
import matplotlib.pyplot as plt
#from low_pass_filter  import *
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from astropy.table import Table
from scipy.special import erf
from tqdm import tqdm
from scipy import constants
from scipy import stats
from scipy.stats import chi2
import etienne_tools as et
from astropy.table import Table
import numpy as np



obj_sci, obj_template = 'TOI-1278','GL846'

outname = 'tbllblrv_{0}_{1}.csv'.format(obj_sci, obj_template)
keys = ['MJDATE', 'EXPTIME', 'AIRMASS', 'FILENAME', 'DATE-OBS', 'BERV', 'TAU_H2O', 'TAU_OTHE', 'ITE_RV', 'SYSTVELO']
doplot = True


if not os.path.isfile(outname):
    # keys to be transfered to the fits table
    scifiles = np.array(glob.glob('xdrv/2*{0}_{1}.xdrv'.format(obj_sci, obj_template)))

    for i in tqdm(range(len(scifiles))):
        hdr = fits.getheader(scifiles[i])

        if i == 0:
            tbl = Table()
            tbl['per_epoch_mean'] = np.zeros_like(scifiles, dtype = float)
            tbl['per_epoch_err'] = np.zeros_like(scifiles, dtype = float)

            for key in keys:
                if key in hdr:
                    tbl[key] = np.zeros_like(np.tile(hdr[key],len(scifiles)))

        for key in keys:
            if key in hdr:
                tbl[i][key] = hdr[key]

        tbl_per_line_ini = fits.getdata(scifiles[i])

        if i == 0:
            g = (tbl_per_line_ini['WAVE_START'] > 900) * (tbl_per_line_ini['WAVE_START'] < 2500)
            rvs = np.zeros([len(scifiles),np.sum(g)])
            dvrms = np.zeros([len(scifiles), np.sum(g)])

        tbl_per_line = tbl_per_line_ini[g]
        rvs[i] = tbl_per_line['RV']
        dvrms[i] = tbl_per_line['DVRMS']

    valid_lines = np.nanmean(np.isfinite(dvrms),axis=0) > 0.8
    rvs = rvs[:,valid_lines]
    dvrms = dvrms[:,valid_lines]
    tbl_per_line_ini = tbl_per_line_ini[valid_lines]

    for i in tqdm(range(len(scifiles))):
        guess,bulk_error  = et.odd_ratio_mean(rvs[i],dvrms[i])
        tbl['per_epoch_mean'][i] = guess
        tbl['per_epoch_err'][i] = bulk_error

    if doplot:
        plt.errorbar(tbl['MJDATE']+0.2, tbl['per_epoch_mean'] - np.nanmean(tbl['per_epoch_mean']), fmt='.g', yerr=tbl['per_epoch_err'], alpha=0.5)

    rv_per_epoch_model = np.reshape(np.repeat(tbl['per_epoch_mean'],rvs.shape[1]),rvs.shape)

    per_line_mean = np.zeros(rvs.shape[1])
    per_line_error = np.zeros(rvs.shape[1])

    for i in tqdm(range(len(per_line_mean))):
        guess,bulk_error  = et.odd_ratio_mean(rvs[:,i]-rv_per_epoch_model[:,i],dvrms[:,i])
        per_line_mean[i] = guess
        per_line_error[i] = bulk_error
    # we keep the per-line mean to zero
    per_line_mean -= (et.odd_ratio_mean(per_line_mean,per_line_error))[0]

    rv_per_line_model = np.reshape(np.tile(per_line_mean,rvs.shape[0]),rvs.shape)

    for i in tqdm(range(len(scifiles))):
        guess,bulk_error  = et.odd_ratio_mean(rvs[i]-rv_per_line_model[i],dvrms[i])
        tbl['per_epoch_mean'][i] = guess
        tbl['per_epoch_err'][i] = bulk_error

    if doplot:
        plt.errorbar(tbl['MJDATE'], tbl['per_epoch_mean']-np.nanmean(tbl['per_epoch_mean']) , fmt='.r', yerr=tbl['per_epoch_err'],alpha = 0.5)
        plt.show()

    bands = ['Y','J','H','K']
    blue_end = [900,1150,1400,1900]
    red_end = [1150,1400,1900,2500]

    for i in tqdm(range(len(scifiles))):
        tmp_rv = rvs[i] - rv_per_line_model[i]
        tmp_err = dvrms[i]

        for iband in range(len(bands)):
            g = (tbl_per_line_ini['WAVE_START'] > blue_end[iband]) * (tbl_per_line_ini['WAVE_START'] < red_end[iband])

            if i == 0:
                tbl['per_epoch_mean_'+bands[iband]] = np.zeros_like(scifiles,dtype = float)
                tbl['per_epoch_err_'+bands[iband]] = np.zeros_like(scifiles,dtype = float)

            guess,bulk_error  = et.odd_ratio_mean(tmp_rv[g],tmp_err[g])

            tbl['per_epoch_mean_' + bands[iband]][i] = guess
            tbl['per_epoch_err_' + bands[iband]][i] = bulk_error

    tbl.write(outname)
