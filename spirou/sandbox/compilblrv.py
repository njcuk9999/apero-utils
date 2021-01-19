import glob
from astropy.io import fits
import os
import matplotlib.pyplot as plt
from tqdm import tqdm
import etienne_tools as et
from astropy.table import Table
import numpy as np


#obj_sci, obj_template = 'GL436','GL436'

def compilblrv(obj_sci, obj_template, doplot = False, force = False):
    outname = 'tbllblrv_{0}_{1}.csv'.format(obj_sci, obj_template)

    # default keywords to be included in the
    keys = ['MJDATE', 'EXPTIME', 'AIRMASS', 'FILENAME', 'DATE-OBS', 'BERV', 'TAU_H2O', 'TAU_OTHE', 'ITE_RV', 'SYSTVELO',
            'TLPDVH2O','TLPDVOTR','CDBWAVE']
    keys = np.array(keys)

    if (not os.path.isfile(outname)) or force:
        # keys to be transfered to the fits table
        scifiles = np.array(glob.glob('lblrv/2*{0}_{1}.lblrv'.format(obj_sci, obj_template)))

        for i in tqdm(range(len(scifiles))):
            hdr = fits.getheader(scifiles[i])
            # if first file, we populate the columns
            if i == 0:
                tbl = Table()
                tbl['per_epoch_mean'] = np.zeros_like(scifiles, dtype = float)
                tbl['per_epoch_err'] = np.zeros_like(scifiles, dtype = float)

                # adding keys from the input file
                for key in keys:
                    if key in hdr:
                        tbl[key] = np.zeros_like(np.tile(hdr[key],len(scifiles)))

            for key in keys:
                # adding header value in output table
                if key in hdr:
                    tbl[i][key] = hdr[key]

            # read line file
            tbl_per_line_ini = fits.getdata(scifiles[i])

            # get lines in the proper domain
            if i == 0:
                g = (tbl_per_line_ini['WAVE_START'] > 900) * (tbl_per_line_ini['WAVE_START'] < 2500)
                rvs = np.zeros([len(scifiles),np.sum(g)])
                dvrms = np.zeros([len(scifiles), np.sum(g)])

            # keep track in two big matrices of the rv and corresponding errors
            tbl_per_line = tbl_per_line_ini[g]
            rvs[i] = tbl_per_line['RV']
            dvrms[i] = tbl_per_line['DVRMS']

        # a line must be present >80% of the line
        valid_lines = np.nanmean(np.isfinite(dvrms),axis=0) > 0.8
        rvs = rvs[:,valid_lines]
        dvrms = dvrms[:,valid_lines]
        tbl_per_line_ini = tbl_per_line_ini[valid_lines]

        # constructing a per-epoch mean velocity
        for i in tqdm(range(len(scifiles))):
            guess,bulk_error  = et.odd_ratio_mean(rvs[i],dvrms[i])
            tbl['per_epoch_mean'][i] = guess
            tbl['per_epoch_err'][i] = bulk_error

        # plot or not
        if doplot:
            plt.errorbar(tbl['MJDATE']+0.2, tbl['per_epoch_mean'] - np.nanmean(tbl['per_epoch_mean']), fmt='.g', yerr=tbl['per_epoch_err'], alpha=0.5)

        # de-biasing line. Matrix that contains a replicated 2d version of the per-epoch mean
        rv_per_epoch_model = np.reshape(np.repeat(tbl['per_epoch_mean'],rvs.shape[1]),rvs.shape)

        # line-by-line mean position
        per_line_mean = np.zeros(rvs.shape[1])
        per_line_error = np.zeros(rvs.shape[1])

        # computing the per-line bias
        for i in tqdm(range(len(per_line_mean))):
            guess,bulk_error  = et.odd_ratio_mean(rvs[:,i]-rv_per_epoch_model[:,i],dvrms[:,i])
            per_line_mean[i] = guess
            per_line_error[i] = bulk_error

        # we keep the per-line mean to zero
        per_line_mean -= (et.odd_ratio_mean(per_line_mean,per_line_error))[0]

        # constructing a 2d model of line biases
        rv_per_line_model = np.reshape(np.tile(per_line_mean,rvs.shape[0]),rvs.shape)

        # updating the table
        for i in tqdm(range(len(scifiles))):
            guess,bulk_error  = et.odd_ratio_mean(rvs[i]-rv_per_line_model[i],dvrms[i])
            tbl['per_epoch_mean'][i] = guess
            tbl['per_epoch_err'][i] = bulk_error

        # plot or not
        if doplot:
            plt.errorbar(tbl['MJDATE'], tbl['per_epoch_mean']-np.nanmean(tbl['per_epoch_mean']) , fmt='.r', yerr=tbl['per_epoch_err'],alpha = 0.5)
            plt.show()

        # keeping track of the per-band RV measurements
        bands = ['Y','J','H','K']
        blue_end = [900,1150,1400,1900]
        red_end = [1150,1400,1900,2500]

        # updating the table with per-band
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



