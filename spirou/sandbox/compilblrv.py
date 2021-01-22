import glob
from astropy.io import fits
import os
from tqdm import tqdm
import etienne_tools as et
import numpy as np

def compilblrv(obj_sci, obj_template = None, doplot = False, force = False):

    if doplot:
        # avoids conflicts if one uses ssh without graphic displays
        import matplotlib.pyplot as plt

    # pass just one object and we assume that the object is it's own template
    if obj_template is None:
        obj_template = obj_sci

    outname = 'tbllblrv_{0}_{1}.csv'.format(obj_sci, obj_template)

    # default keywords to be included in the
    keys = ['MJDATE', 'EXPTIME', 'AIRMASS', 'FILENAME', 'DATE-OBS', 'BERV', 'TAU_H2O', 'TAU_OTHE', 'ITE_RV', 'SYSTVELO',
            'TLPDVH2O','TLPDVOTR','CDBWAVE','OBJECT']
    keys = np.array(keys)

    if (not os.path.isfile(outname)) or force:
        # keys to be transfered to the fits table
        scifiles = np.array(glob.glob('lblrv/2*{0}_{1}.lblrv'.format(obj_sci, obj_template)))

        for i in tqdm(range(len(scifiles))):
            hdr = fits.getheader(scifiles[i])
            # if first file, we populate the columns
            if i == 0:
                tbl = dict()
                tbl['per_epoch_mean'] = np.zeros_like(scifiles, dtype = float)
                tbl['per_epoch_err'] = np.zeros_like(scifiles, dtype = float)
                tbl['per_epoch_DDV'] = np.zeros_like(scifiles, dtype = float)
                tbl['per_epoch_DDVRMS'] = np.zeros_like(scifiles, dtype = float)

                # adding keys from the input file
                for key in keys:
                    if key in hdr:
                        val = hdr[key]
                        if type(val) == str:
                            val =  ' '*100
                        tbl[key] = np.zeros_like(np.tile(val,len(scifiles)))

            for key in keys:
                # adding header value in output table
                if key in hdr:

                    tbl[key][i] = hdr[key]

            # read line file
            tbl_per_line_ini = fits.getdata(scifiles[i])

            # get lines in the proper domain
            if i == 0:
                g = (tbl_per_line_ini['WAVE_START'] > 900) * (tbl_per_line_ini['WAVE_START'] < 2500)
                rvs = np.zeros([len(scifiles),np.sum(g)])
                dvrms = np.zeros([len(scifiles), np.sum(g)])
                DDV = np.zeros([len(scifiles),np.sum(g)])
                DDVRMS = np.zeros([len(scifiles), np.sum(g)])

            # keep track in two big matrices of the rv and corresponding errors
            tbl_per_line = tbl_per_line_ini[g]
            rvs[i] = tbl_per_line['RV']
            dvrms[i] = tbl_per_line['DVRMS']

            DDV[i] = tbl_per_line['DDV']
            DDVRMS[i] = tbl_per_line['DDVRMS']

        # a line must be present >80% of the line
        valid_lines = np.nanmean(np.isfinite(rvs) * np.isfinite(dvrms),axis=0) > 0.8
        rvs = rvs[:,valid_lines]
        dvrms = dvrms[:,valid_lines]
        tbl_per_line_ini = tbl_per_line_ini[valid_lines]
        DDV = DDV[:,valid_lines]
        DDVRMS = DDVRMS[:,valid_lines]

        # constructing a per-epoch mean velocity
        for i in tqdm(range(len(scifiles))):
            guess,bulk_error  = et.odd_ratio_mean(rvs[i],dvrms[i])
            tbl['per_epoch_mean'][i] = guess
            tbl['per_epoch_err'][i] = bulk_error

            guess,bulk_error  = et.odd_ratio_mean(DDV[i],DDVRMS[i])

            tbl['per_epoch_DDV'][i] = guess
            tbl['per_epoch_DDVRMS'][i] = bulk_error

        # plot or not
        #if doplot:
        #    plt.errorbar(tbl['MJDATE']+0.2, tbl['per_epoch_mean'] - np.nanmean(tbl['per_epoch_mean']), fmt='.g', yerr=tbl['per_epoch_err'], alpha=0.5)

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


        # keeping track of the per-band RV measurements
        bands = ['Y','J','H','K']
        blue_end = [900,1150,1400,1900]
        red_end = [1150,1400,1900,2500]
        suffix = ['', '_0-2044', '_2044-4088',  '_1532-2556']


        rvs_matrix = np.zeros([len(scifiles),len(bands),len(suffix)])+np.nan
        err_matrix = np.zeros([len(scifiles),len(bands),len(suffix)])+np.nan

        # updating the table with per-band
        for i in tqdm(range(len(scifiles))):
            tmp_rv = rvs[i] - rv_per_line_model[i]
            tmp_err = dvrms[i]

            for iband in range(len(bands)):
                for reg in range(4):
                    if reg == 0:
                        g = (tbl_per_line_ini['WAVE_START'] > blue_end[iband]) * (tbl_per_line_ini['WAVE_START'] < red_end[iband])
                    if reg == 1:
                        g = (tbl_per_line_ini['WAVE_START'] > blue_end[iband]) * (tbl_per_line_ini['WAVE_START'] < red_end[iband]) * (tbl_per_line_ini['XPIX'] < 2044)
                    if reg == 2:
                        g = (tbl_per_line_ini['WAVE_START'] > blue_end[iband]) * (tbl_per_line_ini['WAVE_START'] < red_end[iband]) * (tbl_per_line_ini['XPIX'] > 2044)
                    if reg == 3:
                        g = (tbl_per_line_ini['WAVE_START'] > blue_end[iband]) * (tbl_per_line_ini['WAVE_START'] < red_end[iband])  * (tbl_per_line_ini['XPIX'] > 1532)  * (tbl_per_line_ini['XPIX'] < 2556)

                    guess,bulk_error  = et.odd_ratio_mean(tmp_rv[g],tmp_err[g])

                    rvs_matrix[i,iband,reg] = guess
                    err_matrix[i,iband,reg] = bulk_error

        for iband in range(len(bands)):
            for reg in range(4):
                tbl['per_epoch_mean_' + bands[iband]+suffix[reg]] = rvs_matrix[:,iband,reg]
                tbl['per_epoch_err_' + bands[iband]+suffix[reg]] = err_matrix[:,iband,reg]

        # plot or not
        if doplot:
            fig, ax = plt.subplots(nrows = 3, ncols = 1,sharex = True)
            ax[0].errorbar(tbl['MJDATE'], tbl['per_epoch_mean_J']-np.nanmedian(tbl['per_epoch_mean_J']) , fmt='.g', yerr=tbl['per_epoch_err_J'],alpha = 0.7,label = 'J')
            ax[0].errorbar(tbl['MJDATE'], tbl['per_epoch_mean_H']-np.nanmedian(tbl['per_epoch_mean_H']) , fmt='.k', yerr=tbl['per_epoch_err_H'],alpha = 0.7,label = 'H')
            ax[2].errorbar(tbl['MJDATE'],tbl['per_epoch_DDV'] ,fmt='.k', yerr=tbl['per_epoch_DDVRMS'], alpha = 0.7)

            diff_JH =  tbl['per_epoch_mean_J']-tbl['per_epoch_mean_H']
            ax[1].errorbar(tbl['MJDATE'],diff_JH - np.nanmedian(diff_JH) , fmt='.g', yerr=np.sqrt(tbl['per_epoch_err_J']**2+tbl['per_epoch_err_H']**2),alpha = 0.5,label = 'J')
            ax[0].legend()
            ax[0].set(xlabel = 'MJDATE', ylabel = 'RV [m/s]',title = 'H velocity')
            ax[1].set(xlabel = 'MJDATE', ylabel = 'RV [m/s]', title = 'J-H velo diff')
            ax[2].set(xlabel = 'MJDATE', ylabel = '2nd deriv')
            plt.show()

        et.td_convert(tbl).write(outname)
