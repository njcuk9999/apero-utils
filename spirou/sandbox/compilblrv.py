import glob
from astropy.io import fits
import os
from tqdm import tqdm
import etienne_tools as et
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.time import Time

def compilblrv(obj_sci, obj_template = None, doplot = True, force = True, common_weights = False,
               get_cumul_plot = False, fcut = 0.95):

    """
    obj_sci = 'TRAPPIST-1'
    obj_template = None
    doplot = False
    force = True
    common_weights = False
    get_cumul_plot = False
    fcut = 0.95
    """

    if doplot or get_cumul_plot:
        # avoids conflicts if one uses ssh without graphic displays
        import matplotlib.pyplot as plt

    # pass just one object and we assume that the object is it's own template
    if obj_template is None:
        obj_template = obj_sci

    suffix = ''

    outname = 'tbllblrv{2}_{0}_{1}.csv'.format(obj_sci, obj_template,suffix)

    # default keywords to be included in the compilation of per-object RVs
    keys = ['MJDATE', 'EXPTIME', 'AIRMASS', 'FILENAME', 'DATE-OBS', 'BERV', 'TAU_H2O', 'TAU_OTHE', 'ITE_RV', 'SYSTVELO',
            'TLPDVH2O','TLPDVOTR','CDBWAVE','OBJECT','SBRHB1_P','SBRHB2_P','SBCDEN_P','SNRGOAL',
            'EXTSN035']
    keys = np.array(keys)

    if (not os.path.isfile(outname)) or force:
        # keys to be transfered to the fits table
        scifiles = np.array(glob.glob('lblrv/2*{0}_{1}_lbl.fits'.format(obj_sci, obj_template)))
        scifiles = scifiles[np.argsort(scifiles)]

        if get_cumul_plot:
            fix, ax = plt.subplots(nrows = 2, ncols =1)
        for i in tqdm(range(len(scifiles))):
            hdr = fits.getheader(scifiles[i])
            # if first file, we populate the columns
            if i == 0:
                tbl = dict()
                tbl['per_epoch_mean'] = np.zeros_like(scifiles, dtype = float)
                tbl['per_epoch_err'] = np.zeros_like(scifiles, dtype = float)
                tbl['per_epoch_DDV'] = np.zeros_like(scifiles, dtype = float)
                tbl['per_epoch_DDVRMS'] = np.zeros_like(scifiles, dtype = float)
                tbl['per_epoch_DDDV'] = np.zeros_like(scifiles, dtype = float)
                tbl['per_epoch_DDDVRMS'] = np.zeros_like(scifiles, dtype = float)
                tbl['LOCAL_FILE_NAME'] = scifiles
                tbl['plot_date'] = np.zeros_like(scifiles, dtype = float) # time for matplotlib

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
            # for date plotting
            tbl['plot_date'][i] = Time(hdr['MJDATE'], format='mjd').plot_date - Time(40588, format='mjd').plot_date


            # read line file
            tbl_per_line_ini = fits.getdata(scifiles[i])

            # get lines in the proper domain
            if i == 0:
                g = (tbl_per_line_ini['WAVE_START'] > 900) * \
                    (tbl_per_line_ini['WAVE_START'] < 2500) * \
                    (tbl_per_line_ini['NPIXLINE'] < 50) * \
                    (tbl_per_line_ini['RMSRATIO'] < 5)
                rvs = np.zeros([len(scifiles),np.sum(g)])
                dvrms = np.zeros([len(scifiles), np.sum(g)])
                DDV = np.zeros([len(scifiles),np.sum(g)])
                DDVRMS = np.zeros([len(scifiles), np.sum(g)])
                DDDV = np.zeros([len(scifiles),np.sum(g)])
                DDDVRMS = np.zeros([len(scifiles), np.sum(g)])

            # keep track in two big matrices of the rv and corresponding errors
            tbl_per_line = tbl_per_line_ini[g]

            rv_per_file = np.array( tbl_per_line['RV'])
            dvrms_per_file = np.array( tbl_per_line['DVRMS'])

            rvs[i] = rv_per_file
            dvrms[i] = dvrms_per_file


            if get_cumul_plot:
                if i ==0:
                    med_velo = np.nanmedian(rv_per_file)
                    xlim = [med_velo-5000,med_velo+5000]
                    best = dvrms_per_file<np.nanpercentile(dvrms_per_file,5)

                vrange = np.arange(xlim[0],xlim[1], 50.0)
                pdf = np.zeros_like(vrange,dtype = float)

                dvrms_per_file = dvrms_per_file[best]
                rv_per_file = rv_per_file[best]

                for ii in (range(len(rv_per_file))):
                    if np.isfinite(rv_per_file[ii])*np.isfinite(dvrms_per_file[ii]):
                        ww = et.exp(-0.5*(vrange-rv_per_file[ii])**2/dvrms_per_file[ii]**2)
                        ww /= dvrms_per_file[ii]
                        pdf += ww
                pdf/=np.nansum(pdf)
                ax[0].plot(vrange/1000,pdf,alpha = 0.2,color = 'grey')

                #cen, ew, amp, zp, slope
                p0 = [med_velo,500,np.max(pdf),0,0]
                #print(p0)
                fit  = et.fit_gauss(vrange,pdf, p0)

                ax[1].plot(vrange/1000,pdf - et.gauss(vrange, *fit),alpha = 0.2,color = 'grey')

            DDV[i] = tbl_per_line['DDV']
            DDVRMS[i] = tbl_per_line['DDVRMS']
            DDDV[i] = tbl_per_line['DDDV']
            DDDVRMS[i] = tbl_per_line['DDDVRMS']

        if get_cumul_plot:
            ax[0].set(xlabel = 'Velocity [km/s]',ylabel = 'Distribution function of RVs')
            ax[1].set(xlabel = 'Velocity [km/s]',ylabel = 'Distribution function of RVs - gaussfit')
            plt.savefig(outname+'_cumul.pdf')
            plt.close()

        # a line must be present >80% of the line
        valid_lines = np.nanmean(np.isfinite(rvs) * np.isfinite(dvrms),axis=0) > fcut
        rvs = rvs[:,valid_lines]
        dvrms = dvrms[:,valid_lines]
        tbl_per_line = tbl_per_line[valid_lines]
        DDV  = DDV[:,valid_lines]
        DDDV = DDDV[:,valid_lines]
        DDVRMS  = DDVRMS[:,valid_lines]
        DDDVRMS = DDDVRMS[:,valid_lines]

        sig =  np.nanpercentile((rvs- np.nanmedian(rvs))/dvrms,[16,84],axis=0)
        sig = (sig[1]-sig[0])/2


        valid_lines = sig<1.2
        rvs = rvs[:,valid_lines]
        dvrms = dvrms[:,valid_lines]
        tbl_per_line = tbl_per_line[valid_lines]
        DDV  = DDV[:,valid_lines]
        DDDV = DDDV[:,valid_lines]
        DDVRMS  = DDVRMS[:,valid_lines]
        DDDVRMS = DDDVRMS[:,valid_lines]



        if common_weights:
            for i in range(rvs.shape[1]):
                dvrms[:, i] *= sig[i]

            # produce a map of median change in error, liked to SNR changes
            err_ratio = np.zeros(dvrms.shape)
            # first guess at median per-line error
            ref = np.nanmedian(dvrms, axis=0)
            for i in range(dvrms.shape[0]):
                err_ratio[i] = et.nanmedian(dvrms[i]/ref)
            # better estimate of median rms
            ref = np.nanmedian(dvrms / err_ratio, axis=0)
            for i in range(dvrms.shape[0]):
                amp = et.nanmedian(dvrms[i]/ref)
                # all lines have the same relative weights
                # but may vary from one spectra to the other while
                # preserving the same relative ratios
                dvrms[i] = amp*ref

        # constructing a per-epoch mean velocity
        for i in tqdm(range(len(scifiles))):
            guess,bulk_error  = et.odd_ratio_mean(rvs[i],dvrms[i])
            tbl['per_epoch_mean'][i] = guess
            tbl['per_epoch_err'][i] = bulk_error

            guess,bulk_error  = et.odd_ratio_mean(DDV[i],DDVRMS[i])

            tbl['per_epoch_DDV'][i] = guess
            tbl['per_epoch_DDVRMS'][i] = bulk_error

            guess,bulk_error  = et.odd_ratio_mean(DDDV[i],DDDVRMS[i])

            tbl['per_epoch_DDDV'][i] = guess
            tbl['per_epoch_DDDVRMS'][i] = bulk_error

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
        bands =   ['Y',  'gapYJ',  'J',    'gapJH',  'H',    'gapHK',   'K',   'redK']
        blue_end =[ 900.0,1113.400,1153.586,1354.422,1462.897,1808.544,1957.792,2343.105]
        red_end = [1113.4,1153.586,1354.422,1462.897,1808.544,1957.792,2343.105,2500.000]
        suffix = ['', '_0-2044', '_2044-4088',  '_1532-2556']


        rvs_matrix = np.zeros([len(scifiles),len(bands),len(suffix)])+np.nan
        err_matrix = np.zeros([len(scifiles),len(bands),len(suffix)])+np.nan


        keep = np.ones_like(scifiles,dtype = bool)
        # updating the table with per-band
        for i in tqdm(range(len(scifiles))):
            tmp_rv = rvs[i] - rv_per_line_model[i]
            tmp_err = dvrms[i]

            for iband in range(len(bands)):
                for reg in range(4):
                    if reg == 0:
                        g = (tbl_per_line['WAVE_START'] > blue_end[iband]) * (tbl_per_line['WAVE_START'] < red_end[iband])
                    if reg == 1:
                        g = (tbl_per_line['WAVE_START'] > blue_end[iband]) * (tbl_per_line['WAVE_START'] < red_end[iband]) * (tbl_per_line['XPIX'] < 2044)
                    if reg == 2:
                        g = (tbl_per_line['WAVE_START'] > blue_end[iband]) * (tbl_per_line['WAVE_START'] < red_end[iband]) * (tbl_per_line['XPIX'] > 2044)
                    if reg == 3:
                        g = (tbl_per_line['WAVE_START'] > blue_end[iband]) * (tbl_per_line['WAVE_START'] < red_end[iband])  * (tbl_per_line['XPIX'] > 1532)  * (tbl_per_line['XPIX'] < 2556)


                    if np.sum(np.isfinite(tmp_err[g]) * np.isfinite(tmp_rv[g]))<np.sum(g)/2:
                        print(et.color('Less than 50% of lines are valid for {0}, band {1}, reg {2}'.format(scifiles[i], bands[iband],reg),'red'))
                        keep[i] = False
                        continue

                    if np.sum(g) < 5:
                        continue
                    guess,bulk_error  = et.odd_ratio_mean(tmp_rv[g],tmp_err[g])

                    rvs_matrix[i,iband,reg] = guess
                    err_matrix[i,iband,reg] = bulk_error

        for iband in range(len(bands)):
            for reg in range(4):
                tbl['per_epoch_mean_' + bands[iband]+suffix[reg]] = rvs_matrix[:,iband,reg]
                tbl['per_epoch_err_' + bands[iband]+suffix[reg]] = err_matrix[:,iband,reg]

        tbl = et.td_convert(tbl)
        tbl = tbl[keep]
        tbl.write(outname, overwrite = True)

        # plot or not
        if doplot:
            fig, ax = plt.subplots(nrows = 2, ncols = 1,sharex = True)
            #ax[0].errorbar(tbl['MJDATE'], tbl['per_epoch_mean_J']-np.nanmedian(tbl['per_epoch_mean_J']) , fmt='.g', yerr=tbl['per_epoch_err_J'],alpha = 0.3,label = 'J')
            ax[0].errorbar(tbl['MJDATE'], tbl['per_epoch_mean_H']-et.nanmedian(tbl['per_epoch_mean_H']) , fmt='.r', yerr=tbl['per_epoch_err_H'],alpha = 0.2,label = 'H')
            ax[0].errorbar(tbl['MJDATE'], tbl['per_epoch_mean']-et.nanmedian(tbl['per_epoch_mean']) , fmt='.k', yerr=tbl['per_epoch_err'],alpha = 0.8,label = 'all')
            #ax[2].errorbar(tbl['MJDATE'],tbl['per_epoch_DDV'] ,fmt='.k', yerr=tbl['per_epoch_DDVRMS'], alpha = 0.7)

            diff_JH =  tbl['per_epoch_mean_J']-tbl['per_epoch_mean_H']
            ax[1].errorbar(tbl['MJDATE'],diff_JH - et.nanmedian(diff_JH) , fmt='.g', yerr=np.sqrt(tbl['per_epoch_err_J']**2+tbl['per_epoch_err_H']**2),alpha = 0.5,label = 'J')
            ax[0].legend()
            ax[0].set(xlabel = 'MJDATE', ylabel = 'RV [m/s]',title = 'H velocity')
            ax[1].set(xlabel = 'MJDATE', ylabel = 'RV [m/s]', title = 'J-H velo diff')
            #ax[2].set(xlabel = 'MJDATE', ylabel = '2nd deriv')
            plt.show()
    else:
        print('File {0} exists, we read it'.format(outname))
        print('You may use force=True to overwrite it')

    return Table.read(outname)

