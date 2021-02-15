import glob
from astropy.io import fits
import os
from tqdm import tqdm
import etienne_tools as et
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.time import Time

def compilbl(obj_sci, obj_template = None, doplot = False, force = True, common_weights = False,
               get_cumul_plot = False):

    """
    obj_sci = 'TRAPPIST-1'
    obj_template = None
    doplot = False
    force = True
    common_weights = False
    get_cumul_plot = False
    """

    if doplot or get_cumul_plot:
        # avoids conflicts if one uses ssh without graphic displays
        import matplotlib.pyplot as plt

    # pass just one object and we assume that the object is it's own template
    if obj_template is None:
        obj_template = obj_sci

    suffix = ''

    outname = 'lbl_{0}_{1}.rdb'.format(obj_sci, obj_template)
    outname2 = 'lbl2_{0}_{1}.rdb'.format(obj_sci, obj_template)

    # default keywords to be included in the compilation of per-object RVs
    keys = ['MJDATE', 'EXPTIME', 'AIRMASS', 'FILENAME',
            'DATE-OBS', 'BERV', 'TAU_H2O', 'TAU_OTHE',
            'ITE_RV', 'SYSTVELO',
            'TLPDVH2O','TLPDVOTR','CDBWAVE','OBJECT',
            'SBRHB1_P','SBRHB2_P','SBCDEN_P','SNRGOAL',
            'EXTSN035','BJD']
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
                tbl['rjd'] = np.zeros_like(scifiles, dtype = float)
                tbl['vrad'] = np.zeros_like(scifiles, dtype = float)
                tbl['svrad'] = np.zeros_like(scifiles, dtype = float)
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
            tbl['rjd'][i] = hdr['BJD'] - 2400000


            for key in keys:
                # adding header value in output table
                if key in hdr:
                    try:
                        tbl[key][i] = hdr[key]
                    except:
                        print(et.color('key {0} not present in file {1}'.format(key,scifiles[i]),'red'))
            # for date plotting
            tbl['plot_date'][i] = Time(hdr['MJDATE'], format='mjd').plot_date - Time(40588, format='mjd').plot_date

            # read line file
            tbl_per_line_ini = fits.getdata(scifiles[i])

            # get lines in the proper domain
            if i == 0:
                g = (tbl_per_line_ini['WAVE_START'] > 900) * \
                    (tbl_per_line_ini['WAVE_START'] < 2500) * \
                    (tbl_per_line_ini['NPIXLINE'] < 50) #* \
                    #(tbl_per_line_ini['RMSRATIO'] < 5)
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

        print('Forcing a stddev of 1 for all lines')
        print('constructing a per-epoch mean velocity')
        for i in tqdm(range(rvs.shape[0])):
            nsig = (rvs[i] - np.nanmedian(rvs[i])) / dvrms[i]
            dvrms[i]*=et.sigma(nsig)

            guess,bulk_error  = et.odd_ratio_mean(rvs[i],dvrms[i])
            tbl['vrad'][i] = guess
            tbl['svrad'][i] = bulk_error

        # de-biasing line. Matrix that contains a replicated 2d version of the per-epoch mean
        rv_per_epoch_model = np.reshape(np.repeat(tbl['vrad'],rvs.shape[1]),rvs.shape)

        # line-by-line mean position
        per_line_mean = np.zeros(rvs.shape[1])
        per_line_error = np.zeros(rvs.shape[1])

        # computing the per-line bias
        for i in tqdm(range(len(per_line_mean))):
            tmp1 = rvs[:,i]-rv_per_epoch_model[:,i]
            err1 = dvrms[:,i]
            try:
                guess,bulk_error  = et.odd_ratio_mean(tmp1,err1)
                per_line_mean[i] = guess
                per_line_error[i] = bulk_error
            except:
                per_line_mean[i] = np.nan
                per_line_error[i] = np.nan

        # we keep the per-line mean to zero
        per_line_mean -= (et.odd_ratio_mean(per_line_mean,per_line_error))[0]

        # constructing a 2d model of line biases
        rv_per_line_model = np.reshape(np.tile(per_line_mean,rvs.shape[0]),rvs.shape)

        # updating the table
        for i in tqdm(range(len(scifiles))):
            guess,bulk_error  = et.odd_ratio_mean(rvs[i]-rv_per_line_model[i],dvrms[i])
            tbl['vrad'][i] = guess
            tbl['svrad'][i] = bulk_error

            guess,bulk_error  = et.odd_ratio_mean(DDV[i],DDVRMS[i])
            tbl['per_epoch_DDV'][i] = guess
            tbl['per_epoch_DDVRMS'][i] = bulk_error

            guess,bulk_error  = et.odd_ratio_mean(DDDV[i],DDDVRMS[i])
            tbl['per_epoch_DDDV'][i] = guess
            tbl['per_epoch_DDDVRMS'][i] = bulk_error

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

                        continue

                    if np.sum(g) < 5:
                        continue

                    guess,bulk_error  = et.odd_ratio_mean(tmp_rv[g],tmp_err[g])

                    rvs_matrix[i,iband,reg] = guess
                    err_matrix[i,iband,reg] = bulk_error

        for iband in range(len(bands)):
            for reg in range(4):
                tbl['vrad_' + bands[iband]+suffix[reg]] = rvs_matrix[:,iband,reg]
                tbl['svrad_' + bands[iband]+suffix[reg]] = err_matrix[:,iband,reg]

        tbl = et.td_convert(tbl)
        tbl.write(outname, overwrite = True)

        # plot or not
        if doplot:
            fig, ax = plt.subplots(nrows = 2, ncols = 1,sharex = True)
            #ax[0].errorbar(tbl['MJDATE'], tbl['per_epoch_mean_J']-np.nanmedian(tbl['per_epoch_mean_J']) , fmt='.g', yerr=tbl['per_epoch_err_J'],alpha = 0.3,label = 'J')
            ax[0].errorbar(tbl['rjd'], tbl['vrad_H']-et.nanmedian(tbl['vrad_H']) , fmt='.r', yerr=tbl['svrad_H'],alpha = 0.2,label = 'H')
            ax[0].errorbar(tbl['rjd'], tbl['vrad']-et.nanmedian(tbl['svrad']) , fmt='.k', yerr=tbl['svrad'],alpha = 0.8,label = 'all')
            #ax[2].errorbar(tbl['MJDATE'],tbl['per_epoch_DDV'] ,fmt='.k', yerr=tbl['per_epoch_DDVRMS'], alpha = 0.7)

            diff_JH =  tbl['vrad_J']-tbl['vrad_H']
            ax[1].errorbar(tbl['rjd'],diff_JH - et.nanmedian(diff_JH) , fmt='.g', yerr=np.sqrt(tbl['svrad_J']**2+tbl['svrad_H']**2),alpha = 0.5,label = 'J')
            ax[0].legend()
            ax[0].set(xlabel = 'rjd', ylabel = 'RV [m/s]',title = 'H velocity')
            ax[1].set(xlabel = 'rjd', ylabel = 'RV [m/s]', title = 'J-H velo diff')
            #ax[2].set(xlabel = 'MJDATE', ylabel = '2nd deriv')
            plt.show()
    else:
        print('File {0} exists, we read it'.format(outname))
        print('You may use force=True to overwrite it')



    udates = np.unique( tbl['DATE-OBS'])

    tbl2 = Table(tbl[0:len(udates)]) # create a table with a per-epoch value

    for i in tqdm(range(len(udates))):
        tbl_date = tbl[udates[i] ==  tbl['DATE-OBS']]
        for key in tbl_date.keys():
            if 'vrad' not in key:
                try:
                    tbl2[key][i] = np.mean(tbl_date[key])
                except:
                    tbl2[key][i] = tbl_date[key][0]

            if key[0:4] == 'vrad':
                rv = tbl_date[key]

                err_rv = tbl_date['s'+key]
                tbl2[key][i] = np.nansum(rv/err_rv**2)/np.nansum(1/err_rv**2)

                tbl2['s'+key][i] = np.sqrt(1/np.nansum(1/err_rv**2))
    tbl2.write(outname2, overwrite = True)


    return Table.read(outname)
