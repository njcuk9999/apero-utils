import numpy as np
import glob
from astropy.io import fits
import os
from astropy.table import Table
from tqdm import tqdm
from scipy import constants
from datetime import datetime
import etienne_tools as et
import warnings
from time import time
from scipy import stats
from scipy.interpolate import InterpolatedUnivariateSpline as ius

def lbl(obj_sci,obj_template = None,doplot_ccf = False,doplot_debug = False, force = False,
          lblrv_path = 'lblrv/',mask_path = 'masks/',template_path = 'templates/',
          science_path = 'tellurics/',ref_blaze_file = '2498F798T802f_pp_blaze_AB.fits',
          noise_model = False,check_fp = False, science_search_string = '*e2dsff*.fits',
          template_file = None):
    """
    if True:
        obj_sci = 'GL436'
        force = True
        obj_template = None
        doplot_ccf = False
        doplot_debug = False
        lblrv_path = 'lblrv/'
        mask_path = 'masks/'
        template_path = 'templates/'
        science_path = 'tellurics/'
        ref_blaze_file = '2498F798T802f_pp_blaze_AB.fits'
        check_fp = False
        noise_model = False
    """
    # pass just one object and we assume that the object is it's own template
    if obj_template is None:
        obj_template = obj_sci

    if doplot_ccf or doplot_debug:
        import matplotlib.pyplot as plt

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Inputs to the code
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    #obj_sci,obj_template = 'TOI-1452', 'GL699'
    #doplot_ccf,doplot_debug,force  = False, False, False

    #lblrv_path = 'lblrv/'  # must exist, path to catalog of lines
    #mask_path = 'masks/' # path to mask files, must have _pos files in there
    #template_path = 'templates/' # stellar templates
    #science_path = 'tellurics/' # direcory with per-object telluric-corect (tcorr) files are hidden


    if template_file == None:
        template_file = template_path+'Template_s1d_'+obj_template+'_sc1d_v_file_AB.fits'

    template = fits.getdata(template_file)

    scifiles = glob.glob(science_path+obj_sci+'/'+science_search_string)


    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    End of inputs to the code
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Start of checks
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # reference table that will always be used afterward as a starting point
    ref_table_name = 'ref_table_{0}.csv'.format(obj_template)
    maskfile = mask_path+obj_template+'_pos.fits'


    if not os.path.isfile(ref_table_name):
        print('File {0} will be created, that''s fine but you have been warned!'.format(ref_table_name))
    else:
        print('File {0} exists, it will be read, if you want to generate it again, delete the file'.format(ref_table_name))

    if not os.path.isdir(lblrv_path):
        raise ValueError('Directory {0} does not exist, that''s really bad, we exit'.format(lblrv_path))
    else:
        print(et.color('Direcotry {0} exists, we are happy'.format(lblrv_path),'green'))

    if not os.path.isdir(mask_path):
        raise ValueError(et.color('Directory {0} does not exist, that''s really bad, we exit'.format(mask_path),'red'))
    else:
        print(et.color('Direcotry {0} exists, we are happy'.format(mask_path),'green'))

    if not os.path.isdir(template_path):
        raise ValueError('Directory {0} does not exist, that''s really bad, we exit'.format(template_path))
    else:
        print(et.color('Direcotry {0} exists, we are happy'.format(template_path),'green'))

    if not os.path.isdir(science_path):
        raise ValueError('Directory {0} does not exist, that''s really bad, we exit'.format(science_path))
    else:
        print(et.color('Direcotry {0} exists, we are happy'.format(science_path),'green'))

    if not os.path.isdir(science_path+'/'+obj_sci):
        raise ValueError('Directory {0} does not exist, that''s really bad, we exit'.format(science_path+'/'+obj_sci))
    else:
        print(et.color('Direcotry {0} exists, we are happy'.format(science_path+'/'+obj_sci),'green'))

    if not os.path.isfile(maskfile):
        raise ValueError('File {0} does not exist, that''s really bad, we exit'.format(maskfile))
    else:
        print(et.color('File {0} exists, we are happy!'.format(maskfile),'green'))

    if not os.path.isfile(ref_blaze_file):
        raise ValueError('File {0} does not exist, that''s really bad, we exit'.format(ref_blaze_file))
    else:
        print(et.color('File {0} exists, we are happy!'.format(ref_blaze_file),'green'))

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    End of checks
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    print(et.color('We get bad odometers','magenta'))
    bad_odo = et.get_bad_odo()

    # sort files by name so that they are consecutive in time
    scifiles.sort()

    all_durations = []

    mask_hdr = fits.getheader(maskfile)

    if os.path.isfile(ref_table_name) == False:
        # we create an empty table that serves as the baseline for all future observations
        # using that template as a starting point. This table is read if it exists.
        wavegrid = et.fits2wave(scifiles[0])
        tbl = et.td_convert(Table.read(maskfile))


        wave_start = []
        wave_end = []
        order = []
        xpix = []
        weight_line = [] # used only to get systemic velocity as a starting point
        for i in range(wavegrid.shape[0]):
            gg = (tbl['ll_mask_s']>np.min(wavegrid[i]))*(tbl['ll_mask_s']<np.max(wavegrid[i]))
            if np.sum(gg) !=0:
                order = np.append(order,np.ones(np.sum(gg)-1)*i)
                wave_start = np.append(wave_start, np.array(tbl['ll_mask_s'][gg][:-1]))
                wave_end = np.append(wave_end, np.array(tbl['ll_mask_s'][gg][1:]))
                weight_line = np.append(weight_line,np.array(tbl['w_mask'][gg][:-1]) )
                xpix = np.append(xpix, ius(wavegrid[i],np.arange(len(wavegrid[i])))(np.array(tbl['ll_mask_s'][gg][:-1])))


        tbl = dict()
        tbl['ORDER'] = np.array(order,dtype = int)
        tbl['WAVE_START'] = wave_start #np.array(wave_start*(1 + objrv * 1000 / constants.c))
        tbl['WAVE_END'] = wave_end#np.array(wave_end* (1 + objrv * 1000 / constants.c))
        tbl['WEIGHT_LINE'] = weight_line # from mask
        tbl['XPIX'] = xpix #pixel position along array
        tbl['RMSRATIO'] = np.zeros_like(xpix) # ratio of expected VS actual RMS in difference of model vs line
        tbl['NPIXLINE'] = np.zeros_like(xpix,dtype = int) # effective number of pixels in line
        tbl['MEANXPIX'] = np.zeros_like(xpix,dtype = float) # mean line position in pixel space
        tbl['MEANBLAZE'] = np.zeros_like(xpix,dtype = float) # blaze value compared to peak for that order
        tbl['AMP_CONTINUUM'] = np.zeros_like(xpix,dtype = float) # blaze value compared to peak for that order

        # Considering the number of pixels, expected and actual RMS, this is the likelihood that the line is
        # acually valid from a Chi2 test point of view
        tbl['CHI2'] = np.zeros_like(xpix)
        # probability of valid considering the chi2 CDF for the number of DOF
        tbl['CHI2_VALID_CDF'] = np.zeros_like(xpix)

        print('We write {0}'.format(ref_table_name))
        et.td_convert(tbl).write(ref_table_name)
    else:
        print('We read {0}'.format(ref_table_name))
        tbl = et.td_convert(Table.read(ref_table_name))
        order = np.array(tbl['ORDER'] )
        wave_start = np.array(tbl['WAVE_START'])
        wave_end = np.array(tbl['WAVE_END'])
        weight_line = np.array(tbl['WEIGHT_LINE'] )

    # we create the spline of the template to be used everywhere further down
    valid = np.isfinite(template['flux'])
    bl = fits.getdata(ref_blaze_file)
    for i in range(bl.shape[0]):
        bl[i]/=np.nanpercentile(bl[i],90)

    # get info on template systvel for splining correctly
    systemic_vel = -1000*mask_hdr['SYSTVEL']

    print('defining all the splines required later')
    flux = np.array(template['flux'])

    flux0 =np.array(flux)


    # we select a scale of 223 km/s
    width = int(223000/np.array(1/np.nanmedian((template['wavelength']/np.gradient(template['wavelength']))/constants.c)))
    if (width % 2) ==0:
        width+=1

    # we high-pass on a scale of ~101 pixels in the e2ds
    flux -= et.lowpassfilter(flux, width=width)

    dflux = np.gradient(flux)/ np.gradient(np.log(template['wavelength'])) / constants.c
    # 2nd derivative
    ddflux = np.gradient(dflux) / np.gradient(np.log(template['wavelength'])) / constants.c
    # 3rd derivative
    dddflux = np.gradient(ddflux) / np.gradient(np.log(template['wavelength'])) / constants.c

    # we create the spline of the template to be used everywhere further down
    valid = np.isfinite(flux) * np.isfinite(dflux) * np.isfinite(ddflux)* np.isfinite(dddflux)

    # template removed from its systemic velocity, the spline is redefined for good
    #
    # Flux, 1st and 2nd derivatives. We could add more derivatives if need be one day
    # Be careful to update the 'valid' mask above

    k_order = 5
    spline0 = ius(et.doppler(template['wavelength'][valid],-systemic_vel),flux0[valid],k=k_order,ext=1)
    spline = ius(et.doppler(template['wavelength'][valid],-systemic_vel),flux[valid],k=k_order,ext=1)
    dspline = ius(et.doppler(template['wavelength'][valid],-systemic_vel),dflux[valid],k=k_order,ext=1)

    ddspline = ius(et.doppler(template['wavelength'][valid],-systemic_vel),ddflux[valid],k=k_order,ext=1)
    dddspline = ius(et.doppler(template['wavelength'][valid],-systemic_vel),dddflux[valid],k=k_order,ext=1)

    # we create a mask to know if the splined point  is valid
    spline_mask = ius(et.doppler(template['wavelength'],-systemic_vel),np.isfinite(template['flux']),k=1,ext=1)

    # all systemic velocities. Used to guesstimate the velocity
    systemic_all = np.zeros_like(scifiles, dtype = float) + np.nan
    mjdate_all = np.zeros_like(scifiles, dtype = float)

    CCF_EWIDTH = None

    # flag to take a completely new RV measurement
    failed_convergence = True
    for ifile in range(len(scifiles)):

        # output name
        outname = lblrv_path + scifiles[ifile].split('.fits')[0].split('/')[-1]+'_'+obj_sci+'_'+obj_template+'_lbl.fits'

        if os.path.isfile(outname) and not force:
            print(et.color('\t\tfile {0} exists and force = {1}'.format(outname, force),'cyan'))
            continue

        # get the science file info
        sp, hdr = fits.getdata(scifiles[ifile], header=True)
        if hdr['INSTRUME'] == 'HARPS':
            print(et.color('This is a HARPS file','cyan'))
            hdr = et.harps2spirou(hdr)

        if 'EXPNUM' in hdr.keys():
            if str(hdr['EXPNUM']) in bad_odo:
                print(et.color('File {} is a known bad file, we skip'.format(scifiles[ifile]),'red'))
                continue
            else:
                print(et.color('File {} is *not* a known bad file'.format(scifiles[ifile]),'cyan'))

        if 'EXTSN035' in hdr.keys():
            if hdr['EXTSN035']<10:
                print(et.color('We have an SNR of <10 in H band, we skip file','red'))
                continue
            else:
                print(et.color('SNR in H is {0:.2f} and >10, all good'.format(hdr['EXTSN035']),'cyan'))

        time_start = time()

        if noise_model:
            rms = get_noise_model(scifiles[ifile])
        else:
            rms = np.zeros_like(sp)

        sp0 = np.array(sp)
        # get the wave grid from science file

        wave = et.fits2wave(hdr)
        for iord in range(sp.shape[0]):

            # we select a scale of 223 km/s
            width = int(
                223000 / np.array(1 / np.nanmedian((wave[iord] / np.gradient(wave[iord])) / constants.c)))
            if (width % 2) == 0:
                width += 1

            sp[iord] -= et.lowpassfilter(sp[iord], width = width)


        if 'FP' not in obj_sci:
            BERV = hdr['BERV']*1000
        else:
            BERV = 0

        if failed_convergence:
            if ('FP' not in obj_sci):
                rv,ewidth = et.get_rough_ccf_rv(wave, sp, wave_start, np.ones_like(weight_line), doplot=doplot_ccf)
                if CCF_EWIDTH == None:
                    CCF_EWIDTH = np.array(ewidth)
                    print(et.color('CCF e-wdith = {:.2f} m/s'.format(CCF_EWIDTH),'blue'))
            else:
                rv = 0
                CCF_EWIDTH = 0
        else:
            rv = systemic_all[np.argmin(hdr['MJDATE'] - mjdate_all)] + BERV# + 1000*np.random.random()

        keep = np.ones_like(tbl['ORDER'],dtype = bool)
        done = False
        ite_convergence = 10

        for ite_rv in range(1,10):
            if done:
                continue

            now = datetime.now()

            if (ite_rv == 1):
                # spline the model onto the wavelenght grid at the specified velocity
                model = np.zeros_like(sp)
                model0 = np.zeros_like(sp)
                dmodel = np.zeros_like(sp)
                ddmodel = np.zeros_like(sp)
                dddmodel = np.zeros_like(sp)
                model_mask = np.zeros_like(sp)

            for ii in range(sp.shape[0]):
                # RV shift the spline and give it the shape of the model
                model[ii] = spline(et.doppler(wave[ii],-rv) ) * bl[ii]

                amp =  np.nansum(model[ii]*sp[ii])/np.nansum(model[ii]**2)
                model[ii]*=amp


                if ite_rv == 1:
                    model0[ii] = spline0(et.doppler(wave[ii],-rv) ) * bl[ii]
                    model0[ii]/=np.nanmedian(model0[ii])
                    model0[ii]*=np.nanmedian(sp0[ii])

                dmodel[ii] = dspline(et.doppler(wave[ii],-rv)) * bl[ii] * amp

                ddmodel[ii] = ddspline(et.doppler(wave[ii],-rv)) * bl[ii]* amp
                dddmodel[ii] = dddspline(et.doppler(wave[ii],-rv)) * bl[ii]* amp

                # if on first iteration, get the low-frequency component out
                model_mask[ii] = spline_mask(et.doppler(wave[ii], -rv))
                g = tbl['ORDER'] == ii

            # keep only legit splined points
            model[model_mask < .99] = np.nan

            for ii in range(sp.shape[0]):
                tmp = sp[ii]-model[ii]

                if not noise_model:
                    ipix = np.arange(0,model.shape[1],100)
                    sig = np.zeros_like(ipix,dtype = float)
                    for ipixi in range(len(ipix)):
                        i1 = ipix[ipixi] - 100
                        i2 = ipix[ipixi] + 100
                        if i1 <0:
                            i1 = 0
                        if i2>model.shape[1]:
                            i2 = model.shape[1]

                        sig[ipixi] = et.sigma(tmp[i1:i2])

                        sig[sig == 0] = np.nan
                        gg = np.isfinite(sig)
                    if np.sum(gg)>2:
                        # we must have at least two bins to do anything useful here
                        rms[ii] = ius(ipix[gg], sig[gg], k=1, ext=3)(np.arange(model.shape[1]))
                    else:
                        rms[ii] = np.nan

            if ite_rv == 1:
                # create a dummy array that will contain velocities and corresponding errors
                dv = np.zeros(len(tbl['WAVE_START']))+np.nan
                ddv = np.zeros(len(tbl['WAVE_START']))+np.nan
                dddv = np.zeros(len(tbl['WAVE_START']))+np.nan

                dvrms = np.zeros(len(tbl['WAVE_START']))+np.nan
                ddvrms = np.zeros(len(tbl['WAVE_START']))+np.nan
                dddvrms = np.zeros(len(tbl['WAVE_START']))+np.nan

            # keep track of which order we are looking at
            iord_prev = -1

            # order of each line in the table
            iord = np.array(tbl['ORDER'])

            # noinspection PyInterpreter
            # loop through all lines
            for i in range(0,len(tbl['ORDER'])):#, leave = False):
                # if the line has been flagged as bad, skip right away
                if (ite_rv !=1) and (keep[i] == False):
                    continue

                # if not the same as previous order, then get residuals for that order
                if iord[i] != iord_prev:
                    iord_prev = iord[i]

                    ww_ord = et.doppler(wave[iord[i]],-rv)
                    sp_ord = sp[iord[i]]

                    wave2pix = ius(ww_ord,np.arange(model.shape[1]))

                    rms_ord = rms[iord[i]]
                    model_ord = model[iord[i]]

                    dmodel_ord = dmodel[iord[i]]
                    ddmodel_ord = ddmodel[iord[i]]
                    dddmodel_ord = dddmodel[iord[i]]

                    bl_ord = bl[iord[i]]

                    # we will perform a line-by-line amplitude correction
                    #
                    if doplot_debug:
                        if hdr['INSTRUME'] == 'HARPS':
                            if (iord[i] == 60)*(ite_rv == 1):
                                plt.plot(ww_ord,model_ord, color = 'grey', linewidth=3, alpha = 0.3,label = 'template')
                        else:
                            if (iord[i] == 35)*(ite_rv == 1):
                                plt.plot(ww_ord,model_ord, color = 'grey', linewidth=3, alpha = 0.3,label = 'template')


                imin, imax = wave2pix([wave_start[i], wave_end[i]])
                imin, imax = int(np.floor(imin)),int(np.ceil(imax))

                if (imax-imin)<5:
                    keep[i] == False
                    continue
                if imin<0:
                    continue
                if imax>( ww_ord.shape[0]-2):
                    continue

                # get weights at the edge of the domain. Pixels inside have a weight of 1, at the edge, it's proportional
                # to the overlap
                weight_mask = np.ones(imax-imin+1)

                if ww_ord[imin]<wave_start[i]:
                    weight_mask[0] = (ww_ord[imin+1]-wave_start[i])/(ww_ord[imin+1]-ww_ord[imin])
                if ww_ord[imax+1] > wave_end[i]:
                    weight_mask[-1] = 1-( ww_ord[imax] - wave_end[i] )/(ww_ord[imax+1]-ww_ord[imax])

                xpix = np.arange(len(weight_mask))+imin
                tbl['MEANXPIX'][i] = et.nansum(weight_mask*xpix)/np.nansum(weight_mask)
                tbl['MEANBLAZE'][i] = bl_ord[(imin+imax)//2]
                # maybe some plots
                if doplot_debug:
                    if hdr['INSTRUME'] == 'HARPS':
                        if (iord[i] == 60)*(ite_rv == 1):
                            color = (['red','green','blue'])[i % 3]
                            plt.plot(ww_ord[imin:imax+1],sp_ord[imin:imax+1],color = color)
                    else:
                        if (iord[i] == 35)*(ite_rv == 1):
                            color = (['red','green','blue'])[i % 3]
                            plt.plot(ww_ord[imin:imax+1],sp_ord[imin:imax+1],color = color)







                # derivative of the segment
                d_segment = dmodel_ord[imin:imax+1]*weight_mask

                # keep track of 2nd derivative
                dd_segment = ddmodel_ord[imin:imax + 1] * weight_mask
                ddd_segment = dddmodel_ord[imin:imax + 1] * weight_mask

                # residuals of the segment

                sp_segment = (sp_ord[imin:imax+1])
                model_segment = (model_ord[imin:imax+1])


                diff_segment = (sp_segment - model_segment)*weight_mask

                sum_weight_mask = et.sum(weight_mask)  # we need this value 4 times
                mean_rms = et.sum(rms_ord[imin:imax+1]*weight_mask)/sum_weight_mask

                # 1st derivative
                # From Bouchy 2001 equation, RV error for each pixel
                dvrms_pix = mean_rms/ d_segment
                # RV error for the line
                dvrms[i] = 1/np.sqrt(np.sum((1/dvrms_pix**2)))
                # feed the monster
                dv[i] = et.sum(diff_segment*d_segment)/et.sum(d_segment**2)

                # 2nd derivative
                # From Bouchy 2001 equation, RV error for each pixel
                with warnings.catch_warnings(record=True) as _:
                    ddvrms_pix = mean_rms / dd_segment
                # RV error for the line
                ddvrms[i] = 1 / np.sqrt(np.sum((1 / ddvrms_pix ** 2)))
                ddv[i] = et.sum(diff_segment * dd_segment) / np.sum(dd_segment ** 2)


                # 3rd derivative
                # From Bouchy 2001 equation, RV error for each pixel
                with warnings.catch_warnings(record=True) as _:
                    dddvrms_pix = mean_rms / ddd_segment
                    # RV error for the line
                    dddvrms[i] = 1 / np.sqrt(np.sum((1 / dddvrms_pix ** 2)))
                    dddv[i] = et.sum(diff_segment * ddd_segment) / np.sum(ddd_segment ** 2)

                #ratio of expected VS actual RMS in difference of model vs line
                tbl['RMSRATIO'][i] = et.nanstd(diff_segment)/mean_rms
                # effective number of pixels in line
                tbl['NPIXLINE'][i] = len(diff_segment)
                # Considering the number of pixels, expected and actual RMS, this is the likelihood that the line is
                # acually valid from a chi2 point of view
                tbl['CHI2'][i] = et.nansum( (diff_segment/mean_rms)**2)


            if doplot_debug:
                if ite_rv ==1:
                    plt.xlabel('Wavelength [nm]')
                    plt.ylabel('Arbitrary flux')
                    plt.show()

            nsig = dv / dvrms
            nsig = nsig[np.isfinite(nsig)]
            nsig = nsig[np.abs(nsig) < 8]
            stddev_nsig = et.sigma(nsig)

            print('\t\tstdev_meas/stdev_pred = {0:.2f}'.format(stddev_nsig))
            # we force to one
            #dvrms *= stddev_nsig

            # get the best estimate of the velocity and update spline
            rv_mean, bulk_error = et.odd_ratio_mean(dv,dvrms)

            # some printout
            print('\t\tite_rv {0}, bulk error {1:.2f} m/s, rv = {2:.2f} m/s '.format(ite_rv,bulk_error, -rv))
            print('\t\tDuration of loop iteration : {0}'.format( datetime.now()-now))

            rv_finale = np.array(dv + rv - BERV)
            rv += rv_mean

            #if ite_rv == 1:
            #    rv1 = np.array(rv_finale)

            print(et.color('\tRV = {0:.2f} m/s, sigma = {1:.2f} m/s'.format(rv_mean,bulk_error),'green'))

            if np.abs(rv_mean) < bulk_error/5:
                ite_convergence = ite_rv
                done = True

        tbl['RV'] = -rv_finale # express to have sign fine relative to convention
        tbl['DVRMS'] = dvrms

        # adding to the fits table the 2nd derivative projection
        tbl['DDV'] = ddv
        tbl['DDVRMS'] = ddvrms
        # adding to the fits table the 3rd derivative projection

        tbl['DDDV'] = dddv
        tbl['DDDVRMS'] = dddvrms


        systemic_all[ifile] = rv-BERV

        mjdate_all[ifile] = hdr['MJDATE']

        if ite_convergence >=9:
            print(et.color('\t\tWe assume that the RV is (probably) bad. Next step we''ll measure it with a CCF','red'))
            failed_convergence = True
        else:
            print(et.color('\t\tConverged in {0} steps'.format(ite_convergence),'green'))
            failed_convergence = False

        print(et.color('\tWe write file {0} [{1}/{2}]'.format(outname,ifile, len(scifiles)),'magenta'))


        hdu1 = fits.PrimaryHDU()

        hdr['ITE_RV'] = ite_convergence,'N of iter to reach 0.1 sigma num accuracy'
        hdr['SYSTVELO'] = -systemic_vel, 'Systemic velocity in m/s'
        hdr['RMSRATIO'] = stddev_nsig, 'RMS vs photon noise'

        for key in hdr.keys():
            if 'ESO ' in key:
                del hdr[key]

        # be sure that the header is healthy
        for key in hdr.keys():
            try:
                val = hdr[key]
            except:
                print(et.color('key {} was bad, it has been padded'.format(key),'red'))
                hdr[key] = ''

        conversion_d2_fwhm = CCF_EWIDTH*1e3 # CCF_EWIDTH is expressed in m/s
        tbl['fwhm'] =  (tbl['DDV']/conversion_d2_fwhm+CCF_EWIDTH)* np.sqrt(np.log(2)*2)*2
        tbl['sig_fwhm'] = (tbl['DDVRMS']/conversion_d2_fwhm)* np.sqrt(np.log(2)*2)*2
        hdr['CCF_WIDTH'] = CCF_EWIDTH,'e-width of LBL CCF in m/s'


        hdu1.header = hdr
        # to handle problematic keys in HARPS headers
        hdu1.verify('fix')

        tbl['CHI2_VALID_CDF'] = 1-stats.chi2.cdf(tbl['CHI2'],tbl['NPIXLINE'])



        # convert back from dictionnary to table and save
        hdu2 = fits.BinTableHDU(et.td_convert(tbl))
        new_hdul = fits.HDUList([hdu1, hdu2])
        new_hdul.writeto(outname, overwrite=True)

        all_durations = np.append(all_durations,time()-time_start)

        nleft = len(scifiles) - ifile - 1

        if len(all_durations) >2:
            print(et.color('\tDuration per file {0:.2f}+-{1:.2f}s'.format(et.nanmean(all_durations),
                                                                          et.nanstd(all_durations)), 'yellow'))
            print(et.color('\tTime left to completion {0}, {1} / {2} files todo/done\n'.format(et.smart_time(et.nanmean(all_durations)*nleft), nleft,ifile), 'white'))
        else:
            print()