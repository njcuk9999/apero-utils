import numpy as np
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
from fits2wave import fits2wave
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.signal import medfilt
from tqdm import tqdm
import os


def sanitize(files, blaze_file, model_s1d_file, prob_bad = 1e-3, Force = False, doplot = False):
    #
    # Program to create *sanitized* files, where outlying points are replaced with
    # the model
    #
    # Sample inputs --
    #
    # Files to be sanitized, needs to be torr
    #    files = glob.glob('Gl699/2*o_pp_e2dsff_tcorr_AB.fits')
    #
    # Blaze file
    #    blaze_file = '2019-04-20_2400404f_pp_blaze_AB.fits'
    #
    # s1d_v model file for the object
    #    model_s1d_file = 'Template_s1d_Gl699_sc1d_v_file_AB.fits'
    #
    # sanitize(files, blaze_file,model_s1d_file)
    #
    # optional inputs --
    #
    #   prob_bad
    #     probability that a point in a spectrum is bad. Leave at 1e-3
    #     if you don't have a good reason to modify
    #
    #   Force
    #     force rewriting of files. Useful if you modify the code
    #   doplot
    #     show debug plots
    #


    # read the object's model
    tbl = Table(fits.getdata(model_s1d_file))
    # read the blaze file
    blaze = fits.getdata(blaze_file)

    # normalize the blaze
    for i in range(49):
        blaze[i] /= np.nanmedian(blaze[i])

    # get the template to be used for comparison
    wave_template = tbl['wavelength']
    flux_template = tbl['flux']

    # keep valid pixels
    g = np.isfinite(flux_template)
    # spline valid template pixels
    template = InterpolatedUnivariateSpline(wave_template[g], flux_template[g], k=1, ext=1)

    # keep track of amplitudes per order
    amps = np.zeros(49)
    fgood = np.zeros(49)

    for i in range(len(files)):
        # loop through files
        print('Sanitizing file #{0} in {1} -> {2}'.format(i+1,len(files),files[i]))
        outname = 'sanitized'.join(files[i].split('tcorr'))

        if os.path.isfile(outname) and (Force == False):
            continue

        data, hdr = fits.getdata(files[i], header = True)

        nanmask = np.isfinite(data) == False

        dv = np.sqrt((1+hdr['BERV']/299792.458) /(1-hdr['BERV']/299792.458) )

        wave = fits2wave(hdr)*dv

        for iord in range(49):
            if np.max(np.isfinite(data[iord])) == False:
                data[iord] = np.nan
                continue


            if amps[iord] == 0:
                amps[iord] = np.nanmedian(data[iord])


            data_tmp =  data[iord]/np.nanmedian(data[iord])
            model_tmp = template(wave[iord])*blaze[iord]
            model_tmp[model_tmp == 0] = np.nan


            ker = np.exp(-.5*(np.arange(121)-60)**2/(20)**2)
            ker /=np.nansum(ker)
            for ite in range(1):

                ratio = np.nanmedian(data_tmp/model_tmp)#,)201)
                #ratio[]

                #mask = np.isfinite(ratio)
                #ratio[~mask] = 0
                #mask[~mask] = 0

                #w = np.convolve(mask*1.0,ker,mode = 'same')

                #ratio[w!=0] = np.convolve(ratio,ker,mode = 'same')[w!=0]/w[w!=0]
                #ratio[blaze[iord]<0.5] = 0

                #fig, ax = plt.subplots(nrows=3,ncols=1,sharex=  True)
                #ax[0].plot(data_tmp)
                #ax[0].plot(model_tmp)
                #ax[1].plot(ratio)

                data_tmp /=ratio#[ratio != 0]

                #plt.show()
            #stop

            residual = data_tmp - model_tmp

            sig = (np.nanpercentile(residual,84)-np.nanpercentile(residual,16))/2

            nsig = residual / sig

            prob_good = np.exp( -0.5*nsig**2)

            pp = prob_good / (prob_bad+prob_good)

            data[iord] = data_tmp*pp+model_tmp*(1-pp)

            data[iord][np.isnan(data_tmp)] = model_tmp[np.isnan(data_tmp)]

            min_valid = np.min( np.where(np.isfinite(data_tmp)) )
            max_valid = np.max( np.where(np.isfinite(data_tmp)) )

            tmp = data[iord]
            tmp[0:min_valid] = np.nan
            tmp[max_valid:] = np.nan
            data[iord] = tmp

            if doplot:
                plt.plot(wave[iord],data_tmp,'r', label = 'observed')
                plt.plot(wave[iord],model_tmp, label = 'model' )
                plt.plot(wave[iord],data[iord],'k', label = 'sanitized')
                plt.ylabel('normalized flux')
                plt.xlabel('wavelength (nm)')
                plt.legend()
                plt.show()

            data[iord]/=np.nanmedian(data[iord])
            data[iord]*=amps[iord]

            fgood[iord] = np.nanmean(pp)
            print('order #{0}, frac input {1:4.2f}%'.format(iord,fgood[iord]*100))

        print('We write file {0}'.format(outname))

        if doplot:
            plt.plot(1-fgood,'go')
            plt.xlabel('Nth order')
            plt.ylabel('Model padding fraction')
            plt.show()
        fits.writeto(outname, data, hdr, overwrite = True)

