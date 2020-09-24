import numpy as np
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
from fits2wave import fits2wave
from scipy.interpolate import InterpolatedUnivariateSpline
from tqdm import tqdm
import os
from tqdm import tqdm
import warnings
#from sigma import *
from scipy.signal import medfilt

def sigma(im):
    return (np.nanpercentile(im,0.841344746068543) - np.nanpercentile(im,0.15865525393145702))/2



def lowpassfilter(v,w):
    # low-pass filtering of vector v using a running median of width w. This code differs from medfilt in proprely
    # handling NaNs
    v3 = []
    i3 = []
    #
    for i in range(0,len(v),w//3):

        i0 = i-w//2
        i1 = i+w//2

        if i0 <0:
            i0 = 0
        if i1 > (len(v)):
            i1 = len(v)
        vv =v[i0:i1]

        if np.max(np.isfinite(vv)):
            v3.append(np.nanmedian(vv))
            i3.append(i)

    v3 = np.array(v3)
    i3 = np.array(i3)

    if len(i3)<5:
        return np.zeros_like(v)+np.nan


    spline = InterpolatedUnivariateSpline(i3,v3, k=2, ext=3)
    return spline(np.arange(len(v)))


def running_sigma(v, w):
    # running stddev witnin a box of width w
    if (w % 2) == 0:
        w+=1

    res = v - lowpassfilter(v, w)
    rms = lowpassfilter(np.abs(res), w)

    rms = rms*sigma(res)/np.nanmedian(np.abs(res))

    return rms


def sanitize(files, blaze_file, model_s1d_file, prob_bad = 1e-3, force = False, doplot = -1):
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
    #    blaze_file = 'varia/2019-04-20_2400404f_pp_blaze_AB.fits'
    #
    # s1d_v model file for the object
    #    model_s1d_file = 'varia/Template_s1d_Gl699_sc1d_v_file_AB.fits'
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
    flux_template[flux_template <0] = 0

    # keep valid pixels
    g = np.isfinite(flux_template)

    # spline valid template pixels
    template = InterpolatedUnivariateSpline(wave_template[g], flux_template[g], k=1, ext=1)
    mask_template = InterpolatedUnivariateSpline(wave_template, np.isfinite(flux_template), k=1, ext=1)

    # keep track of amplitudes per order
    amps = np.zeros(49)
    fgood = np.zeros(49)

    for i in (range(len(files))):
        plt.close()

        # loop through files
        print('Sanitizing file #{0} in {1} -> {2}'.format(i+1,len(files),files[i]))
        outname = 'sani'.join(files[i].split('tcorr'))

        if os.path.isfile(outname) and (force == False):
            print(outname+' exisits, we skip')
            continue

        data, hdr = fits.getdata(files[i], header = True)

        with warnings.catch_warnings(record=True) as _:
            data[blaze < 0.3] = np.nan


        prob_map = np.zeros_like(data)

        dv = np.sqrt((1+hdr['BERV']/299792.458) /(1-hdr['BERV']/299792.458) )
        wave = fits2wave(hdr)*dv
        model = np.zeros_like(data)

        # we project the template onto the wavelength grid and mask pixels for which we do not have a reliable
        # template
        for iord in tqdm(range(0,49), leave = False):
            mask = mask_template(wave[iord])
            tmp = template(wave[iord]) * blaze[iord]
            tmp[mask < 0.99] = np.nan
            model[iord] = tmp

        # keep track of the median value of each order
        med = np.zeros(49)
        for iord in range(49):
            med[iord] = np.nanmedian(data[iord])
            lowf = medfilt(data[iord]/model[iord],101)
            data[iord]/=lowf

        # loop through orders and replace large sigma excursions by template
        for iord in tqdm(range(0,49), leave = False):
            data_tmp = np.array(data[iord])

            if np.max(np.isfinite(data[iord])) == False:
                data[iord] = np.nan
                continue

            model_tmp = np.array(model[iord])


            residual = data_tmp - model_tmp
            sig = running_sigma(residual, 101)
            with warnings.catch_warnings(record=True) as _:
                nsig = residual / sig
            prob_good = np.exp( -0.5*nsig**2)
            pp = prob_good * (1+prob_bad) / (prob_bad+prob_good)

            data[iord] = data_tmp*pp+model_tmp*(1-pp)
            data[iord][np.isnan(data_tmp)] = model_tmp[np.isnan(data_tmp)]

            min_valid = np.min( np.where(np.isfinite(data_tmp)) )
            max_valid = np.max( np.where(np.isfinite(data_tmp)) )
            nanmask = np.zeros_like(data_tmp) + np.nan
            nanmask[min_valid:max_valid] = 0

            residual_after = data[iord] - model_tmp
            with warnings.catch_warnings(record=True) as _:
                nsig_after = residual_after / sig

            residual += nanmask
            sig += nanmask
            nsig += nanmask
            nsig_after += nanmask
            prob_good += nanmask
            pp += nanmask
            data[iord] += nanmask
            model_tmp += nanmask


            if iord in doplot:
                fig, ax = plt.subplots(nrows =3, ncols = 1, sharex = True)
                ax[1].plot(wave[iord],data_tmp,color = 'red', label = 'observed')
                ax[1].plot(wave[iord],data[iord], color = 'blue', label = 'sanitized')
                ax[0].set(xlabel = 'wavelength (nm)', ylabel = 'nsig')
                ax[1].plot(wave[iord],model_tmp,color = 'orange', label = 'model',alpha = 0.5)
                ax[1].set(xlabel = 'wavelength (nm)', ylabel = 'flux')
                ax[1].legend()
                ax[0].plot(wave[iord],nsig,label = 'nsig  before', alpha = .8)
                ax[0].plot(wave[iord],nsig_after  ,label = 'nsig  after',alpha = .8)

                ax[2].plot(wave[iord][min_valid:max_valid],pp[min_valid:max_valid], color = 'green', label = 'prob' )
                ax[2].legend()

                plt.show()
                plt.close()

            tmp = data[iord]
            tmp[0:min_valid] = np.nan
            tmp[max_valid:] = np.nan
            data[iord] = tmp
            prob_map[iord] = pp

            fgood[iord] = np.nanmean(pp[min_valid:max_valid])
            str_iord = str(iord)

            if len(str_iord) <2:
                str_iord = '0'+str_iord

            if np.isfinite(fgood[iord]) == 0:
                fgood[iord] = 0
            hdr['PSANI'+str_iord] = fgood[iord],'Order {0} mean(p) sanitize'.format(iord)

        if (doplot[0] !=  (-1)):
            plt.plot(1-fgood,'go')
            plt.xlabel('Nth order')
            plt.ylabel('Model padding fraction')
            plt.show()

        prob_map[np.isfinite(prob_map) == False] = 0
        model[~np.isfinite(data)] = np.nan


        for iord in range(49):
            data[iord] /=np.nanmedian(data[iord])
            data[iord]*=med[iord]

        # produce a template e2ds
        for iord in range(49):
            if np.max(np.isfinite(data[iord])) == False:
                continue
            model[iord]/=np.nanmedian(model[iord])
            model[iord]*=np.nanmedian(data[iord])


        data[np.isfinite(data) == False] = np.nan

        fits.writeto(outname, data, hdr, overwrite = True)
        # we also write the e2ds file for the model. This can be used as a sanity check and to determine telluric
        # correction accuracy.
        fits.writeto('model'.join(outname.split('sani')), model, hdr, overwrite = True)




"""
if True:
    #obj = 'TOI-1452'
    #obj = 'Gl514'
    obj = 'Gl699'
    files = np.array(glob.glob(obj+'/2??????o_pp_e2dsff_tcorr_AB.fits'))
    files = files[np.argsort(np.random.random(len(files)))]
    blaze_file = 'varia/2019-04-20_2400404f_pp_blaze_AB.fits'
    force = False
    doplot = [-1]#[5,20,35,42,48]
    prob_bad = 1e-4
    model_s1d_file = 'varia/Template_s1d_'+obj+'_sc1d_v_file_AB.fits'
    sanitize(files, blaze_file, model_s1d_file, force = force, doplot = doplot, prob_bad = prob_bad)
"""