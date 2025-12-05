import numpy as np
from astropy.io import fits
from astropy.table import Table
from fits2wave import fits2wave
from scipy.interpolate import InterpolatedUnivariateSpline
import os
import warnings

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


def sed_match(files, blaze_file, model_s1d_file, filter_width = 101, force = False, verbose = False):
    #
    # Program to create *sed* files, these are identical to e2ds_tcorr, but the SED it forces to
    # match that of template
    #
    # Sample inputs --
    #
    # Files to be sed_corrected, needs to be torr
    #    files = glob.glob('Gl699/2*o_pp_e2dsff_tcorr_AB.fits')
    #
    # Blaze file
    #    blaze_file = 'varia/2019-04-20_2400404f_pp_blaze_AB.fits'
    #
    # s1d_v model file for the object
    #    model_s1d_file = 'varia/Template_s1d_Gl699_sc1d_v_file_AB.fits'
    #
    # filter_width is the width over which we smooth the spectrum to match its SED to that of the template
    #
    # sed_match(files, blaze_file,model_s1d_file)
    #
    # optional inputs --
    #
    #   Force
    #     force rewriting of files. Useful if you modify the code

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


    for i in (range(len(files))):
        # loop through files
        if verbose:
            print('Matching SEDs #{0} in {1} -> {2}'.format(i+1,len(files),files[i]))
        outname = 'sed'.join(files[i].split('tcorr'))

        if os.path.isfile(outname) and (force == False):
            if verbose:
                print(outname+' exists, we skip')
            continue

        data, hdr = fits.getdata(files[i], header = True)

        with warnings.catch_warnings(record=True) as _:
            data[blaze < 0.3] = np.nan

        dv = np.sqrt((1+hdr['BERV']/299792.458) /(1-hdr['BERV']/299792.458) )
        wave = fits2wave(hdr)*dv
        model = np.zeros_like(data)

        # we project the template onto the wavelength grid and mask pixels for which we do not have a reliable
        # template
        for iord in range(0,49):
            mask = mask_template(wave[iord])
            tmp = template(wave[iord]) * blaze[iord]
            tmp[mask < 0.99] = np.nan
            model[iord] = tmp

        for iord in range(49):

            with warnings.catch_warnings(record=True) as _:
                lowf = lowpassfilter(data[iord]/model[iord],filter_width)

            key = 'SEDCOR'+str(iord).zfill(2)

            med = np.nanmedian(lowf)
            if np.isfinite(med) == False:
                med = 0.0
            hdr[key] = med

            data[iord]/=lowf

        # put the SED-corrected spectrum back into pace
        data[np.isfinite(data) == False] = np.nan

        fits.writeto(outname, data, hdr, overwrite = True)