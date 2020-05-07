import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline
from fits2wave import fits2wave

def get_dvrms(template_file, blaze_file,sample_file):
    # Provide a template and get the estimate of the RV accuracy per order
    #
    # This vector can be used as a constant weighting factor for CCFs
    #
    # Inputs :
    #
    # Template_file -> Template in s1d format for that object
    #
    # blaze_file -> any blaze file that is for that fiber. Only very
    #               weakly dependent on the shape of the blaze
    #
    # sample_file -> tcorr file for this object. HDR used to get the wavelength
    #                grid in the computation and proper colors to get the
    #                observed (weighted by instrument throughput) SED
    #
    # returns:
    #
    # A per-order estimate of the RV accuracy assuming a photon-noise-limited regime.
    #
    # This value can be used as a weight for the CCF using weight = 1/dv**2
    #
    blaze, hdr = fits.getdata(blaze_file, header = True)
    sample_sp, sample_hdr = fits.getdata(sample_file, header = True)
    wave = fits2wave(sample_hdr)
    template = Table(fits.getdata(template_file))

    c = 299792458 # in m/s

    # per-order rv accuracy. If not valid, will be set to inf
    dv = np.zeros(wave.shape[0]) +  np.inf

    for iord in range(wave.shape[0]):
        # normalize blaze
        blaze[iord] /= np.nanpercentile(blaze[iord],95)
        blaze_tmp = np.array(blaze[iord])
        blaze_tmp[~np.isfinite(blaze_tmp)] = 0

        # if no more than 50% of pixels above a blaze of 0.3 are valid in the
        # template, then skip order
        if np.nanmean(np.isfinite(sample_sp[iord][blaze_tmp>0.3]))<0.5:
            continue

        # get min/max of wavelength domain
        wmin = np.min(wave[iord])
        wmax = np.max(wave[iord])
        # valid points in the tempalte
        g = np.isfinite(template['flux']) & (template['wavelength']>wmin) &  (template['wavelength']<wmax)
        spline = InterpolatedUnivariateSpline(template['wavelength'][g],template['flux'][g], k = 1, ext=3)

        # spline onto that file's wavelength grid. Note that we don't really need to
        # correct for the BERV. The BERV shifts a tad bit the domain, but this is
        # only the integral, which doesn't really change much when you shift a bit.
        tmp = spline(wave[iord])*blaze[iord]
        tmp /= np.nanmedian(tmp)

        # normalize the model and scale to observation
        tmp*=np.nanmean(sample_sp[iord])

        # Equations (5) and (6) in Bouchy 2001
        dv_pix = c*np.sqrt(tmp)/(  wave[iord]*np.gradient(tmp)/ np.gradient(wave[iord]))
        dv[iord] = 1/np.sqrt(np.nansum(1/dv_pix**2))

    return dv

