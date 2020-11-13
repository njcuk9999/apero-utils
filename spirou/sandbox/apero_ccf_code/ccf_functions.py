#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CODE DESCRIPTION HERE

Created on 2020-11-2020-11-05 10:21

@author: cook
"""
import argparse
from astropy import constants as cc
from astropy import units as uu
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import copy
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from scipy.interpolate import InterpolatedUnivariateSpline
import sys
from typing import Union
import warnings
import yaml


# =============================================================================
# Define variables
# =============================================================================
__NAME__ = 'ccf_functions.py'
# Get version and author
__version__ = '0.6.132-ccf2'
__author__ = 'Neil Cook'
__date__ = '2020-11-09'
# Speed of light
# noinspection PyUnresolvedReferences
speed_of_light_ms = cc.c.to(uu.m / uu.s).value
# noinspection PyUnresolvedReferences
speed_of_light = cc.c.to(uu.km / uu.s).value
# set up time
now = Time.now()


# =============================================================================
# Define general functions
# =============================================================================
class AperoError(Exception):
    pass


def get_params(**kwargs):
    # start arg parser
    desc = ('CCF Code to emulate APERO cal_ccf_spirou.py')
    parser = argparse.ArgumentParser(description=desc)
    # add all keys
    parser.add_argument('--config', dest='config', default='ccf.yaml',
                        help='The configuration yaml to use.')
    # parse arguments
    args = parser.parse_args()
    # get config file
    config_file = kwargs.get('config', args.config)
    # deal with bad config_file
    if not os.path.exists(config_file):
        config_file = 'ccf.yaml'
    # load params from constants
    with open(config_file, 'r') as yfile:
        params = yaml.load(yfile, Loader=yaml.FullLoader)
    # loop around and add kwargs
    for key in params:
        # deal with key in kwargs
        if key in kwargs:
            params[key] = kwargs[key]
    # add config file to args
    params['CONFIG'] = config_file
    # return params
    return params


def get_berv(params, header):
    estimate = params['EXT_BERV_EST_ACC']
    bprops = dict()
    try:
        bprops['BERV'] = float(header['BERV'])
        bprops['BJD'] = float(header['BJD'])
        bprops['BERV_MAX'] = float(header['BERVMAX'])
        bprops['DBERV'] = float(header['DBERV'])
        bprops['BERV_SOURCE'] = header['BERVSRCE']
        bprops['BERV_EST'] = float(header['BERV_EST'])
        bprops['BJD_EST'] = float(header['BJD_EST'])
        bprops['BERV_MAX_EST'] = float(header['BERVMAXE'])
        bprops['DBERV_EST'] = float(header['DBERVE'])
        bprops['OBS_TIME'] = header['BERVOBST']
        bprops['OBS_TIME_METHOD'] = header['BERVOBSM']
    except Exception as e:
        wlog(params, 'error', str(e))
    # -------------------------------------------------------------------------
    # need to decide which values should be used (and report if we are using
    #   estimate)
    cond = (bprops['BERV'] is not None and np.isfinite(bprops['BERV']))
    cond &= (bprops['BJD'] is not None and np.isfinite(bprops['BJD']))
    cond &= (bprops['BERV_MAX'] is not None and np.isfinite(bprops['BERV_MAX']))

    if not cond or (bprops['BERV_SOURCE'] == 'pyasl'):
        # log warning that we are using an estimate
        wmsg = ("No barycorrpy berv values found, using pyasl estimate. "
                "\n\t Only good to  +/- {0:.3f} m/s "
                "\n\t For precise RV must have barycorrpy working "
                "(see extraction log)")
        wlog(params, 'warning', wmsg.format(estimate))
        # set parameters
        bprops['USE_BERV'] = bprops['BERV_EST']
        bprops['USE_BJD'] = bprops['BJD_EST']
        bprops['USE_BERV_MAX'] = bprops['BERV_MAX_EST']
        bprops['USED_ESTIMATE'] = True
    # Case 3: Barycorrpy used
    else:
        # set parameters
        bprops['USE_BERV'] = bprops['BERV']
        bprops['USE_BJD'] = bprops['BJD']
        bprops['USE_BERV_MAX'] = bprops['BERV_MAX']
        bprops['USED_ESTIMATE'] = False
    # return bprops
    return bprops


def get_wavesolution(params, header, fiber='AB'):
    wprops = dict()
    # deal with fiber
    wavebasefile = params['WAVEFILE_{0}'.format(fiber)]
    # construct wave file
    wavefile = os.path.join(params['WORKING_DIR'], wavebasefile)
    # set source
    wavesource = 'input'
    # deal with not existing
    if not os.path.exists(wavefile):
        # get working directory
        wdir = params['WORKING_DIR']
        # construct wave file
        wavefile = os.path.join(wdir, header['CDBWAVE'])
        # set source
        wavesource = 'header'
    # deal with still not existing
    if not os.path.exists(wavefile):
        emsg = 'Wave file {0} not found (source={1})'
        wlog(params, 'error', emsg.format(wavefile, wavesource))
    # load map
    wavemap, wheader = fits.getdata(wavefile, header=True)
    # get shape
    nbo, nbpix = wavemap.shape
    # fill wprops
    wprops['WAVEFILE'] = wavebasefile
    wprops['WAVEINIT'] = wavebasefile
    wprops['WAVESOURCE'] = wavesource
    wprops['WFP_DRIFT'] = float(wheader['WFPDRIFT'])
    wprops['NBO'] = nbo
    wprops['DEG'] = int(wheader['WAVEDEGN'])
    wprops['COEFFS'] = None
    wprops['WAVEMAP'] = wavemap
    wprops['WAVEINST'] = None
    wprops['WAVETIME'] = float(wheader['MJDMID'])
    # return wprops
    return wprops


def get_blaze(params, header):
    # construct wave file
    blazefile = os.path.join(params['WORKING_DIR'], params['BLAZEFILE'])
    # set source
    blazesource = 'input'
    # deal with not existing
    if not os.path.exists(blazefile):
        # get working directory
        wdir = params['WORKING_DIR']
        # construct wave file
        blazefile = os.path.join(wdir, header['CDBWAVE'])
        # set source
        blazesource = 'header'
    # deal with not existing
    if not os.path.exists(blazefile):
        emsg = 'Blaze file {0} not found (source={1})'
        wlog(params, 'error', emsg.format(blazefile, blazesource))
    # load map
    blaze = fits.getdata(blazefile)
    # return blaze file
    return blazefile, blaze


def wlog(params, level, msg):
    # get colours
    green = params['COLORS.GREEN']
    yellow = params['COLORS.YELLOW']
    blue = params['COLORS.BLUE']
    red = params['COLORS.RED']
    endc = params['COLORS.ENDC']
    # add time stamp
    now = str(Time.now().iso).split(' ')[-1]
    msg_fmt = '{0}{1}|{2}|{3}{4}'
    # deal with error
    if level == 'error':
        message = msg_fmt.format(red, now, '!', msg, endc)
        print(message)
        raise AperoError(msg)
    # deal with warning
    elif level == 'warning':
        message = msg_fmt.format(yellow, now, '@', msg, endc)
        print(message)
    # deal with info
    elif level == 'info':
        message = msg_fmt.format(blue, now, '*', msg, endc)
        print(message)
    # deal with all other messages
    else:
        message = msg_fmt.format(green, now, '', msg, endc)
        print(message)


def locate_reference_file(params):
    # get infile from params
    infile = os.path.join(params['WORKING_DIR'], params['INFILE_C'])
    # check recon file exists
    if not os.path.exists(infile):
        emsg = 'Cannot load infile C file {0}'
        wlog(params, 'error', emsg.format(infile))
    # return filename
    return infile


def get_ccf_mask(params, filename, mask_width, mask_units='nm'):
    func_name = __NAME__ + '.get_ccf_mask()'
    # load table
    table = Table.read(filename, format=params['CCF_MASK_FMT'])
    colnames = ['ll_mask_s', 'll_mask_e', 'w_mask']
    oldcols = table.colnames
    for c_it, col in enumerate(colnames):
        table[oldcols[c_it]].name = col
    # convert to floats
    ll_mask_e = np.array(table['ll_mask_e']).astype(float)
    ll_mask_s = np.array(table['ll_mask_s']).astype(float)
    # calculate the difference in mask_e and mask_s
    ll_mask_d = ll_mask_e - ll_mask_s
    ll_mask_ctr = ll_mask_s + ll_mask_d * 0.5
    # if mask_width > 0 ll_mask_d is multiplied by mask_width/c
    if mask_width > 0:
        ll_mask_d = mask_width * ll_mask_s / speed_of_light
    # make w_mask an array
    w_mask = np.array(table['w_mask']).astype(float)
    # ----------------------------------------------------------------------
    # deal with the units of ll_mask_d and ll_mask_ctr
    # must be returned in nanometers
    # ----------------------------------------------------------------------
    # get unit object from mask units string
    try:
        unit = getattr(uu, mask_units)
    except Exception as e:
        # log error
        eargs = [mask_units, type(e), e, func_name]
        emsg = ("Fiber={0} is not valid for 'reference'. (Required '{1}') "
                "\n\t File {2}: {3} \n\t Function = {4}")
        wlog(params, 'error', emsg.format(*eargs))
        return None, None, None
    # add units
    ll_mask_d = ll_mask_d * unit
    ll_mask_ctr = ll_mask_ctr * unit
    # convert to nanometers
    ll_mask_d = ll_mask_d.to(uu.nm).value
    ll_mask_ctr = ll_mask_ctr.to(uu.nm).value
    # ----------------------------------------------------------------------
    # return the size of each pixel, the central point of each pixel
    #    and the weight mask
    return ll_mask_d, ll_mask_ctr, w_mask


def relativistic_waveshift(dv, units='km/s'):
    """
    Relativistic offset in wavelength

    default is dv in km/s

    :param dv: float or numpy array, the dv values
    :param units: string or astropy units, the units of dv
    :return:
    """
    # get c in correct units
    # noinspection PyUnresolvedReferences
    if units == 'km/s' or units == uu.km / uu.s:
        c = speed_of_light
    # noinspection PyUnresolvedReferences
    elif units == 'm/s' or units == uu.m / uu.s:
        c = speed_of_light_ms
    else:
        raise ValueError("Wrong units for dv ({0})".format(units))
    # work out correction
    corrv = np.sqrt((1 + dv / c) / (1 - dv / c))
    # return correction
    return corrv


class NanSpline:
    def __init__(self, emsg: str, x: Union[np.ndarray, None] = None,
                 y: Union[np.ndarray, None] = None, **kwargs):
        """
        This is used in place of scipy.interpolateInterpolatedUnivariateSpline
        (Any spline following this will return all NaNs)

        :param emsg: str, the error that means we have to use the NanSpline
        """
        self.emsg = str(emsg)
        self.x = copy.deepcopy(x)
        self.y = copy.deepcopy(y)
        self.kwargs = copy.deepcopy(kwargs)

    def __repr__(self) -> str:
        """
        String representation of the class

        :return: string representation
        """
        return 'NanSpline: \n' + self.emsg

    def __str__(self) -> str:
        """
        String representation of the class

        :return: string representation
        """
        return self.__repr__()

    def __call__(self, x: np.ndarray, nu: int = 0,
                 ext: Union[int, None] = None) -> np.ndarray:
        """
        Return a vector of NaNs (this means the spline failed due to
        less points

        This is used in place of scipy.interpolateInterpolatedUnivariateSpline

        Parameters
        ----------
        x : array_like
            A 1-D array of points at which to return the value of the smoothed
            spline or its derivatives. Note: x can be unordered but the
            evaluation is more efficient if x is (partially) ordered.
        nu  : int
            The order of derivative of the spline to compute.
        ext : int
            Controls the value returned for elements of ``x`` not in the
            interval defined by the knot sequence.

            * if ext=0 or 'extrapolate', return the extrapolated value.
            * if ext=1 or 'zeros', return 0
            * if ext=2 or 'raise', raise a ValueError
            * if ext=3 or 'const', return the boundary value.

            The default value is 0, passed from the initialization of
            UnivariateSpline.

        """
        return np.repeat(np.nan, len(x))


def iuv_spline(x, y, **kwargs):
    # copy x and y
    x, y = np.array(x), np.array(y)
    # find all NaN values
    nanmask = ~np.isfinite(y)

    # deal with dimensions error (on k)
    #   otherwise get   dfitpack.error: (m>k) failed for hidden m
    if kwargs.get('k', None) is not None:
        # deal with to few parameters in x
        if len(x[~nanmask]) < (kwargs['k'] + 1):
            # raise exception if len(x) is bad
            emsg = ('IUV Spline len(x) < k+1 '
                    '\n\tk={0}\n\tlen(x) = {1}'
                    '\n\tx={2}\n\ty={3}')
            eargs = [kwargs['k'], len(x), str(x)[:70], str(y)[:70]]
            # return a nan spline
            return NanSpline(emsg.format(*eargs), x=x, y=y, **kwargs)
        if len(y[~nanmask]) < (kwargs['k'] + 1):
            # raise exception if len(x) is bad
            emsg = ('IUV Spline len(y) < k+1 '
                    '\n\tk={0}\n\tlen(y) = {1}'
                    '\n\tx={2}\n\ty={3}')
            eargs = [kwargs['k'], len(y), str(x)[:70], str(y)[:70]]
            # return a nan spline
            return NanSpline(emsg.format(*eargs), x=x, y=y, **kwargs)

    if np.sum(~nanmask) < 2:
        y = np.repeat(np.nan, len(x))
    elif np.sum(~nanmask) == 0:
        y = np.repeat(np.nan, len(x))
    else:
        # replace all NaN's with linear interpolation
        badspline = InterpolatedUnivariateSpline(x[~nanmask], y[~nanmask],
                                                 k=1, ext=1)
        y[nanmask] = badspline(x[nanmask])
    # return spline
    return InterpolatedUnivariateSpline(x, y, **kwargs)


def fwhm(sigma=1.0):
    """
    Get the Full-width-half-maximum value from the sigma value (~2.3548)

    :param sigma: float, the sigma, default value is 1.0 (normalised gaussian)
    :return: 2 * sqrt(2 * log(2)) * sigma = 2.3548200450309493 * sigma
    """
    return 2 * np.sqrt(2 * np.log(2)) * sigma


def gauss_function(x, a, x0, sigma, dc):
    """
    A standard 1D gaussian function (for fitting against)

    :param x: numpy array (1D), the x data points
    :param a: float, the amplitude
    :param x0: float, the mean of the gaussian
    :param sigma: float, the standard deviation (FWHM) of the gaussian
    :param dc: float, the constant level below the gaussian

    :return gauss: numpy array (1D), size = len(x), the output gaussian
    """
    return a * np.exp(-0.5 * ((x - x0) / sigma) ** 2) + dc


def fitgaussian(x, y, weights=None, guess=None, return_fit=True,
                return_uncertainties=False):
    """
    Fit a single gaussian to the data "y" at positions "x", points can be
    weighted by "weights" and an initial guess for the gaussian parameters

    :param x: numpy array (1D), the x values for the gaussian
    :param y: numpy array (1D), the y values for the gaussian
    :param weights: numpy array (1D), the weights for each y value
    :param guess: list of floats, the initial guess for the guassian fit
                  parameters in the following order:

                  [amplitude, center, fwhm, offset from 0 (in y-direction)]

    :param return_fit: bool, if True also calculates the fit values for x
                       i.e. yfit = gauss_function(x, *pfit)

    :param return_uncertainties: bool, if True also calculates the uncertainties
                                 based on the covariance matrix (pcov)
                                 uncertainties = np.sqrt(np.diag(pcov))

    :return pfit: numpy array (1D), the fit parameters in the
                  following order:

                [amplitude, center, fwhm, offset from 0 (in y-direction)]

    :return yfit: numpy array (1D), the fit y values, i.e. the gaussian values
                  for the fit parameters, only returned if return_fit = True

    """

    # if we don't have weights set them to be all equally weighted
    if weights is None:
        weights = np.ones(len(x))
    weights = 1.0 / weights
    # if we aren't provided a guess, make one
    if guess is None:
        guess = [np.nanmax(y), np.nanmean(y), np.nanstd(y), 0]
    # calculate the fit using curve_fit to the function "gauss_function"
    with warnings.catch_warnings(record=True) as _:
        pfit, pcov = curve_fit(gauss_function, x, y, p0=guess, sigma=weights,
                               absolute_sigma=True)
    if return_fit and return_uncertainties:
        # calculate the fit parameters
        yfit = gauss_function(x, *pfit)
        # work out the normalisation constant
        chis, _ = chisquare(y, f_exp=yfit)
        norm = chis / (len(y) - len(guess))
        # calculate the fit uncertainties based on pcov
        efit = np.sqrt(np.diag(pcov)) * np.sqrt(norm)
        # return pfit, yfit and efit
        return pfit, yfit, efit
    # if just return fit
    elif return_fit:
        # calculate the fit parameters
        yfit = gauss_function(x, *pfit)
        # return pfit and yfit
        return pfit, yfit
    # if return uncertainties
    elif return_uncertainties:
        # calculate the fit parameters
        yfit = gauss_function(x, *pfit)
        # work out the normalisation constant
        chis, _ = chisquare(y, f_exp=yfit)
        norm = chis / (len(y) - len(guess))
        # calculate the fit uncertainties based on pcov
        efit = np.sqrt(np.diag(pcov)) * np.sqrt(norm)
        # return pfit and efit
        return pfit, efit
    # else just return the pfit
    else:
        # return pfit
        return pfit


def write_file(params, props, infile, maskname, header, wprops):
    # ----------------------------------------------------------------------
    # construct out file name
    inbasename = os.path.basename(infile).split('.')[0]
    maskbasename = os.path.basename(maskname).split('.')[0]
    inpath = os.path.dirname(infile)
    outfile = 'CCFTABLE_{0}_{1}.fits'.format(inbasename, maskbasename)
    outpath = os.path.join(inpath, outfile)
    # ----------------------------------------------------------------------
    # produce CCF table
    table1 = Table()
    table1['RV'] = props['RV_CCF']
    for order_num in range(len(props['CCF'])):
        table1['ORDER{0:02d}'.format(order_num)] = props['CCF'][order_num]
    table1['COMBINED'] = props['MEAN_CCF']
    # ----------------------------------------------------------------------
    # produce stats table
    table2 = Table()
    table2['ORDERS'] = np.arange(len(props['CCF'])).astype(int)
    table2['NLINES'] = props['CCF_LINES']
    # get the coefficients
    coeffs = props['CCF_FIT_COEFFS']
    table2['CONTRAST'] = np.abs(100 * coeffs[:, 0])
    table2['RV'] = coeffs[:, 1]
    table2['FWHM'] = coeffs[:, 2]
    table2['DC'] = coeffs[:, 3]
    table2['SNR'] = props['CCF_SNR']
    table2['NORM'] = props['CCF_NORM']
    table2['DVRMS_SP'] = props['ORD_SPEC_RMS']
    table2['DVRMS_CC'] = props['CCF_NOISE']

    # ----------------------------------------------------------------------
    # add to the header
    # ----------------------------------------------------------------------
    # add results from the CCF
    header['CCFMNRV'] = (props['MEAN_RV'],
                         'Mean RV calc. from the mean CCF [km/s]')
    header['CCFMCONT'] = (props['MEAN_CONTRAST'],
                          'Mean contrast (depth of fit) from mean CCF')
    header['CCFMFWHM'] = (props['MEAN_FWHM'],
                          'Mean FWHM from mean CCF')
    header['CCFMRVNS'] = (props['MEAN_RV_NOISE'],
                          'Mean RV Noise from mean CCF')
    header['CCFTLINE'] = (props['TOT_LINE'],
                          'Total no. of mask lines used in CCF')
    # ----------------------------------------------------------------------
    # add constants used to process
    header['CCFMASK'] = (props['CCF_MASK'], 'CCF mask file used')
    header['CCFSTEP'] = (props['CCF_STEP'], 'CCF step used [km/s]')
    header['CCFWIDTH'] = (props['CCF_WIDTH'], 'CCF width used [km/s]')
    header['CCFTRGRV'] = (props['TARGET_RV'],
                          'CCF central RV used in CCF [km/s]')
    header['CCFSIGDT'] = (props['CCF_SIGDET'],
                          'Read noise used in photon noise calc. in CCF')
    header['CCFBOXSZ'] = (props['CCF_BOXSIZE'],
                          'Size of bad px used in photon noise calc. in CCF')
    header['CCFMAXFX'] = (props['CCF_MAXFLUX'],
                          'Flux thres for bad px in photon noise calc. in CCF')
    header['CCFORDMX'] = (props['CCF_NMAX'],
                          'Last order used in mean for mean CCF')
    header['CCFMSKMN'] = (props['MASK_MIN'],
                          'Minimum weight of lines used in the CCF mask')
    header['CCFMSKWD'] = (props['MASK_WIDTH'],
                          'Width of lines used in the CCF mask')
    header['CCFMUNIT'] = (props['MASK_UNITS'], 'Units used in CCF Mask')
    # ----------------------------------------------------------------------
    header['RV_WAVFN'] = (os.path.basename(wprops['WAVEFILE']),
                          'RV wave file used')
    header['RV_WAVTM'] = (wprops['WAVETIME'],
                          'RV wave file time [mjd]')
    header['RV_WAVTD'] = (header['MJDMID'] - wprops['WAVETIME'],
                          'RV timediff [days] btwn file and wave solution')
    header['RV_WAVFP'] = ('None', 'RV measured from wave sol FP CCF [km/s]')
    header['RV_SIMFP'] = ('None', 'RV measured from simultaneous FP CCF [km/s]')
    header['RV_DRIFT'] = ('None',
                          'RV drift between wave sol and sim. FP CCF [km/s]')
    header['RV_OBJ'] = (props['MEAN_RV'],
                        'RV calc in the object CCF (non corr.) [km/s]')
    header['RV_CORR'] = ('None', 'RV corrected for FP CCF drift [km/s]')
    # ----------------------------------------------------------------------
    # log where we are writing the file to
    wlog(params, '', 'Writing file to {0}'.format(outpath))
    # construct hdus
    hdu = fits.PrimaryHDU()
    t1 = fits.BinTableHDU(table1, header=header)
    t2 = fits.BinTableHDU(table2, header=header)
    # construct hdu list
    hdulist = fits.HDUList([hdu, t1, t2])
    # write hdulist
    hdulist.writeto(outpath, overwrite=True)


# =============================================================================
# Define rv functions
# =============================================================================
def remove_telluric_domain(params, infile):
    # get parameters from params/kwargs
    ccf_tellu_thres = params['CCF_TELLU_THRES']
    # get the infile image
    image = fits.getdata(infile)
    # get the image
    reconfile = os.path.join(params['WORKING_DIR'], params['RECONFILE'])
    # check recon file exists
    if not os.path.exists(reconfile):
        eargs = [infile, 'recon', reconfile]
        emsg = ("Cannot remove telluric domain for file: {0} "
                "\n\t The {1} file was not found: {2}")
        wlog(params, 'error', emsg.format(*eargs))
    # read recon file
    recondata = fits.getdata(reconfile)
    # find all places below threshold
    with warnings.catch_warnings(record=True) as _:
        keep = recondata > ccf_tellu_thres
    # set all bad data to NaNs
    image[~keep] = np.nan
    # return in file
    return image


def compute_ccf_science(params, infile, image, header, blaze, wavemap,
                        bprops, fiber):
    func_name = __NAME__ + '.compute_ccf()'
    # get parameters from params/kwargs
    noise_sigdet = params['CCF_NOISE_SIGDET']
    noise_size = params['CCF_NOISE_BOXSIZE']
    noise_thres = params['CCF_NOISE_THRES']
    mask_min = params['CCF_MASK_MIN_WEIGHT']
    mask_width = params['CCF_MASK_WIDTH']
    mask_units = params['CCF_MASK_UNITS']
    fit_type = params['CCF_FIT_TYPE']
    ccfnmax = params['CCF_N_ORD_MAX']
    null_targetrv = params['CCF_OBJRV_NULL_VAL']
    maxwsr = params['CCF_MAX_CCF_WID_STEP_RATIO']

    # get image size
    nbo, nbpix = image.shape
    # get parameters from inputs
    ccfstep = params['CCF_STEP']
    ccfwidth = params['CCF_WIDTH']
    targetrv = params['CCF_RV']
    # ----------------------------------------------------------------------
    # TODO: eventually this should come from object database (so that each
    # TODO: object has a constant target rv
    # need to deal with no target rv step
    if np.isnan(targetrv):
        targetrv = header.get('OBJRV', None)
        # targetrv must be a float (could be a "NaN")
        if targetrv is not None:
            targetrv = float(targetrv)
        # set target rv to zero if we don't have a value
        if targetrv is None:
            wargs = [params['KW_INPUTRV'][0], infile]
            wmsg = ("RV not found in header. Setting input RV to 0.0 "
                    "\n\t Key={0}) \n\t File = {1}")
            wlog(params, 'warning', wmsg.format(*wargs))
            targetrv = 0.0
        elif np.abs(targetrv) > null_targetrv:
            wargs = [params['KW_INPUTRV'][0], null_targetrv, infile]
            wmsg = ("RV null value found in header. Setting input RV to 0.0 "
                    "\n\t Key={0} Nullvalue={1} \n\t File = {2}")
            wlog(params, 'warning', wmsg.format(*wargs))
            targetrv = 0.0
    # ----------------------------------------------------------------------
    # need to deal with mask coming from inputs
    ccfmask = os.path.join(params['MASK_DIR'], params['CCF_MASK_AB'])
    # get the berv
    berv = bprops['USE_BERV']

    # ----------------------------------------------------------------------
    # Need some sanity checking on width and step
    # ----------------------------------------------------------------------
    if ccfstep > (ccfwidth / maxwsr):
        eargs = [ccfwidth, ccfstep, maxwsr, func_name]
        emsg = ("CCF step much not be greater than CCF width / {2}. "
                "\n\t Current CCF step = {0} \n\t Current CCF width = {1} "
                "\n\t Function = {3}")
        wlog(params, 'error', emsg.format(*eargs))

    # ----------------------------------------------------------------------
    # Check we are using correct fiber
    # ----------------------------------------------------------------------
    sfiber, rfiber = params['FIBER_CCF']
    if fiber != sfiber:
        # log that the science fiber was not correct
        eargs = [fiber, sfiber, 'INFILE_AB', infile]
        emsg = ("Fiber={0} is not valid for 'science'. (Required '{1}') "
                "\n\t File {2}: {3} \n\t Function = {4}")
        wlog(params, 'error', emsg.format(*eargs))

    # ----------------------------------------------------------------------
    # Compute photon noise uncertainty for reference file
    # ----------------------------------------------------------------------
    # set up the arguments for DeltaVrms2D
    dkwargs = dict(spe=image, wave=wavemap, sigdet=noise_sigdet,
                   size=noise_size, threshold=noise_thres)
    # run DeltaVrms2D
    dvrmsref, wmeanref, wmeanrefo = delta_v_rms_2d(**dkwargs)
    # log the estimated RV uncertainty
    wargs = [fiber, wmeanref]
    wmsg = "On fiber {0} estimated RV uncertainty on spectrum is {1:.3f}"
    wlog(params, 'info', wmsg.format(*wargs))
    # ----------------------------------------------------------------------
    # Do the CCF calculations
    # ----------------------------------------------------------------------
    # get the mask parameters
    mkwargs = dict(filename=ccfmask, mask_width=mask_width,
                   mask_units=mask_units)
    ll_mask_d, ll_mask_ctr, w_mask = get_ccf_mask(params, **mkwargs)

    # calculate the CCF
    props = ccf_calculation(params, image, blaze, wavemap, berv, targetrv,
                            ccfwidth, ccfstep, ll_mask_ctr, w_mask,
                            fit_type, fiber)
    # ----------------------------------------------------------------------
    # Reference plots
    # ----------------------------------------------------------------------
    if params['PLOT_ALL']:
        # the image vs wavelength for an order
        plot_ccf_swave_ref(wavemap=wavemap, image=image, fiber=fiber,
                           nbo=nbo)

    if params['PLOT']:
        # the image vs wavelength for an order
        plot_ccf_swave_ref(wavemap=wavemap, image=image, fiber=fiber,
                           nbo=nbo, order=35)
        # the photon noise uncertainty plot
        plot_ccf_photon_uncert(x=np.arange(nbo), y_sp=wmeanrefo,
                               y_cc=props['CCF_NOISE'])
    # ----------------------------------------------------------------------
    # Calculate the mean CCF
    # ----------------------------------------------------------------------
    # get the average ccf
    mean_ccf = np.nanmean(props['CCF'][: ccfnmax], axis=0)

    # get the fit for the normalized average ccf
    mean_ccf_coeffs, mean_ccf_fit = fit_ccf(params, 'mean', props['RV_CCF'],
                                            mean_ccf, fit_type=fit_type)
    # get the max cpp
    # TODO: How do we calculate max_cpp and what is it? Do we need it?
    # max_cpp = mp.nansum(props['CCF_MAX']) / mp.nansum(props['PIX_PASSED_ALL'])
    # get the RV value from the normalised average ccf fit center location
    ccf_rv = float(mean_ccf_coeffs[1])
    # get the contrast (ccf fit amplitude)
    ccf_contrast = np.abs(100 * mean_ccf_coeffs[0])
    # get the FWHM value
    ccf_fwhm = mean_ccf_coeffs[2] * fwhm()
    # ----------------------------------------------------------------------
    #  combined CCF_NOISE uncertainty
    rv_noise = 1.0 / np.sqrt(np.nansum(1.0 / props['CCF_NOISE'] ** 2))
    # ----------------------------------------------------------------------
    # log the stats
    wargs = [ccf_contrast, float(mean_ccf_coeffs[1]), rv_noise, ccf_fwhm]
    wmsg = ("Correlation: C={0:1f}[%] RV={1:.5f}[km/s] RV_NOISE={2:.5f}[km/s] "
            "FWHM={3:.4f}[km/s]")
    wlog(params, 'info', wmsg.format(*wargs))
    # ----------------------------------------------------------------------
    # add to output array
    props['TOT_SPEC_RMS'] = wmeanref
    props['ORD_SPEC_RMS'] = wmeanrefo
    props['MEAN_CCF'] = mean_ccf
    props['MEAN_RV'] = ccf_rv
    props['MEAN_CONTRAST'] = ccf_contrast
    props['MEAN_FWHM'] = ccf_fwhm
    props['MEAN_CCF_RES'] = mean_ccf_coeffs
    props['MEAN_CCF_FIT'] = mean_ccf_fit
    props['MEAN_RV_NOISE'] = rv_noise
    # add constants to props
    props['CCF_MASK'] = ccfmask
    props['CCF_STEP'] = ccfstep
    props['CCF_WIDTH'] = ccfwidth
    props['TARGET_RV'] = targetrv
    props['CCF_SIGDET'] = noise_sigdet
    props['CCF_BOXSIZE'] = noise_size
    props['CCF_MAXFLUX'] = noise_thres
    props['CCF_NMAX'] = ccfnmax
    props['MASK_MIN'] = mask_min
    props['MASK_WIDTH'] = mask_width
    props['MASK_UNITS'] = mask_units
    # ------------------------------------------------------------------
    # rv ccf plot
    # ------------------------------------------------------------------
    # loop around every order
    if params['PLOT_ALL']:
        plot_ccf_rv_fit(params=params, x=props['RV_CCF'],
                        y=props['CCF'], yfit=props['CCF_FIT'], kind='SCIENCE',
                        rv=props['CCF_FIT_COEFFS'][:, 1], ccfmask=ccfmask,
                        orders=np.arange(len(props['CCF'])), order=None)

    if params['PLOT']:
        plot_ccf_rv_fit(params=params, x=props['RV_CCF'],
                        y=mean_ccf, yfit=mean_ccf_fit, kind='MEAN SCIENCE',
                        rv=ccf_rv, ccfmask=ccfmask,
                        orders=None, order=None)
    # ------------------------------------------------------------------
    # return property dictionary
    return props


def compute_ccf_fp(params, image, blaze, wavemap, fiber):
    # get constants from params/kwargs
    noise_sigdet = params['WAVE_CCF_NOISE_SIGDET']
    noise_size = params['WAVE_CCF_NOISE_BOXSIZE']
    noise_thres = params['WAVE_CCF_NOISE_THRES']
    ccfstep = params['WAVE_CCF_STEP']
    ccfwidth = params['WAVE_CCF_WIDTH']
    targetrv = params['WAVE_CCF_TARGET_RV']
    ccfmask = os.path.join(params['MASK_DIR'], params['CCF_MASK_C'])
    ccfnmax = params['WAVE_CCF_N_ORD_MAX']
    mask_min = params['WAVE_CCF_MASK_MIN_WEIGHT']
    mask_width = params['WAVE_CCF_MASK_WIDTH']
    mask_units = params['WAVE_CCF_MASK_UNITS']
    # set the berv to zero (fp have no berv)
    berv = 0
    # the fit type must be set to 1 (for emission lines)
    fit_type = 1
    # ------------------------------------------------------------------
    # Compute photon noise uncertainty for FP
    # ------------------------------------------------------------------
    # set up the arguments for DeltaVrms2D
    dkwargs = dict(spe=image, wave=wavemap, sigdet=noise_sigdet,
                   size=noise_size, threshold=noise_thres)
    # run DeltaVrms2D
    dvrmsref, wmeanref, wmeanrefo = delta_v_rms_2d(**dkwargs)
    # log the estimated RV uncertainty
    wargs = [fiber, wmeanref]
    wmsg = "On fiber {0} estimated RV uncertainty on spectrum is {1:.3f}"
    wlog(params, 'info', wmsg.format(*wargs))
    # ----------------------------------------------------------------------
    # Do the CCF calculations
    # ----------------------------------------------------------------------
    # get the mask parameters
    mkwargs = dict(filename=ccfmask, mask_width=mask_width,
                   mask_units=mask_units)
    ll_mask_d, ll_mask_ctr, w_mask = get_ccf_mask(params, **mkwargs)
    # calculate the CCF
    props = ccf_calculation(params, image, blaze, wavemap, berv, targetrv,
                            ccfwidth, ccfstep, ll_mask_ctr, w_mask, fit_type,
                            fiber)
    # ----------------------------------------------------------------------
    # Calculate the mean CCF
    # ----------------------------------------------------------------------
    # get the average ccf
    mean_ccf = np.nanmean(props['CCF'][: ccfnmax], axis=0)

    # get the fit for the normalized average ccf
    mean_ccf_coeffs, mean_ccf_fit = fit_ccf(params, 'mean', props['RV_CCF'],
                                            mean_ccf, fit_type=fit_type)
    # get the max cpp
    # TODO: How do we calculate max_cpp and what is it? Do we need it?
    # max_cpp = mp.nansum(props['CCF_MAX']) / mp.nansum(props['PIX_PASSED_ALL'])
    # get the RV value from the normalised average ccf fit center location
    ccf_rv = float(mean_ccf_coeffs[1])
    # get the contrast (ccf fit amplitude)
    ccf_contrast = np.abs(100 * mean_ccf_coeffs[0])
    # get the FWHM value
    ccf_fwhm = mean_ccf_coeffs[2] * fwhm()
    # ----------------------------------------------------------------------
    #  combined CCF_NOISE uncertainty
    rv_noise = 1.0 / np.sqrt(np.nansum(1.0 / props['CCF_NOISE'] ** 2))
    # ----------------------------------------------------------------------
    # log the stats
    wargs = [ccf_contrast, float(mean_ccf_coeffs[1]), rv_noise, ccf_fwhm]
    wmsg = ("Correlation: C={0:1f}[%] RV={1:.5f}[km/s] RV_NOISE={2:.5f}[km/s] "
            "FWHM={3:.4f}[km/s]")
    wlog(params, 'info', wmsg.format(*wargs))
    # ----------------------------------------------------------------------
    # add to output array
    props['TOT_SPEC_RMS'] = wmeanref
    props['ORD_SPEC_RMS'] = wmeanrefo
    props['MEAN_CCF'] = mean_ccf
    props['MEAN_RV'] = ccf_rv
    props['MEAN_CONTRAST'] = ccf_contrast
    props['MEAN_FWHM'] = ccf_fwhm
    props['MEAN_CCF_COEFFS'] = mean_ccf_coeffs
    props['MEAN_CCF_FIT'] = mean_ccf_fit
    props['MEAN_RV_NOISE'] = rv_noise
    # add constants to props
    props['CCF_MASK'] = ccfmask
    props['CCF_STEP'] = ccfstep
    props['CCF_WIDTH'] = ccfwidth
    props['TARGET_RV'] = targetrv
    props['CCF_SIGDET'] = noise_sigdet
    props['CCF_BOXSIZE'] = noise_size
    props['CCF_MAXFLUX'] = noise_thres
    props['CCF_NMAX'] = ccfnmax
    props['MASK_MIN'] = mask_min
    props['MASK_WIDTH'] = mask_width
    props['MASK_UNITS'] = mask_units

    # ----------------------------------------------------------------------
    # rv ccf plot
    # ----------------------------------------------------------------------
    if params['PLOT_ALL']:
        # loop around every order
        plot_ccf_rv_fit(params=params, x=props['RV_CCF'],
                        y=props['CCF'], yfit=props['CCF_FIT'],
                        kind='FP fiber={0}'.format(fiber),
                        rv=props['CCF_FIT_COEFFS'][:, 1], ccfmask=ccfmask,
                        orders=np.arange(len(props['CCF'])), order=None)
    if params['PLOT']:
        # the mean ccf
        plot_ccf_rv_fit(params=params, x=props['RV_CCF'],
                        y=mean_ccf, yfit=mean_ccf_fit,
                        kind='MEAN FP fiber={0}'.format(fiber),
                        rv=props['MEAN_CCF_COEFFS'][1], ccfmask=ccfmask,
                        orders=None, order=None)
    # return the rv props
    return props


def delta_v_rms_2d(spe, wave, sigdet, threshold, size):
    """
    Compute the photon noise uncertainty for all orders (for the 2D image)

    :param spe: numpy array (2D), the extracted spectrum
                size = (number of orders by number of columns (x-axis))
    :param wave: numpy array (2D), the wave solution for each pixel
    :param sigdet: float, the read noise (sigdet) for calculating the
                   noise array
    :param threshold: float, upper limit for pixel values, above this limit
                      pixels are regarded as saturated
    :param size: int, size (in pixels) around saturated pixels to also regard
                 as bad pixels

    :return dvrms2: numpy array (1D), the photon noise for each pixel (squared)
    :return weightedmean: float, weighted mean photon noise across all orders
    """
    # flag (saturated) fluxes above threshold as "bad pixels"
    with warnings.catch_warnings(record=True) as _:
        flag = spe < threshold
    # flag all fluxes around "bad pixels" (inside +/- size of the bad pixel)
    for i_it in range(1, 2 * size, 1):
        flag[:, size:-size] *= flag[:, i_it: i_it - 2 * size]
    # get the wavelength normalised to the wavelength spacing
    nwave = wave[:, 1:-1] / (wave[:, 2:] - wave[:, :-2])
    # get the flux + noise array
    sxn = (spe[:, 1:-1] + sigdet ** 2)
    # get the flux difference normalised to the flux + noise
    nspe = (spe[:, 2:] - spe[:, :-2]) / sxn
    # get the mask value
    maskv = flag[:, 2:] * flag[:, 1:-1] * flag[:, :-2]
    # get the total per order
    tot = np.nansum(sxn * ((nwave * nspe) ** 2) * maskv, axis=1)
    # convert to dvrms2
    with warnings.catch_warnings(record=True) as _:
        dvrms2 = (speed_of_light_ms ** 2) / abs(tot)
    # weighted mean of dvrms2 values
    weightedmean = 1. / np.sqrt(np.nansum(1.0 / dvrms2))
    # per order value
    weightedmeanorder = np.sqrt(dvrms2)
    # return dv rms and weighted mean
    return dvrms2, weightedmean, weightedmeanorder


def ccf_calculation(params, image, blaze, wavemap, berv, targetrv, ccfwidth,
                    ccfstep, mask_centers, mask_weights, fit_type, fiber,
                    **kwargs):
    # get properties from params
    blaze_norm_percentile = params['CCF_BLAZE_NORM_PERCENTILE']
    blaze_threshold = params['BLAZE_THRES']
    # get rvmin and rvmax
    rvmin = targetrv - ccfwidth
    if 'rvmin' in kwargs:
        rvmin = kwargs['rvmin']
    else:
        rvmin = params.get('RVMIN', rvmin)

    rvmax = targetrv + ccfwidth + ccfstep
    if 'rvmax' in kwargs:
        rvmax = kwargs['rvmax']
    else:
        rvmax = params.get('RVMAX', rvmax)

    # get the dimensions
    nbo, nbpix = image.shape
    # create a rv ccf range
    rv_ccf = np.arange(rvmin, rvmax, ccfstep)
    # storage of the ccf
    ccf_all = []
    ccf_noise_all = []
    ccf_all_fit = []
    ccf_all_results = []
    ccf_lines = []
    ccf_all_snr = []
    ccf_norm_all = []

    # switch normalization across all weights
    if params['MASK_NORM'] == 'all':
        mask_weights = mask_weights / np.nanmean(np.abs(mask_weights))

    # ----------------------------------------------------------------------
    # loop around the orders
    for order_num in range(nbo):
        # log the process
        wmsg = "Processing CCF fiber {0} for Order {1}"
        wlog(params, '', wmsg.format(fiber, order_num))
        # ------------------------------------------------------------------
        # get this orders values
        wa_ord = np.array(wavemap[order_num])
        sp_ord = np.array(image[order_num])
        bl_ord = np.array(blaze[order_num])

        # we express sp_ord as a flux in photons per km/s
        grad = speed_of_light * np.gradient(wa_ord) / wa_ord
        sp_ord = sp_ord / grad

        # normalize per-ord blaze to its peak value
        # this gets rid of the calibration lamp SED
        bl_ord /= np.nanpercentile(bl_ord, blaze_norm_percentile)
        # change NaNs in blaze to zeros
        bl_ord[~np.isfinite(bl_ord)] = 0.0
        # mask on the blaze
        with warnings.catch_warnings(record=True) as _:
            blazemask = bl_ord > blaze_threshold
        # get order mask centers and mask weights
        min_ord_wav = np.nanmin(wa_ord[blazemask])
        max_ord_wav = np.nanmax(wa_ord[blazemask])
        # adjust for rv shifts
        # min_ord_wav = min_ord_wav * (1 - rvmin / speed_of_light)
        # max_ord_wav = max_ord_wav * (1 - rvmax / speed_of_light)
        # mask the ccf mask by the order length
        mask_wave_mask = (mask_centers > min_ord_wav)
        mask_wave_mask &= (mask_centers < max_ord_wav)
        omask_centers = mask_centers[mask_wave_mask]
        omask_weights = mask_weights[mask_wave_mask]

        # ------------------------------------------------------------------
        # find any places in spectrum or blaze where pixel is NaN
        nanmask = np.isnan(sp_ord) | np.isnan(bl_ord)
        # ------------------------------------------------------------------
        # deal with no valid lines
        if np.sum(mask_wave_mask) == 0:
            # log all NaN
            wargs = [order_num, min_ord_wav, max_ord_wav]
            wmsg = ("Order {0} CCF mask for this order has no lines. CCF set "
                    "to NaNs. \n\t order wavelength: {1:.3f} to {2:.3f}")
            wlog(params, 'warning', wmsg.format(*wargs))
            # set all values to NaN
            ccf_all.append(np.repeat(np.nan, len(rv_ccf)))
            ccf_all_fit.append(np.repeat(np.nan, len(rv_ccf)))
            ccf_all_results.append(np.repeat(np.nan, 4))
            ccf_noise_all.append(np.nan)
            ccf_lines.append(0)
            ccf_all_snr.append(np.nan)
            ccf_norm_all.append(np.nan)
            continue
        # ------------------------------------------------------------------
        # deal with all nan
        if np.sum(nanmask) == nbpix:
            # log all NaN
            wargs = [order_num]
            wmsg = "Order {0} has no finite values. CCF set to NaNs."
            wlog(params, 'warning', wmsg.format(*wargs))
            # set all values to NaN
            ccf_all.append(np.repeat(np.nan, len(rv_ccf)))
            ccf_all_fit.append(np.repeat(np.nan, len(rv_ccf)))
            ccf_all_results.append(np.repeat(np.nan, 4))
            ccf_noise_all.append(np.nan)
            ccf_lines.append(0)
            ccf_all_snr.append(np.nan)
            ccf_norm_all.append(np.nan)
            continue
        # ------------------------------------------------------------------
        # set the spectrum or blaze NaN pixels to zero (dealt with by divide)
        sp_ord[nanmask] = 0
        bl_ord[nanmask] = 0
        # now every value that is zero is masked (we don't want to spline these)
        good = (sp_ord != 0) & (bl_ord != 0)
        weight_ord = np.array(good, dtype=float)
        # ------------------------------------------------------------------
        # spline the spectrum and the blaze
        spline_sp = iuv_spline(wa_ord[good], sp_ord[good], k=5, ext=1)
        spline_bl = iuv_spline(wa_ord[good], bl_ord[good], k=5, ext=1)
        spline_weight = iuv_spline(wa_ord, weight_ord, k=1, ext=1)
        # ------------------------------------------------------------------
        # set up the ccf for this order
        ccf_ord = np.zeros_like(rv_ccf)
        # ------------------------------------------------------------------
        # get the wavelength shift (dv) in relativistic way
        wave_shifts = relativistic_waveshift(rv_ccf - berv)
        # ------------------------------------------------------------------
        # propagating the extreme wave shifts to see if any lines fall off
        #  the domain that is considered valid for the spline

        # find the wave grid for the first shift
        wave_tmp_start = omask_centers * wave_shifts[0]
        # find the wave grid for the last shift
        wave_tmp_end = omask_centers * wave_shifts[-1]
        # find the valid lines within these limits
        # (ext=1 puts 0 when point is beyond domain)
        valid_lines_start = spline_bl(wave_tmp_start) != 0
        valid_lines_end = spline_bl(wave_tmp_end) != 0
        # combine the valid masks for start and end
        keep = valid_lines_start & valid_lines_end
        # ------------------------------------------------------------------
        # deal with no valid lines
        if np.sum(keep) == 0:
            # log all NaN
            wargs = [order_num]
            wmsg = ("Order {0} CCF mask for this order only has lines in "
                    "invalid blaze locations. CCF set to NaNs")
            wlog(params, 'warning', wmsg.format(*wargs))
            # set all values to NaN
            ccf_all.append(np.repeat(np.nan, len(rv_ccf)))
            ccf_all_fit.append(np.repeat(np.nan, len(rv_ccf)))
            ccf_all_results.append(np.repeat(np.nan, 4))
            ccf_noise_all.append(np.nan)
            ccf_lines.append(0)
            ccf_all_snr.append(np.nan)
            ccf_norm_all.append(np.nan)
            continue
        # ------------------------------------------------------------------
        # apply masks to centers and weights
        omask_centers = omask_centers[keep]
        omask_weights = omask_weights[keep]
        # normalise omask weights by
        if params['MASK_NORM'] == 'order':
            omask_weights = omask_weights / np.nanmean(omask_weights)

        # Number of photons at line centers for 1 CCF step
        sweights = spline_weight(omask_centers)
        nphot = spline_sp(omask_centers) * sweights / ccfstep

        # Poisson noise is a bit bigger because of weights
        wsum = np.sum(nphot * omask_weights)
        wsum2 = np.sum(nphot * omask_weights ** 2)
        # we can't calculate wnoise for negative values --> set to inf
        if (wsum <= 0) or (wsum2 <= 0):
            wargs = [order_num]
            wmsg = ("Order {0} nphot negative (poor spectrum) cannot "
                    "calculate CCF noise – setting noise and SNR to NaNs")
            wlog(params, 'warning', wmsg.format(*wargs))
            wsum, wnoise = 0.0, np.inf
        else:
            wnoise = np.sqrt(wsum2)
        # ------------------------------------------------------------------
        # set number of valid lines used to zero
        numlines = 0
        # loop around the rvs and calculate the CCF at this point
        part3 = spline_bl(omask_centers)
        for rv_element in range(len(rv_ccf)):
            wave_tmp = omask_centers * wave_shifts[rv_element]
            part1 = spline_sp(wave_tmp)
            part2 = spline_bl(wave_tmp)
            part4 = spline_weight(wave_tmp)
            numlines = np.sum(spline_bl(wave_tmp) != 0)
            # CCF is the division of the sums
            with warnings.catch_warnings(record=True) as _:
                ccf_element = ((part1 * part3) / part2) * omask_weights * part4
                ccf_ord[rv_element] = np.nansum(ccf_element)
        # ------------------------------------------------------------------
        # deal with NaNs in ccf
        if np.sum(np.isnan(ccf_ord)) > 0:
            # log all NaN
            wargs = [order_num]
            wmsg = ("Order {0} CCF generated a non-finite value. CCF set "
                    "to NaNs.")
            wlog(params, 'warning', wmsg.format(*wargs))
            # set all values to NaN
            ccf_all.append(np.repeat(np.nan, len(rv_ccf)))
            ccf_all_fit.append(np.repeat(np.nan, len(rv_ccf)))
            ccf_all_results.append(np.repeat(np.nan, 4))
            ccf_noise_all.append(np.nan)
            ccf_lines.append(0)
            ccf_all_snr.append(np.nan)
            ccf_norm_all.append(np.nan)
            continue
        # ------------------------------------------------------------------
        # TODO -- check that its ok to remove the normalization
        # TODO -- this should preserve the stellar flux weighting
        # normalise each orders CCF to median
        # TODO -- keep track of the norm factor, write a look-up table
        # TODO -- with reasonable mid-M values and use these values for
        # TODO -- all stars. At some point, have a temperature-dependent
        # TODO -- LUT of weights.
        ccf_norm = np.nanmedian(ccf_ord)
        # ccf_ord = ccf_ord / ccf_norm
        # ------------------------------------------------------------------
        # fit the CCF with a gaussian
        fargs = [order_num, rv_ccf, ccf_ord, fit_type]
        ccf_coeffs_ord, ccf_fit_ord = fit_ccf(params, *fargs)
        # ------------------------------------------------------------------
        # get the RV accuracy from Bouchy 2001 equation
        dv_pix = (np.gradient(ccf_ord) / np.gradient(rv_ccf)) / wnoise
        # set the bad values for ccf noise and ccf snr --> NaN value is bad
        if wsum == 0:
            ccf_noise = np.nan
            ccf_snr = np.nan
        else:
            ccf_noise = 1 / np.sqrt(np.nansum(dv_pix ** 2))
            # ge the snr
            ccf_snr = wsum / wnoise
        # ------------------------------------------------------------------
        # append ccf to storage
        ccf_all.append(ccf_ord)
        ccf_all_fit.append(ccf_fit_ord)
        ccf_all_results.append(ccf_coeffs_ord)
        ccf_noise_all.append(ccf_noise)
        ccf_lines.append(numlines)
        ccf_all_snr.append(ccf_snr)
        ccf_norm_all.append(ccf_norm)

    # store outputs in param dict
    props = dict()
    props['RV_CCF'] = rv_ccf
    props['CCF'] = np.array(ccf_all)
    props['CCF_LINES'] = np.array(ccf_lines)
    props['TOT_LINE'] = np.sum(ccf_lines)
    props['CCF_NOISE'] = np.array(ccf_noise_all) * 1000  # [m/s]
    props['CCF_SNR'] = np.array(ccf_all_snr)
    props['CCF_FIT'] = np.array(ccf_all_fit)
    props['CCF_FIT_COEFFS'] = np.array(ccf_all_results)
    props['CCF_NORM'] = np.array(ccf_norm_all)
    # return props
    return props


def fit_ccf(params, order_num, rv, ccf, fit_type):
    """
    Fit the CCF to a guassian function

    :param params:
    :param order_num:
    :param rv: numpy array (1D), the radial velocities for the line
    :param ccf: numpy array (1D), the CCF values for the line
    :param fit_type: int, if "0" then we have an absorption line
                          if "1" then we have an emission line

    :return result: numpy array (1D), the fit parameters in the
                    following order:

                [amplitude, center, fwhm, offset from 0 (in y-direction)]

    :return ccf_fit: numpy array (1D), the fit values, i.e. the gaussian values
                     for the fit parameters in "result"
    """
    func_name = __NAME__ + '.fit_ccf()'
    # deal with inconsistent lengths
    if len(rv) != len(ccf):
        eargs = [len(rv), len(ccf), func_name]
        emsg = ("Length of rv vector (len={0}) and ccf vector (len={1}) are "
                "not the same. \n\t Function = {2}")
        wlog(params, 'error', emsg.format(*eargs))

    # deal with all nans
    if np.sum(np.isnan(ccf)) == len(ccf):
        # log warning about all NaN ccf
        wargs = [order_num]
        wmsg = ("All NaN slice found when fitting CCF (Order = {0}). "
                "Returning NaNs instead of fit.")
        wlog(params, 'warning', wmsg.format(*wargs))
        # return NaNs
        result = np.zeros(4) * np.nan
        ccf_fit = np.zeros_like(ccf) * np.nan
        return result, ccf_fit

    # get constants
    max_ccf, min_ccf = np.nanmax(ccf), np.nanmin(ccf)
    argmin, argmax = np.nanargmin(ccf), np.nanargmax(ccf)
    diff = max_ccf - min_ccf
    rvdiff = rv[1] - rv[0]
    # set up guess for gaussian fit
    # if fit_type == 0 then we have absorption lines
    if fit_type == 0:
        if np.nanmax(ccf) != 0:
            a = np.array([-diff / max_ccf, rv[argmin], 4 * rvdiff, 0])
        else:
            a = np.zeros(4)
    # else (fit_type == 1) then we have emission lines
    else:
        a = np.array([diff / max_ccf, rv[argmax], 4 * rvdiff, 1])
    # normalise y
    y = ccf / max_ccf - 1 + fit_type
    # x is just the RVs
    x = rv
    # uniform weights
    w = np.ones(len(ccf))
    # get gaussian fit
    nanmask = np.isfinite(y)
    y[~nanmask] = 0.0
    # fit the gaussian
    try:
        with warnings.catch_warnings(record=True) as _:
            result, fit = fitgaussian(x, y, weights=w, guess=a)
    except RuntimeError:
        result = np.repeat(np.nan, 4)
        fit = np.repeat(np.nan, len(x))

    # scale the ccf
    ccf_fit = (fit + 1 - fit_type) * max_ccf

    # return the best guess and the gaussian fit
    return result, ccf_fit


# =============================================================================
# Plot functions
# =============================================================================
def plot_ccf_rv_fit(**kwargs):
    # ------------------------------------------------------------------
    # get the arguments from kwargs
    x = kwargs['x']
    y = kwargs['y']
    yfit = kwargs['yfit']
    kind = kwargs['kind']
    found_rv = kwargs['rv']
    ccf_mask = kwargs['ccfmask']
    order = kwargs.get('order', None)
    orders = kwargs.get('orders', None)
    # ------------------------------------------------------------------
    # if orders is None we assume x, y and yfit are just one vector
    if orders is None and order is None:
        # set y, yfit to length-1 lists of themselves
        y, yfit = [y], [yfit]
        found_rv = [found_rv]
        # order gen is just set to a list containing the first
        #     (and only) element
        order_gen = [0]
    # else if orders is None we get all orders
    elif order is None:
        order_gen = orders
    # else we just deal with the order specified
    else:
        order_gen = [order]
    # ------------------------------------------------------------------
    # deal with plot style
    black = 'black'
    # ------------------------------------------------------------------
    # loop around orders
    for order_num in order_gen:
        # get this orders values
        y_i = y[order_num]
        yfit_i = yfit[order_num]
        # work out the residuals
        res = y_i - yfit_i
        # ------------------------------------------------------------------
        # set up plot
        gs = dict(height_ratios=[2, 1])
        # noinspection PyTypeChecker
        fig, frames = plt.subplots(nrows=2, ncols=1, sharex=True,
                                   gridspec_kw=gs)
        # ------------------------------------------------------------------
        # plot x vs y and yfit
        frames[0].plot(x, y_i, label='data', marker='x', ls='None',
                       color=black)
        frames[0].plot(x, yfit_i, label='fit')
        # plot residuals
        frames[1].plot(x, res, label='residuals')
        # plot legends
        frames[0].legend(loc=0)
        frames[1].legend(loc=0)
        # set labels and title
        targs = ['({0})'.format(kind), found_rv[order_num], ccf_mask]
        title = 'CCF plot {0}\n RV={1:.5f} km/s Mask={2}'.format(*targs)
        if orders is not None:
            title = 'RV Fit plot. Order {0}'.format(order_num)
        frames[0].set(ylabel='CCF', title=title)
        frames[1].set(xlabel='RV [km/s]', ylabel='CCF')
        # ------------------------------------------------------------------
        # adjust size
        fig.subplots_adjust(hspace=0, wspace=0)
        # ------------------------------------------------------------------
        plt.show()
        plt.close()


def plot_ccf_swave_ref(**kwargs):
    # ------------------------------------------------------------------
    # get the arguments from kwargs
    wavemap = kwargs['wavemap']
    image = kwargs['image']
    fiber = kwargs['fiber']
    nbo = kwargs['nbo']
    # optional arguments
    order = kwargs.get('order', None)
    orders = kwargs.get('orders', None)
    # ------------------------------------------------------------------
    if order is None and orders is None:
        order_gen = np.arange(nbo)
    # else we check whether orders is set
    elif orders is not None:
        order_gen = list(orders)
    # else we just deal with the order specified
    elif order is not None:
        order_gen = [order]
    else:
        order_gen = [0]
    # ------------------------------------------------------------------
    # loop around orders
    for order_num in order_gen:
        # ------------------------------------------------------------------
        # set up plot
        fig, frame = plt.subplots(nrows=1, ncols=1)
        # ------------------------------------------------------------------
        # plot fits
        frame.plot(wavemap[order_num], image[order_num])
        # set title labels limits
        title = 'spectral order {0} fiber {1}'
        frame.set(xlabel='Wavelength [nm]', ylabel='flux',
                  title=title.format(order_num, fiber))
        # ------------------------------------------------------------------
        plt.show()
        plt.close()


def plot_ccf_photon_uncert(**kwargs):
    # ------------------------------------------------------------------
    # get the arguments from kwargs
    x = kwargs.get('x')
    y_sp = kwargs.get('y_sp')
    y_ccf = kwargs.get('y_cc')
    # get max/min points
    with warnings.catch_warnings(record=True) as _:
        ymin = np.nanmin(y_ccf)
        ymax = np.nanmax(y_ccf)
        if not np.isfinite(ymin):
            ymin = np.nanmin(y_sp)
        if not np.isfinite(ymax):
            ymax = np.nanmax(y_sp)
    # ------------------------------------------------------------------
    # set up plot
    fig, frame = plt.subplots(nrows=1, ncols=1)
    # ------------------------------------------------------------------
    # plot fits
    frame.plot(x, y_sp, label='DVRMS Spectrum', marker='x', linestyle='None')
    # plot ccf noise (unless all NaNs)
    if np.sum(np.isnan(y_ccf)) != len(y_ccf):
        frame.plot(x, y_ccf, label='DVRMS CCF', marker='o', linestyle='None')
    # set title labels limits
    title = 'Photon noise uncertainty versus spectral order'
    frame.set(xlabel='Order number', ylabel='Photon noise uncertainty [m/s]',
              title=title)
    # deal with limits (may be NaN)
    if np.isfinite(ymin) and np.isfinite(ymax):
        frame.set_ylim(bottom=ymin, top=ymax)
    elif np.isfinite(ymin):
        frame.set_ylim(bottom=ymin)
    else:
        frame.set_ylim(top=ymax)

    frame.legend(loc=0)
    # ------------------------------------------------------------------
    plt.show()
    plt.close()


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    print('Hello World')

# =============================================================================
# End of code
# =============================================================================
