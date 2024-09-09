#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
import glob
import os
from typing import Any, Dict

import numpy as np
from astropy.io import fits
from tqdm import tqdm
from scipy.interpolate import InterpolatedUnivariateSpline

from apero_checks.core import apero_functions
from apero_checks.core import misc


def estimate_sigma(sp):
    # returns the standard deviation of a spectrum
    n1, p1 = np.nanpercentile(sp, [16, 84])
    return (p1 - n1) / 2


def lowpassfilter(input_vect: np.ndarray, width: int = 101,
                  k=1, frac_valid_min=0) -> np.ndarray:
    """
    Computes a low-pass filter of an input vector.

    This is done while properly handling NaN values, but at the same time
    being reasonably fast.

    Algorithm:

    provide an input vector of an arbitrary length and compute a running NaN
    median over a box of a given length (width value). The running median is
    NOT computed at every pixel but at steps of 1/4th of the width value.
    This provides a vector of points where the nan-median has been computed
    (ymed) and mean position along the input vector (xmed) of valid (non-NaN)
    pixels. This xmed/ymed combination is then used in a spline to recover a
    vector for all pixel positions within the input vector.

    When there are no valid pixel in a 'width' domain, the value is skipped
    in the creation of xmed and ymed, and the domain is splined over.

    :param input_vect: numpy 1D vector, vector to low pass
    :param width: int, width (box size) of the low pass filter
    :param k: int, order of the spline interpolation
    :param frac_valid_min: float, minimum fraction of valid pixels in a
                           'width' domain to compute the low pass filter.
                           If the fraction of valid pixels is below this value,
                           the low pass filter is not computed and the value
                           is interpolated over.

    :return: np.array, the low-pass of the input_vector
    """
    # set function name
    # _ = display_func('lowpassfilter', __NAME__)
    # indices along input vector
    index = np.arange(len(input_vect))
    # placeholders for x and y position along vector
    xmed = []
    ymed = []
    # loop through the lenght of the input vector
    for it in np.arange(-width // 2, len(input_vect) + width // 2, width // 4):
        # if we are at the start or end of vector, we go 'off the edge' and
        # define a box that goes beyond it. It will lead to an effectively
        # smaller 'width' value, but will provide a consistent result at edges.
        low_bound = it
        high_bound = it + int(width)
        # deal with lower bounds out of bounds --> set to zero
        if low_bound < 0:
            low_bound = 0
        # deal with upper bounds out of bounds --> set to max
        if high_bound > (len(input_vect) - 1):
            high_bound = (len(input_vect) - 1)
        # get the pixel bounds
        pixval = index[low_bound:high_bound]
        # do not low pass if not enough points
        if len(pixval) < 3:
            continue
        # if no finite value, skip
        if np.mean(np.isfinite(input_vect[pixval])) <= frac_valid_min:
            continue
        # mean position along vector and NaN median value of
        # points at those positions
        xmed.append(np.nanmean(pixval))
        ymed.append(np.nanmedian(input_vect[pixval]))
    # convert to arrays
    xmed = np.array(xmed, dtype=float)
    ymed = np.array(ymed, dtype=float)
    # we need at least 3 valid points to return a
    # low-passed vector.
    if len(xmed) < 3:
        return np.zeros_like(input_vect) + np.nan
    # low pass with a mean
    if len(xmed) != len(np.unique(xmed)):
        xmed2 = np.unique(xmed)
        ymed2 = np.zeros_like(xmed2)
        for i in range(len(xmed2)):
            ymed2[i] = np.mean(ymed[xmed == xmed2[i]])
        xmed = xmed2
        ymed = ymed2
    # splining the vector
    spline = InterpolatedUnivariateSpline(xmed, ymed, k=k, ext=3)
    lowpass = spline(np.arange(len(input_vect)))
    # return the low pass filtered input vector
    return lowpass


# =============================================================================
# Define variables
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)

# =============================================================================
# Define functions
# =============================================================================
def test(params: Dict[str, Any], obsdir: str, log=False) -> bool:
    """
    Test for excess modal noise in telluric stars. It checks for tcorr files
    of vetted telluric stars and computes the rms of the pixel to pixel vs
    the rms with a stride a 20 pixels on a sample order in H band. The pixel
    to pixel rms should be smaller than the rms with a stride of 20 pixels. The
    pixel-to-pixel samples photon noise while the other one samples the modal
    noise. The test is passed if the pixel-to-pixel rms is smaller than the
    modal noise rms by a known thershold.
    """

    # featureless telluric stars
    vetted_stars = ['15_PEG', 'HR1903', 'HR3117', 'HR3131', 'HR3314', 'HR4023',
                    'HR4467', 'HR4468', 'HR4889', 'HR5671', 'HR6743', 'HR7590',
                    'HR8709', 'ZETVIR']

    sample_order = 58  # order in the middle of H band, clean from tellurics
    threshold_HA = 0.015  # excess noise in HA mode
    threshold_HE = 0.010  # same in HE mode

    # update apero-profile
    apero_params = apero_functions.update_apero_profile(params)
    # tmp directory
    red_dir = apero_params['DRS_DATA_REDUC']
    # directory to check
    obsdir_path = os.path.join(red_dir, obsdir)

    # -------------------------------------------------------------------------
    if log:
        msg = 'Analysing observation directory: {0}'
        margs = [obsdir]
        misc.log_msg(msg.format(*margs), level='')
    # -------------------------------------------------------------------------

    # check if directory exists
    if not os.path.exists(obsdir_path):
        if log:
            print('red directory {} does not exist'.format(obsdir))
        return True

    # list of all the files in observation directory
    files = glob.glob(os.path.join(obsdir_path, '*e2ds*tcorr_A.fits'))

    # check if there are files in the directory
    if len(files) == 0:
        if log:
            print('No torr files in directory {}'.format(obsdir))
        return True

    passed = True
    failed_msg = ''

    no_telluric_stars = True
    # check pixel shift header keys for all files in obsdir
    for filename in tqdm(files, leave=False):
        hdr = fits.getheader(filename)
        if hdr['DRSOBJN'] not in vetted_stars:
            continue

        no_telluric_stars = False

        sp = fits.getdata(filename)[sample_order]
        # keep only the central 50% of the 4088 spectrum
        sp = sp[2048 - 1024:2048 + 1024]
        sp /= lowpassfilter(sp, 101, k=1)
        rms_pixel_to_pixel = estimate_sigma(sp - np.roll(sp, 1))
        # rms with a strid of 20 pixel
        rms_pixel_to_pixel_20 = estimate_sigma(sp - np.roll(sp, 20))
        # quadratic subtraction of the two
        if rms_pixel_to_pixel_20 > rms_pixel_to_pixel:
            rms_pixel_to_pixel_lf = np.sqrt(rms_pixel_to_pixel_20 ** 2 - rms_pixel_to_pixel ** 2)
        else:
            rms_pixel_to_pixel_lf = 0

        if hdr['DRSMODE'] == 'HA':
            threshold = threshold_HA
        elif hdr['DRSMODE'] == 'HE':
            threshold = threshold_HE
        else:
            # no threshold for other modes
            return True

        rmsmsg = '(rms20 = {:.4f}, rms1 = {:.4f}), rmslf = {:4,f})'
        margs = [rms_pixel_to_pixel_20, rms_pixel_to_pixel,
                 rms_pixel_to_pixel_lf]

        if rms_pixel_to_pixel_lf > threshold:
            passed = False
            if log:
                msg = ('Excess modal noise detected in {0} {1}\nfile={2}\n\n')
                fargs = [hdr['DRSOBJN'], rmsmsg.format(*margs), filename]
                failed_msg += msg.format(*fargs)
        else:
            msg = ('No excess noise detected in {0} {1}\nfile={2}\n\n')
            if log:
                fargs = [hdr['DRSOBJN'], rmsmsg.format(*margs), filename]
                failed_msg += msg.format(*fargs)

    if no_telluric_stars:
        if log:
            print('No telluric stars in directory {}'.format(obsdir))
        return True

    if log:
        if len(failed_msg) > 0:
            print(failed_msg)
        else:
            print('No excess noise detected in telluric stars')

    return passed


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # provide a working example of a test
    _params = dict()
    # define the observation directory
    _obsdir = '2021-03-15'
    # run the test
    test(_params, _obsdir, log=True)

# =============================================================================
# End of code
# =============================================================================
