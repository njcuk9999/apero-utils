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

from apero_checks.core import apero_functions
from apero_checks.core import misc


def estimate_sigma(sp):
    # returns the standard deviation of a spectrum
    n1, p1 = np.nanpercentile(sp, [16, 84])
    return (p1 - n1) / 2


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
    threshold_HA = 0.012  # excess noise in HA mode
    threshold_HE = 0.007  # same in HE mode

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
            print('tmp directory {} does not exist'.format(obsdir))
        return False

    # list of all the files in observation directory
    files = glob.glob(os.path.join(obsdir_path, '*e2ds*tcorr_A.fits'))

    # check if there are files in the directory
    if len(files) == 0:
        if log:
            print('No torr files in directory {}'.format(obsdir))
        return False

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
        sp /= np.median(sp)
        rms_pixel_to_pixel = estimate_sigma(sp - np.roll(sp, 1))
        # rms with a strid of 20 pixel
        rms_pixel_to_pixel_20 = estimate_sigma(sp - np.roll(sp, 20))
        # quadratic subtraction of the two
        if rms_pixel_to_pixel_20 > rms_pixel_to_pixel:
            rms_pixel_to_pixel_20 = np.sqrt(rms_pixel_to_pixel_20 ** 2 - rms_pixel_to_pixel ** 2)
        else:
            rms_pixel_to_pixel_20 = 0

        if hdr['DRSMODE'] == 'HA':
            threshold = threshold_HA
        elif hdr['DRSMODE'] == 'HE':
            threshold = threshold_HE
        else:
            # no threshold for other modes
            return True

        if rms_pixel_to_pixel_20 > threshold:
            passed = False
            msg = 'Excess modal noise detected in {} (rms = {:.4f})\n'
            failed_msg += msg.format(hdr['DRSOBJN'], rms_pixel_to_pixel_20)
        else:
            msg = 'No excess noise detected in {} (rms = {:.4f})\n'
            if log:
                print(msg.format(hdr['DRSOBJN'], rms_pixel_to_pixel_20))

    if no_telluric_stars:
        if log:
            print('No telluric stars in directory {}'.format(obsdir))
        return False

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
