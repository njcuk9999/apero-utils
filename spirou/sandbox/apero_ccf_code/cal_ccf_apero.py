#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2019-09-19 at 13:16

@author: cook
"""
from astropy.io import fits
import numpy as np
import os

import ccf_functions as cf

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = 'cal_ccf_apero.py'
__INSTRUMENT__ = 'SPIROU'
# Get version and author
__version__ = '0.6.132-ccf1'
__author__ = 'Neil Cook'
__date__ = '2020-11-05'


# =============================================================================
# Define functions
# =============================================================================
# All recipe code goes in _main
#    Only change the following from here:
#     1) function calls  (i.e. main(arg1, arg2, **kwargs)
#     2) fkwargs         (i.e. fkwargs=dict(arg1=arg1, arg2=arg2, **kwargs)
#     3) config_main  outputs value   (i.e. None, pp, reduced)
# Everything else is controlled from recipe_definition
def main(config=None, **kwargs):
    """
    Main function for cal_ccf_spirou.py

    :param directory: string, the night name sub-directory
    :param files: list of strings or string, the list of files to process
    :param kwargs: any additional keywords

    :type directory: str
    :type files: list[str]

    :keyword debug: int, debug level (0 for None)

    :returns: dictionary of the local space
    :rtype: dict
    """
    # function
    # assign function calls (must add positional)
    fkwargs = dict(config=config, **kwargs)
    # get parameters
    params = cf.get_params(**fkwargs)
    recipe = None
    # title
    cf.wlog(params, '', '=' * 50)
    targs = [__NAME__, __version__, __date__]
    cf.wlog(params, '', '\t{0} ({1}, {2})'.format(*targs))
    cf.wlog(params, '', '=' * 50)
    # print inputs
    cf.wlog(params, 'info', 'Main input parameters:')
    pparams = ['CONFIG', 'WORKING_DIR', 'INFILE_AB', 'WAVEFILE_AB', 'BLAZEFILE',
               'CCF_MASK_AB', 'CCF_STEP', 'CCF_WIDTH', 'CCF_RV',
               'INFILE_C', 'WAVEFILE_C', 'CCF_MASK_C', 'WAVE_CCF_STEP',
               'WAVE_CCF_WIDTH', 'WAVE_CCF_TARGET_RV']
    for pparam in pparams:
        cf.wlog(params, 'info', '\t--{0}: {1}'.format(pparam.lower(),
                                                      params[pparam]))
    cf.wlog(params, '', '=' * 50)
    # ----------------------------------------------------------------------
    # run main bulk of code (catching all errors)
    llmain, success = __main__(recipe, params)
    # ----------------------------------------------------------------------
    # End Message
    # ----------------------------------------------------------------------
    return llmain


def __main__(recipe, params):
    """
    Main code: should only call recipe and params (defined from main)

    :param recipe:
    :param params:
    :return:
    """
    # ----------------------------------------------------------------------
    # Main Code
    # ----------------------------------------------------------------------
    mainname = __NAME__ + '.__main__()'
    # get files
    infile = os.path.join(params['WORKING_DIR'], params['INFILE_AB'])
    # get header from file instance
    header = fits.getheader(infile)
    # ------------------------------------------------------------------
    # check that file has valid DPRTYPE
    # ------------------------------------------------------------------
    dprtype = header['DPRTYPE']
    # if dprtype is incorrect skip
    if dprtype not in params['CCF_ALLOWED_DPRTYPES']:
        # join allowed dprtypes
        allowed_dprtypes = ', '.join(params.listp('CCF_ALLOWED_DPRTYPES'))
        # log that we are skipping
        wargs = [dprtype, recipe.name, allowed_dprtypes, infile.basename]
        wmsg = ('DPRTYPE = ‘{0}’ is not valid for {1}. \n\t Allowed DPRTYPES '
                'are: {2} \n\t Skipping filename = {3}')
        cf.wlog(params, 'error', wmsg.format(*wargs))
    # flag whether calibration fiber is FP
    has_fp = dprtype.upper().split('_')[1] == 'FP'
    # ------------------------------------------------------------------
    # get fiber from infile
    fiber = header['FIBER']
    # ----------------------------------------------------------------------
    # Check we are using correct fiber
    # ----------------------------------------------------------------------
    sfiber, rfiber = params['FIBER_CCF']
    if fiber != sfiber:
        # log that the science fiber was not correct
        eargs = [fiber, sfiber, infile.name, infile.filename]
        emsg = ("Fiber={0} is not valid for 'science'. (Required '{1}') "
                "\n\t File {2}: {3} \n\t Function = {4}'")
        cf.wlog(params, 'error', emsg.format(*eargs))

    # ------------------------------------------------------------------
    # Get barycentric corrections (BERV)
    # ------------------------------------------------------------------
    bprops = cf.get_berv(params, header)
    # ------------------------------------------------------------------
    # load wavelength solution for this fiber
    wprops = cf.get_wavesolution(params, header, fiber='AB')
    # ------------------------------------------------------------------
    # Get blaze
    # ------------------------------------------------------------------
    blazefile, blaze = cf.get_blaze(params, header)

    # ------------------------------------------------------------------
    #   Remove domain with telluric > 50%
    # ------------------------------------------------------------------
    outtype = header['DRSOUTID']

    if outtype in params['CCF_CORRECT_TELLU_TYPES']:
        # remove telluric domain below a defined threshold
        #    and return the infile (with infile.data updated)
        image = cf.remove_telluric_domain(params, infile)
    else:
        image = np.array(infile.data)

    # ------------------------------------------------------------------
    # Compute CCF on science channel
    # ------------------------------------------------------------------
    # log progress: Computing CCF on fiber
    cf.wlog(params, '', "Computing CCF on fiber {0}".format(fiber))
    # compute ccf
    cargs = [infile, image, header, blaze, wprops['WAVEMAP'], bprops, fiber]
    rv_props1 = cf.compute_ccf_science(params, *cargs)

    # ------------------------------------------------------------------
    # Compute CCF on reference fiber (FP only)
    # ------------------------------------------------------------------
    if has_fp:
        # find the c fiber file
        infile_r = cf.locate_reference_file(params)
        image_r, header_r = fits.getdata(infile_r, header=True)
        # get the wave solution associated with this file
        wprops_r = cf.get_wavesolution(params, header_r, fiber='C')
        # get c fiber file time
        filetime_r = header['MJDMID']
        # --------------------------------------------------------------
        # deal with differing wavelength solutions (between science and
        #    reference)
        if wprops['WAVETIME'] != wprops_r['WAVETIME']:
            # log warning
            wargs = [wprops_r['WAVETIME'], wprops['WAVETIME'],
                     wprops_r['WAVEFILE'], wprops['WAVEFILE']]
            wmsg = ("Reference wave solution time did not match science wave "
                    "solution time. Using science wave solution for reference "
                    "fiber. \n\t Reference wave time = {0}   science wave "
                    "time = {1} \n\t Reference wave file: {2} \n\t Science "
                    "wave file: {3}")
            cf.wlog(params, 'warning', wmsg.format(*wargs))
            # set the reference wave solution to the science wave solution
            wprops_r = wprops
        # log progress: Computing CCF on fiber
        cf.wlog(params, '', "Computing CCF on fiber {0}".format(fiber))
        # --------------------------------------------------------------
        # Compute CCF on reference channel
        cargs = [image_r, blaze, wprops_r['WAVEMAP'], rfiber]
        rv_props2 = cf.compute_ccf_fp(params, *cargs)
        # get the time difference (between file and wave)
        timediff = filetime_r - wprops_r['WAVETIME']
        # --------------------------------------------------------------
        # compute the rv output stats
        # --------------------------------------------------------------
        # need to deal with no drift from wave solution
        if wprops_r['WFP_DRIFT'] is None:
            rv_wave_fp = np.nan
            rv_simu_fp = rv_props2['MEAN_RV']
            rv_drift = rv_simu_fp
            rv_obj = rv_props1['MEAN_RV']
            rv_corrected = rv_obj - rv_drift
        # else we have drift from wave solution
        else:
            rv_wave_fp = wprops_r['WFP_DRIFT']
            rv_simu_fp = rv_props2['MEAN_RV']
            rv_drift = rv_wave_fp - rv_simu_fp
            rv_obj = rv_props1['MEAN_RV']
            rv_corrected = rv_obj - rv_drift
    # need to deal with no drift from wave solution and no simultaneous FP
    elif wprops['WFP_DRIFT'] is None:
        # set rv_props2
        rv_props2 = dict()
        # compute the stats
        rv_wave_fp = np.nan
        rv_simu_fp = np.nan
        rv_drift = 0.0
        rv_obj = rv_props1['MEAN_RV']
        rv_corrected = rv_obj - rv_drift
        infile_r = None
        header_r = None
        wprops_r = None
        timediff = np.nan
    # need way to deal no simultaneous FP
    else:
        # set rv_props2
        rv_props2 = dict()
        # compute the stats
        rv_wave_fp = 0.0
        rv_simu_fp = 0.0
        rv_drift = 0.0
        rv_obj = rv_props1['MEAN_RV']
        rv_corrected = rv_obj - rv_drift
        infile_r = None
        header_r = None
        wprops_r = None
        timediff = np.nan
    # ------------------------------------------------------------------
    # add rv stats to properties
    rv_props1['RV_WAVEFILE'] = wprops['WAVEFILE']
    rv_props1['RV_WAVETIME'] = wprops['WAVETIME']
    rv_props1['RV_WAVESRCE'] = wprops['WAVESOURCE']
    rv_props1['RV_TIMEDIFF'] = timediff
    rv_props1['RV_WAVE_FP'] = rv_wave_fp
    rv_props1['RV_SIMU_FP'] = rv_simu_fp
    rv_props1['RV_DRIFT'] = rv_drift
    rv_props1['RV_OBJ'] = rv_obj
    rv_props1['RV_CORR'] = rv_corrected
    # add the fp fiber properties
    if has_fp:
        rv_props2['RV_WAVEFILE'] = wprops_r['WAVEFILE']
        rv_props2['RV_WAVETIME'] = wprops_r['WAVETIME']
        rv_props2['RV_WAVESRCE'] = wprops_r['WAVESOURCE']
        rv_props2['RV_TIMEDIFF'] = timediff
        rv_props2['RV_WAVE_FP'] = rv_wave_fp
        rv_props2['RV_SIMU_FP'] = rv_simu_fp
        rv_props2['RV_DRIFT'] = rv_drift
        rv_props2['RV_OBJ'] = rv_obj
        rv_props2['RV_CORR'] = rv_corrected
    # ------------------------------------------------------------------
    # Quality control
    # ------------------------------------------------------------------
    # set passed variable and fail message list
    fail_msg, qc_values, qc_names, qc_logic, qc_pass = [], [], [], [], []
    # no quality control currently
    qc_values.append('None')
    qc_names.append('None')
    qc_logic.append('None')
    qc_pass.append(1)
    # ------------------------------------------------------------------
    # finally log the failed messages and set QC = 1 if we pass the
    # quality control QC = 0 if we fail quality control
    if np.sum(qc_pass) == len(qc_pass):
        wmsg = "QUALITY CONTROL SUCCESSFUL - Well Done -"
        cf.wlog(params, 'info', wmsg)
        passed = 1
    else:
        for farg in fail_msg:
            wmsg = "QUALITY CONTROL FAILED:" + farg
            cf.wlog(params, 'warning', wmsg)
        passed = 0
    # store in qc_params
    qc_params = [qc_names, qc_values, qc_logic, qc_pass]

    # ------------------------------------------------------------------
    # archive ccf from science fiber
    # ------------------------------------------------------------------
    cf.write_file(params, rv_props1, infile, params['CCF_MASK_AB'], header,
                  wprops)

    # ------------------------------------------------------------------
    # archive ccf from reference fiber
    # ------------------------------------------------------------------
    if has_fp:
        cf.write_file(params, rv_props2, infile_r, params['CCF_MASK_C'],
                      header_r, wprops_r)

    # ----------------------------------------------------------------------
    # End of main code
    # ----------------------------------------------------------------------
    cf.wlog(params, 'info', '=' * 50)
    msg = 'Recipe {0} has been successfully completed'
    cf.wlog(params, 'info', msg.format(__NAME__))
    cf.wlog(params, 'info', '=' * 50)
    # return
    return locals(), True


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # run main with no arguments (get from command line - sys.argv)
    # ll = main()

    ll = main(config='ccf.yaml',
              INFILE_AB='1234568o_pp_e2dsff_tcorr_AB.fits',
              INFILE_C='1234568o_pp_e2dsff_C.fits',
              )

# =============================================================================
# End of code
# =============================================================================
