"""
Created on Oct. 14 2020

Author: Thomas Vandal <thomas.vandal@umontreal.ca>

The spirou_ccf_analysis script by Eder Martioli was used as a starting point
for this script (https://github.com/edermartioli/spirou-ccf-analysis).

The tools used under the hood to calculate RV timeseries from CCFs were written
by Étienne Artigau and edited by Eder Martioli and Thomas Vandal.
"""
import os
import glob

import numpy as np

from ccf2rv import get_object_rv
import getrv_utils as gu


# CLI arguments
clargs = gu.parse_clargs()

cfg = gu.get_config(clargs)

# If no pattern, build it from other information
# In both cases, we load a list of CCF files
if cfg['pattern'] is None:
    ccf_files = gu.get_ccf_files(
            cfg['path'],
            cfg['object'],
            cfg['mask'],
            cfg['sanitize']
            )

else:
    if cfg['verbose']:
        print('WARNING: pattern overrides mask, object, and sanitize')
    pattern = os.path.join(cfg['path'], cfg['pattern'])
    ccf_files = sorted(glob.glob(pattern))

# Get dictionnary where each key is the "ID" of a timeseries with given
# object, mask, sanitize, DRS version. The keys have format
# object__ mask__sani__version.
ccf_dict = gu.get_ccf_dict(ccf_files)

bandpass = ''.join(np.unique(cfg['bandpass']))

for run_id in ccf_dict.keys():

    list_of_files = ccf_dict[run_id]

    if cfg['verbose']:
        print('Processing batch {0} containing {1} files'
              .format(run_id, len(list_of_files))
              )

    tbl = get_object_rv(
            ccf_files=list_of_files,
            outdir=cfg['outdir'],
            run_id=run_id,
            method=cfg['method'],
            exclude_orders=cfg['exclude_orders'],
            force=cfg['force'],
            snr_min=cfg['snr_min'],
            bandpass=bandpass,
            velocity_window=cfg['velocity_window'],
            dvmax_per_order=cfg['dvmax_per_order'],
            save_ccf_cube=cfg['save_cube'],
            save_result_table=cfg['save_table'],
            save_rv_timeseries=cfg['save_rv'],
            save_weight_table=cfg['save_weights'],
            bin_rv_timeseries=cfg['bin'],
            median_bin=cfg['median'],
            doplot=cfg['do_plots'],
            showplots=cfg['show_plots'],
            saveplots=cfg['save_plots'],
            )
