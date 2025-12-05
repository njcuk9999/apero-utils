#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-11-24

@author: cook
"""
import glob
import numpy as np
from tqdm import tqdm


# =============================================================================
# Define variables
# =============================================================================
workspace = '/spirou2/drs-data/full_211018/msg/processing/APEROG-PID-00016365840209758570-CQXV_apero_processing_group'

recipe = 'apero_extract_spirou'

# =============================================================================
# Define functions
# =============================================================================


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    logfiles = glob.glob('*{0}.log'.format(recipe))

    files_with_errors = []

    for logfile in tqdm(logfiles):

        # read log file
        with open(logfile, 'r') as logf:
            lines = logf.readlines()

        for line in lines:
            if '!!|' in line:
                files_with_errors.append(logfile)
                break

    calib_error_files = []
    calib_obs_dirs = []
    target_names = []

    for logfile in files_with_errors:
        # read log file
        with open(logfile, 'r') as logf:
            lines = logf.readlines()

        # find OBS_DIR
        obs_dir = None
        for line in lines:
            if '--OBS_DIR:' in line:
                obs_dir = line.split('--OBS_DIR: ')[-1].replace('\n', '')

        # find target name
        target_name = None
        for line in lines:
            if 'Identified object as ' in line:
                target_name = line.split('Identified object as ')[-1].split(' with')[0]

        for line in lines:
            if '09-002-00004' in line:
                calib_error_files.append(logfile)
                calib_obs_dirs.append(obs_dir)
                target_names.append(target_name)
                break

    # counters
    count_obs_dir = dict()
    count_target_names = dict()

    for item in calib_obs_dirs:
        if item in count_obs_dir:
            count_obs_dir[item] += 1
        else:
            count_obs_dir[item] = 1

    for item in target_names:
        if item in count_target_names:
            count_target_names[item] += 1
        else:
            count_target_names[item] = 1

    for logfile in files_with_errors:
        # read log file
        with open(logfile, 'r') as logf:
            lines = logf.readlines()

        minpos = np.inf
        maxpos = -np.inf
        for it, line in enumerate(lines):
            if '!!|' in line:
                minpos = int(np.min([minpos, it]))
                maxpos = int(np.max([maxpos, it]))

        print()
        print('=' * 50)
        print('Log file: {0}'.format(logfile))
        print('=' * 50)

        printlines = lines[minpos - 5: maxpos + 5]
        for printline in printlines:
            print(printline.replace('\n', ''))

# =============================================================================
# End of code
# =============================================================================
