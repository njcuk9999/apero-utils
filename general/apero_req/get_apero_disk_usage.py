#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-05-02 at 10:21

@author: cook
"""
import os
import socket
from typing import Tuple

import numpy as np
import yaml
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
EXCLUDE_DIRS = ['raw', 'lbl', 'objects']

PROFILES = dict()
PROFILES['spirou_offline'] = dict()
PROFILES['spirou_offline']['machine'] = 'titan.astro.umontreal.ca'
PROFILES['spirou_offline']['dir'] = '/cosmos99/spirou/apero-data/spirou_offline'


PROFILES['spirou_minidata'] = dict()
PROFILES['spirou_minidata']['machine'] = 'jupiter.astro.umontreal.ca'
PROFILES['spirou_minidata']['dir'] = '/scratch2/spirou/drs-data/spirou_minidata2_07286_jupiter'

PROFILES['spirou_xxs'] = dict()
PROFILES['spirou_xxs']['machine'] = 'jupiter.astro.umontreal.ca'
PROFILES['spirou_xxs']['dir'] = '/scratch2/spirou/drs-data/spirou_xxs'

PROFILES['nirps_he_online'] = dict()
PROFILES['nirps_he_online']['machine'] = 'rali.astro.umontreal.ca'
PROFILES['nirps_he_online']['dir'] = '/cosmos99/nirps/apero-data/nirps_he_online'

PROFILES['nirps_ha_online'] = dict()
PROFILES['nirps_ha_online']['machine'] = 'rali.astro.umontreal.ca'
PROFILES['nirps_ha_online']['dir'] = '/cosmos99/nirps/apero-data/nirps_ha_online'

PROFILES['nirps_he_xxs'] = dict()
PROFILES['nirps_he_xxs']['machine'] = 'jupiter.astro.umontreal.ca'
PROFILES['nirps_he_xxs']['dir'] = '/scratch2/nirps/drs-data/nirps_he_xxs'

PROFILES['nirps_ha_xxs'] = dict()
PROFILES['nirps_ha_xxs']['machine'] = 'jupiter.astro.umontreal.ca'
PROFILES['nirps_ha_xxs']['dir'] = '/scratch2/nirps/drs-data/nirps_ha_xxs'

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def get_total_size_gb(directory) -> Tuple[float, int]:


    files = []
    visited_inodes = set()
    print('\t\tFinding files')
    # get all files
    for dirpath, dirnames, filenames in os.walk(directory, followlinks=True):

        # Get and track the inode of the current directory
        try:
            dir_stat = os.stat(dirpath, follow_symlinks=True)
            if dir_stat.st_ino in visited_inodes:
                continue
            visited_inodes.add(dir_stat.st_ino)
        except OSError:
            continue  # skip inaccessible dirs

        for filename in filenames:
            files.append(os.path.join(dirpath, filename))
    print('\t\tFound {} files'.format(len(files)))
    print('\t\tCalculating size')
    total_size = 0

    for filename in tqdm(files, leave=False):
        try:
            stat = os.stat(filename, follow_symlinks=True)
            if stat.st_ino not in visited_inodes:
                visited_inodes.add(stat.st_ino)
                total_size += stat.st_size
        except OSError:
            # skip files that can't be accessed
            pass

    size_gb = total_size / (1024 ** 3)
    print('\t\tSize: {:.2f} GB'.format(size_gb))

    return size_gb, len(files)


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":


    machine = socket.gethostname()
    stats = dict()
    # ----------------------------------------------------------------------
    # loop around profiles
    for profile_name in PROFILES:
        # get profile settings
        profile = PROFILES[profile_name]
        # skip other machines
        if profile['machine'] != machine:
            continue
        # print profile name
        print('*' * 50)
        print(profile_name)
        print('*' * 50)
        # storage
        stats[profile_name] = dict()
        # get all directories in the profile directory
        listdir = os.listdir(profile['dir'])
        directories = dict()
        for dirname in listdir:
            directories[dirname] = os.path.join(profile['dir'], dirname)

        # ---------------------------------------------------------------------
        # Get stats on "the raw"
        # ---------------------------------------------------------------------
        if 'raw' in directories:
            # print progress
            print('Calculating size of raw files')
            # get the total size and number of files
            raw_size, raw_num = get_total_size_gb(directories['raw'])

            print(f'Raw size: {raw_size:.2f} Gb')
            print(f'Raw num: {raw_num}')
            # save raw size and raw num
            stats[profile_name]['raw_num'] = raw_num
            stats[profile_name]['raw_size'] = raw_size

        # ---------------------------------------------------------------------
        # Get stats on "the rest"
        # ---------------------------------------------------------------------
        rest_size = 0
        rest_num = 0
        stats[profile_name]['dir'] = dict()
        # loop around all directories
        for dirname in directories:
            # skip raw directory
            if dirname in EXCLUDE_DIRS:
                continue
            # print progress
            print(f'Calculating size of {dirname} files')
            # get the total size and number of files
            dir_size, dir_num = get_total_size_gb(directories[dirname])
            # add to stats
            stats[profile_name]['dir'][f'{dirname}_num'] = dir_size
            stats[profile_name]['dir'][f'{dirname}_size'] = dir_num
            # add to the total
            rest_size += dir_size
            rest_num += dir_num

        print(f'Rest size: {rest_size:.2f} Gb')
        print(f'rest num: {rest_num}')


        stats[profile_name]['rest_num'] = rest_num
        stats[profile_name]['rest_size'] = rest_size

    # make sure all size numbers are rounded to 2 decimal places
    # make sure all nums are integers
    for profile_name in stats:
        for key in stats[profile_name]:

            if key == 'dir':
                for key in stats[profile_name]['dir']:

                    value = stats[profile_name]['dir'][key]

                    if key.endswith('size'):
                        value = float(np.round(value, 4))
                    if key.endswith('num'):
                        value = int(value)

                    stats[profile_name]['dir'][key] = value

            elif key.endswith('size'):
                value = stats[profile_name][key]
                value = float(np.round(value, 4))
                stats[profile_name][key] = value
            elif key.endswith('num'):
                value = stats[profile_name][key]
                stats[profile_name][key] = value



    # save stats as a single yaml file
    with open('apero_disk_stats.yaml', 'w') as yfile:
        yaml.dump(stats, yfile, default_flow_style=False)



# =============================================================================
# End of code
# =============================================================================
