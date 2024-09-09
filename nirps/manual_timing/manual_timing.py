#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-02-12 at 18:36

@author: cook
"""
import os
from astropy.table import Table
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Define variables
# =============================================================================
LOG_OVERRIDE = r'H:\nirps_he\misc\log'
# -----------------------------------------------------------------------------
LOG_NAME = 'nirps_he_profile_online.log'
PROFILE_NAME = 'nirps he online'


# =============================================================================
# Define functions
# =============================================================================
def read_log(log_name: str) -> Table:
    # construct log path
    homedir = os.path.expanduser('~')
    # construct log path
    logpath = os.path.join(homedir, '.apero', 'manual_trigger')
    logpath = LOG_OVERRIDE
    # read log file
    logname = os.path.join(logpath, LOG_NAME)
    # read csv file
    if not os.path.exists(logname):
        raise FileNotFoundError('Log name {0} does not exist'.format(logname))
    # otherwise we read the log file
    logtable = Table.read(logname, format='ascii.csv', data_start=0,
                          names=['TIMESTAMP', 'PROFILE', 'PROCESS', 'OBSDIRS',
                                 'COMMENT'])
    # convert timestamp to astropy time
    logtable['TIMESTAMP'] = Time(logtable['TIMESTAMP'])
    # sort by time
    logtable.sort('TIMESTAMP')
    # return log table
    return logtable


def get_stats(logtable: Table, profile_name: str) -> Table:
    # get unique profiles
    profiles = list(set(logtable['PROFILE']))
    # get unique obsdirs
    all_obsdirs = list(set(logtable['OBSDIRS']))
    # remove obsidrs with "|" in and '*'
    obsdirs = []
    for obsdir in all_obsdirs:
        if '|' in obsdir:
            continue
        if '*' in obsdir:
            continue
        obsdirs.append(obsdir)
    # get mask for this profile
    profile_mask = logtable['PROFILE'] == profile_name
    # check for no entries
    if np.sum(profile_mask) == 0:
        return Table()
    # get profile table
    profile_table = logtable[profile_mask]
    # create dict for profile stats
    profile_stats = dict()
    profile_stats['OBSDIR'] = []
    profile_stats['MANUAL_START'] = []
    profile_stats['MANUAL_END'] = []
    profile_stats['MANUAL_DURATION'] = []
    profile_stats['APERO_START'] = []
    profile_stats['APERO_END'] = []
    profile_stats['APERO_DURATION'] = []
    profile_stats['ARI_START'] = []
    profile_stats['ARI_END'] = []
    profile_stats['ARI_DURATION'] = []
    # loop around each obsdir
    for obsdir in obsdirs:
        # get mask for this obsdir
        obsdir_mask = profile_table['OBSDIRS'] == obsdir
        # check for no entries
        if np.sum(obsdir_mask) == 0:
            continue
        # get the table
        obsdir_table = profile_table[obsdir_mask]
        # only continue if we have all stats
        has_all_stats = True
        # storage for the start, end and duration
        timing_dict = dict()
        # loop around and check pairs of keys
        for key in ['MANUAL', 'APERO', 'ARI']:
            # find the first entry for 'MANUAL_START'
            start_mask = obsdir_table['PROCESS'] == f'{key}_START'
            end_mask = obsdir_table['PROCESS'] == f'{key}_END'
            # check for no entries
            if np.sum(start_mask) == 0 or np.sum(end_mask) == 0:
                has_all_stats = False
                break
            # get the first start time
            start_time = obsdir_table['TIMESTAMP'][start_mask][0]
            end_time = obsdir_table['TIMESTAMP'][end_mask][-1]
            # get the duration
            duration = end_time - start_time
            # convert to hours
            duration = duration.to('hour').value
            # push into timing dict
            timing_dict[f'{key}_START'] = start_time
            timing_dict[f'{key}_END'] = end_time
            timing_dict[f'{key}_DURATION'] = duration
        # if we have all stats add to profile stats
        if has_all_stats:
            profile_stats['OBSDIR'].append(obsdir)
            for key in ['MANUAL', 'APERO', 'ARI']:
                for subkey in ['START', 'END', 'DURATION']:
                    pkey = f'{key}_{subkey}'
                    profile_stats[pkey].append(timing_dict[pkey])
    # convert obsdir key to a astropy time array
    profile_stats['OBSDIR'] = Time(profile_stats['OBSDIR'])
    # finally return stats as a table per profile
    return Table(profile_stats)


def plot_stats(stats_dict: Table, profile_name: str):

    fig, frame = plt.subplots(ncols=1, nrows=1)

    frame.set(xlabel='Observation directory', ylabel='Hours',
              title=profile_name)

    for key in ['MANUAL', 'APERO', 'ARI']:
        # get the start, end and duration
        obsdir = stats_dict['OBSDIR']
        duration = stats_dict[f'{key}_DURATION']
        # plot the start and end
        frame.plot_date(obsdir.plot_date, duration, label=f'{key} duration')
    # add legend
    frame.legend()
    # add a grid
    frame.grid()

    plt.show()


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # get log table
    _logtable = read_log(LOG_NAME)
    # get the stats
    _stats = get_stats(_logtable, PROFILE_NAME)
    # plot the stats
    plot_stats(_stats, PROFILE_NAME)


    print(_stats)

# =============================================================================
# End of code
# =============================================================================
