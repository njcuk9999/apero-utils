#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-03-30

@author: cook
"""
from astropy.table import Table
from astropy.time import Time, TimezoneInfo
from astropy import units as uu
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from tqdm import tqdm


# =============================================================================
# Define variables
# =============================================================================
EST = TimezoneInfo(-4 * uu.hour)

TODAY = '2021-03-31'

WORKSPACE = '/spirou/cook/db_test'

LOG_DIR = '/spirou/drs-data/mini-data-07000/msg/processing/APEROG*'

LOG_STR = '-!|'

# =============================================================================
# Define functions
# =============================================================================


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    # file for processes
    PROCESS_FILE = os.path.join(WORKSPACE, 'process.txt')
    # file for connections
    CONNECT_FILE = os.path.join(WORKSPACE, 'connect.txt')

    # load process_file
    process_table = Table.read(PROCESS_FILE, format='ascii.csv')
    connect_table = Table.read(CONNECT_FILE, format='ascii.csv')

    ptime = process_table['unix_time']
    pcount = process_table['count']
    ccount = connect_table['count']

    # get log files
    log_files = glob.glob(os.path.join(LOG_DIR, '*.log'))

    error_times = []
    for log_file in tqdm(log_files):
        with open(log_file, 'r') as lfile:
            lines = lfile.readlines()
        for line in lines:
            if LOG_STR in line:
                raw_time = line.split(LOG_STR)[0]
                human_time = '{0} {1}'.format(TODAY, raw_time)
                error_times.append(Time(human_time, format='iso').unix)
                break

    # get the shortest array - for when we stop between the two counters
    minlen = np.min([len(pcount), len(ccount)])
    ptime = ptime[:minlen]
    pcount = pcount[:minlen]
    ccount = ccount[:minlen]

    fig, frames = plt.subplots(ncols=1, nrows=2, sharex='all')

    for tx in range(len(error_times)):
        frames[0].axvline(error_times[tx], alpha=0.1, color='red')
        frames[1].axvline(error_times[tx], alpha=0.1, color='red')
    frames[0].plot(ptime, pcount)
    frames[1].plot(ptime, ccount)

    frames[0].set(xlabel='unix time', ylabel='Number of connections to "spirou" database')
    frames[1].set(xlabel='unix_time', ylabel='Number of "netstat -nt" entries')

    title = 'Max processes: {0}   Max netstat: {1}'
    targs = [np.max(pcount), np.max(ccount)]
    plt.suptitle(title.format(*targs))

    plt.show()
    plt.close()

# =============================================================================
# End of code
# =============================================================================
