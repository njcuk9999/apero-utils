#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on ${CURRENT_YEAR}-${CURRENT_MONTH}-${CURRENT_DATE}

@author: cook
"""
from astropy.time import Time, TimeDelta
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
PATH = '/cosmos99/spirou/apero-data/spirou_008/log/tool/other/APEROL-PID-00017500796841586900-CRX7_apero_processing.log'

# =============================================================================
# Define functions
# =============================================================================



# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # load log file
    with open(PATH, 'r') as file:
        log_content = file.readlines()

    # loop around each line and get out the timestamp
    timestamps = []
    for line in tqdm(log_content):
        # split the line by the pipe character
        parts = line.split('|')
        # if the line has at least 2 parts, the first part is the timestamp
        if len(parts) >= 2:
            timestamp = parts[0]
            # we actually want an mjd date which just has the date and hour
            # so we split the timestamp by the space character
            date_part = parts[0].split(' ')[0]
            hour_part = parts[0].split(' ')[1].split(':')[0]
            # convert into a mjd using astropy Time
            iso_time = f"{date_part} {hour_part}:00:00"
            # append the iso time to the timestamps list
            timestamps.append(iso_time)

    # convert timestampes to MJD
    print('Converting {0} timestamps to MJD...'.format(len(timestamps)))
    times = Time(timestamps, format='iso', scale='utc')

    # now we have all the times lets plot them as a histogram with
    # bins of hours between minimum and maximum times

    # get the minimum and maximum times
    min_time = times.min()
    max_time = times.max()
    # Pad the end to include the final hour fully
    max_time += TimeDelta(3600, format='sec')
    # Create bins: list of Time objects spaced by 1 hour
    bin_edges = min_time + TimeDelta(np.arange(0, (max_time - min_time).sec + 3600, 3600), format='sec')

    # Convert to datetime for matplotlib
    t_datetime = times.datetime
    bin_edges_datetime = bin_edges.datetime

    # plot the histogram 
    plt.hist(t_datetime, bins=bin_edges_datetime, edgecolor='black', alpha=0.7)
    plt.xlabel('MJD Time')
    plt.ylabel('Count')
    plt.title('Histogram of times from log file (binned by hour)')
    plt.xticks(rotation=45)
    plt.grid(axis='y', alpha=0.75)
    # log the y-axis
    plt.yscale('log')
    plt.gcf().autofmt_xdate()

    plt.tight_layout()
    plt.show()  


# =============================================================================
# End of code
# =============================================================================
