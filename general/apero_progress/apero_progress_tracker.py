#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on ${CURRENT_YEAR}-${CURRENT_MONTH}-${CURRENT_DATE}

@author: cook
"""
import os
from astropy.time import Time, TimeDelta
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
CASE = 'spirou_offline'

if CASE == 'spirou_offline':
    # the PID of the APERO processing run
    # Get the PID from here: /cosmos99/spirou/apero-data/spirou_offline/msg/tool/other
    APERO_PID = 'PID-00017395509988958830-L4FG'
    # the working directory where the log and report files are located
    WORKING_DIR = '/cosmos99/spirou/apero-data/spirou_offline/'
    # the paths to the log and report files
    LOG_DIR = 'msg'  # or "log"
    # Has date
    has_date = False
elif CASE == 'spirou_008':
    # the PID of the APERO processing run
    # Get the PID from here: /cosmos99/spirou/apero-data/spirou_offline/msg/tool/other
    APERO_PID = 'PID-00017500796841586900-CRX7'
    # the working directory where the log and report files are located
    # WORKING_DIR = '/cosmos99/spirou/apero-data/spirou_offline/'
    WORKING_DIR = '/cosmos99/spirou/apero-data/spirou_008/'
    # the paths to the log and report files
    LOG_DIR = 'log'  # or "msg"
    # Has date
    has_date = True
else:
    raise ValueError("Invalid CASE. Choose 'spirou_offline' or 'spirou_008'.")
# the paths to the log and report files
PATH_TO_LOG = os.path.join(WORKING_DIR, LOG_DIR, 'tool', 'other')
LOG_FILE = 'APEROL-{APERO_PID}_apero_processing.log'
# the path to the report file
PATH_TO_REPORT = os.path.join(WORKING_DIR, LOG_DIR, 'report', 'processing')
REPORT_FILE = '{APERO_PID}_apero_processing_ids.txt'

# =============================================================================
# Define functions
# =============================================================================



# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # get log file
    log_filename = os.path.join(PATH_TO_LOG, LOG_FILE.format(APERO_PID=APERO_PID))
    # load log file
    with open(log_filename, 'r') as file:
        log_content = file.readlines()

    counter = 0
    # get creation date of log file
    creation_time = os.path.getctime(log_filename)

    ctime = Time(creation_time, format='unix', scale='utc').fits
    date_part = ctime.split('T')[0]
    running_hour = int(ctime.split('T')[1].split(':')[0])

    # loop around each line and get out the timestamp
    timestamps = []
    for line in tqdm(log_content):
        # split the line by the pipe character
        parts = line.split('|')
        # if the line has at least 2 parts, the first part is the timestamp
        if len(parts) >= 2:
            timestamp = parts[0]
            # if has_date is False we need to get the date
            if has_date:
                # we actually want an mjd date which just has the date and hour
                # so we split the timestamp by the space character
                date_part = parts[0].split(' ')[0]
                hour_part = parts[0].split(' ')[1].split(':')[0]
            else:
                # need to figure out if the hour part has jumped to the next day
                hour_part = int(timestamp.split(':')[0])

                if hour_part < running_hour:
                    # increment the date part by one day
                    date_part = Time(date_part, format='iso').mjd + 1
                    date_part = Time(date_part, format='mjd').fits.split('T')[0]
                # update the running hour
                running_hour = int(hour_part)

            # convert into a mjd using astropy Time
            iso_time = f"{date_part} {str(hour_part).zfill(2)}:00:00"
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
    plt.xlabel('Time')
    plt.ylabel('Number of messages printed per hour')
    plt.title('Histogram of messages from log file (binned by hour)')
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
