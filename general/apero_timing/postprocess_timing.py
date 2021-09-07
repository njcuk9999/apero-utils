#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-08-31

@author: cook
"""
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import textwrap
from tqdm import tqdm


# =============================================================================
# Define variables
# =============================================================================
case = 3


if case == 1:
    # minimum time to filter by
    min_time = 2
    # log dir
    workspace = '/data/spirou/data/07000/msg/processing/'
    logdir = 'APEROG-PID-00016300669142438500-K5ZK_apero_processing_group'
    workspace = os.path.join(workspace, logdir)
    # whether plots are log y
    log_y = False
elif case == 2:
    # minimum time to filter by
    min_time = 120
    # log dir
    workspace = '/spirou/drs-data/full_210802/msg/processing/'
    logdir = 'APEROG-PID-00016280226609722440-67AO_apero_processing_group'
    workspace = os.path.join(workspace, logdir)
    # whether plots are log y
    log_y = True
elif case == 3:
    # minimum time to filter by
    min_time = 60
    # log dir
    workspace = '/spirou/drs-data/full_210802/msg/processing/'
    logdir = 'APEROG-PID-00016304381051960400-5LOP_apero_processing_group'
    workspace = os.path.join(workspace, logdir)
    # whether plots are log y
    log_y = False


# =============================================================================
# Define functions
# =============================================================================
class LogProfile:
    def __init__(self, filename):
        self.filename = filename
        self.times = []
        self.lines = []
        self.deltas = [0]
        self.start = 0

    def __str__(self):
        return 'LogProfile[{0}]'.format(os.path.basename(self.filename))

    def __repr__(self):
        return self.__str__()

    def get_start(self):
        basename = os.path.basename(self.filename)
        unix = float(basename.split('PID-')[-1].split('-')[0])/1e7
        self.start = unix

    def get(self):

        self.get_start()

        with open(self.filename, 'r') as lfile:
            self.lines = lfile.readlines()
        # get timings
        for line in self.lines:
            raw_time = line.split('-')[0]
            self.times.append(pd.Timestamp(raw_time))
        # find deltas
        for it in range(1, len(self.lines)):

            time0 = self.times[it - 1]
            time1 = self.times[it]

            dt = (time1-time0).total_seconds()
            if dt < 0:
                dt += 24*3600
            self.deltas.append(dt)

    def longest_steps(self, num: int, return_pos: bool = False):

        # sort largest to smallest
        argsort = np.argsort(self.deltas)[::-1]
        # for return
        positions = []
        # get "num" rows that are the longest
        for it in range(num):
            pos = argsort[it]


            if not return_pos:
                print('='*50)
                print('\tdt = {0}'.format(self.deltas[pos]))
                print('='*50)
                print(self.lines[pos-1])
                print(self.lines[pos])
            else:
                positions.append(pos)

        return positions


def clean_log(text):

    text = text[31:]
    # try to remove odocode if present
    if 'o_pp' in text:
        before, after = text.split('_pp')
        # remove odocode
        before = before[:-8]

        text = '{0}{{ODOCODE}}_pp{1}'.format(before, after)
    elif '_pp' in text:
        before, after = text.split('_pp')
        # remove odocode
        before = before[:-10]

        text = '{0}{{HASHCODE}}_pp{1}'.format(before, after)
    elif '_e2ds' in text:
        before, after = text.split('_e2ds')
        # remove odocode
        before = before[:-10]
        text = '{0}{{HASHCODE}}_e2ds{1}'.format(before, after)

    text = text.replace('\t', '')
    text = text.replace('\n', '')

    if 'EXT=' in text:
        before, after = text.split('EXT=')

        while after[0].isdigit():
            after = after[1:]

        text = '{0}EXT=N{1}'.format(before, after)


    return '\n'.join(textwrap.wrap(text, 40))


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # get the log files
    logfiles = glob.glob(os.path.join(workspace, '*apero_postprocess*'))
    # -------------------------------------------------------------------------
    # make log classes
    # -------------------------------------------------------------------------
    logclasses = []
    # loop around all files
    for it in tqdm(range(len(logfiles))):
        # get log file
        logfile = logfiles[it]
        # make class
        logclass = LogProfile(logfile)
        logclass.get()
        # save
        logclasses.append(logclass)

    for logclass in logclasses:
        normx = np.linspace(0, 1, len(logclass.deltas))
        plt.plot(normx, logclass.deltas)
    # -------------------------------------------------------------------------
    # get all events over 120s in length
    # -------------------------------------------------------------------------
    befores, afters = [], []
    combineds, deltas = [], []
    for logclass in logclasses:
        for jt in range(len(logclass.lines)):
            if logclass.deltas[jt] > min_time:
                deltas.append(logclass.deltas[jt])
                before = clean_log(logclass.lines[jt - 1])
                after = clean_log(logclass.lines[jt])
                combine = '>>{0}\n>>{1}'.format(before, after)
                befores.append(before)
                afters.append(after)
                combineds.append(combine)
    # numpy array lists
    combineds = np.array(combineds)
    deltas = np.array(deltas)
    # -------------------------------------------------------------------------
    # get unique counts and mean times
    # -------------------------------------------------------------------------
    utimes_med, ucounts = dict(), dict()
    utimes_min, utimes_max = dict(), dict()
    utimes_mean = dict()
    for utag in set(combineds):
        mask = utag == combineds
        utime_med = np.median(deltas[mask])
        utime_min = np.min(deltas[mask])
        utime_max = np.max(deltas[mask])
        utime_mean = np.mean(deltas[mask])
        if utime_med > min_time:
            ucounts[utag] = np.sum(mask)
            utimes_med[utag] = utime_med
            utimes_min[utag] = utime_min
            utimes_max[utag] = utime_max
            utimes_mean[utag] = utime_mean
    # -------------------------------------------------------------------------
    # plot 1: the utimes vs ustats + utimes vs ucounts
    # -------------------------------------------------------------------------
    plt.close()
    fig, frames = plt.subplots(ncols=2, nrows=1)

    frames[0].bar(list(utimes_max.keys()), list(utimes_max.values()),
                  label='max', color='orange')
    frames[0].bar(list(utimes_max.keys()), list(utimes_min.values()),
                  label='min', color='g')
    frames[0].bar(list(utimes_max.keys()), list(utimes_med.values()),
                  label='median', edgecolor='r', color='None')
    frames[0].bar(list(utimes_max.keys()), list(utimes_mean.values()),
                  label='mean', edgecolor='b', color='None')
    frames[0].legend(loc=0)
    frames[0].set(ylabel='Time / s')
    if log_y:
        frames[0].set_yscale('log')
    plt.setp(frames[0].xaxis.get_majorticklabels(), rotation=90, fontsize=10)

    frames[1].bar(list(utimes_max.keys()), list(ucounts.values()))
    frames[1].set(ylabel='Number of counts')
    if log_y:
        frames[1].set_yscale('log')
    plt.setp(frames[1].xaxis.get_majorticklabels(), rotation=90, fontsize=10)

    title = 'Log lines with dt>{0} s \n Between line A and A+1'
    plt.suptitle(title.format(min_time))
    plt.subplots_adjust(bottom=0.5, left=0.05, right=0.975, top=0.9, wspace=0.15)
    plt.show()
    plt.close()
    # -------------------------------------------------------------------------
    # print longest duration lines
    # -------------------------------------------------------------------------
    for logclass in logclasses:
        for jt in range(len(logclass.times)):
            if logclass.deltas[jt] > 20000:
                print("-" * 50)

                if jt > 4:
                    print(logclass.lines[jt - 4])
                if jt > 3:
                    print(logclass.lines[jt - 3])
                if jt > 2:
                    print(logclass.lines[jt - 2])
                if jt > 1:
                    print(logclass.lines[jt - 1])
                print(':::GAP:::\n')
                print(logclass.lines[jt])
                if jt < len(logclass.times) - 1:
                    print(logclass.lines[jt + 1])
                if jt < len(logclass.times) - 2:
                    print(logclass.lines[jt + 2])
                if jt < len(logclass.times) - 3:
                    print(logclass.lines[jt + 3])
                if jt < len(logclass.times) - 4:
                    print(logclass.lines[jt + 4])
    # -------------------------------------------------------------------------
    # get total times and start times as lists
    # -------------------------------------------------------------------------
    total_times = []
    start_times = []
    for logclass in logclasses:
        total_times.append(np.sum(logclass.deltas))
        start_times.append(logclass.start)

    med_time = np.median(total_times)
    std_time = np.std(total_times)
    print('Median time: {0} +/- {1}'.format(med_time, std_time))

    start_times = np.array(start_times) - np.min(start_times)
    # -------------------------------------------------------------------------
    # plot 1: the utimes vs ustats + utimes vs ucounts
    # -------------------------------------------------------------------------
    fig, frame = plt.subplots()
    frame.scatter(start_times, total_times)
    frame.set(xlabel='Time since first [s]', ylabel='Total time taken [s]')
    plt.show()
    plt.close()

# =============================================================================
# End of code
# =============================================================================
