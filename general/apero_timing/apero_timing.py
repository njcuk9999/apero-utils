#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-07-21

@author: cook
"""
from astropy.time import Time
from datetime import datetime, timedelta
from astropy import units as uu
import glob
import numpy as np
from tqdm import tqdm
from typing import Any, Dict, List, Tuple, Union


# =============================================================================
# Define variables
# =============================================================================
loggroup = 'APEROG-PID-00016267204261398790-URSX_apero_processing_group'


# =============================================================================
# Define functions
# =============================================================================
class LogFile:
    pid: Union[str, None]         # the log pid for this log file
    start: Union[Time, None]      # start time
    end: Union[Time, None]        # end time
    duration: Union[float, None]  # duration in seconds
    recipe: str                   # recipe name
    run_id: Union[str, None]      # run id if present
    success: bool = False

    def __init__(self, logfile: str):
        self.logfile = logfile
        self.pid = None
        self.unixtime = None
        # astrop
        self.start = None
        # astropy time (end time)
        self.end = None
        # duration in seconds
        self.duration = None
        # recipe and run id
        self.recipe = 'Unknown'
        self.run_id = None
        # whether recipe was successful
        self.success = False

    def __str__(self) -> str:
        if self.run_id is not None:
            return '{0}[{1}]:{2}'.format(self.recipe, self.run_id, self.pid)
        else:
            return '{0}:{1}'.format(self.recipe, self.pid)

    def __repr__(self) -> str:
        return self.__str__()

    def process(self):
        # open file
        with open(self.logfile, 'r') as logfile:
            lines = logfile.readlines()
        # ---------------------------------------------------------------------
        # get pid and unix time
        self.pid, self.unixtime = _get_pid(lines)
        # ---------------------------------------------------------------------
        # get times for all lines (date time as its way quicker)
        linetimes = _get_times(self.unixtime, lines)
        # get the start and end time
        self.start = Time(min(linetimes))
        self.end = Time(max(linetimes))
        self.duration = (self.end - self.start).to(uu.s).value
        # ---------------------------------------------------------------------
        # get recipe name and run id (if set)
        self.recipe = _get_recipe(lines)
        self.run_id = _get_run_id(lines)
        # get whether run successfully completed
        self.success = _get_success(lines)


# =============================================================================
# Functions to get stuff from log file
# =============================================================================
def _get_pid(lines: List[str]) -> Tuple[Union[str, None], Union[Time, None]]:
    # set up default values
    pid, unixtime = None, None
    # process pid and unix time
    for line in lines:
        if '@PID' in line:
            rawpid = line.split('@PID-')[-1].split()[0]
            pid = 'PID-{0}'.format(rawpid)
            strtime = rawpid.split('-')[0]
            strtime = strtime.lstrip('0')
            # raw float unix time
            rawunix = int(strtime) / 1e7
            # make into a unix time
            unixtime = Time(rawunix, format='unix')
            break
    return pid, unixtime


def _get_times(unixtime: Time, lines: List[str]) -> List[datetime]:
    """
    This is a hack we don't have the full time stamp for each line just
    the HH:MM:SS - have to guess the day from the PID - assumes PID
    is always generated before any time stamp in the log

    :param lines: list of strings, the lines from log file

    :return: list of times, the full time stamp for each line
    """
    # set current time to the unix time found before
    current_time = Time(unixtime)
    current_hour = current_time.datetime.hour
    current_day = current_time.iso.split()[0]
    # store date times
    datetimes = []
    # loop around each line
    for row, line in enumerate(lines):
        # get raw time string (no date)
        rawtimestr = line[:12]
        rawtimehour = line[:2].lstrip('0')
        if len(rawtimehour) == 0:
            rawtimehour = 0
        else:
            rawtimehour = int(rawtimehour)
        # work out guess at time stamp from raw time str and current day
        targs = [current_day, rawtimestr]
        tguess = datetime.fromisoformat('{0} {1}'.format(*targs))
        if rawtimehour < current_hour:
            dtime = tguess + timedelta(days=1)
        else:
            dtime = tguess

        current_hour = dtime.hour
        # append date times
        datetimes.append(dtime)

    # return list of times
    return datetimes


def _get_recipe(lines) -> str:
    recipe = lines[0].split('|')[1].split('[')[0]
    return recipe


def _get_run_id(lines: List[str]) -> Union[int, None]:
    # get program
    program = lines[0].split('|')[1]
    # look for id and case into int if possible
    if '[' in program and ']' in program:
        # noinspection PyBroadException
        try:
            return int(program.split('[')[1].split(']')[0].lstrip('0'))
        except Exception as _:
            return None
    else:
        return None


def _get_success(lines: List[str]) -> bool:
    if 'successfully completed' in lines[-2]:
        return True
    else:
        return False


# =============================================================================
# Stats functions
# =============================================================================
def get_stats(logs: List[LogFile]) -> Dict[str, Any]:
    durations = list(map(lambda x: x.duration, logs))
    start_times = list(map(lambda x: x.start, logs))
    end_times = list(map(lambda x: x.end, logs))
    stats = dict()
    stats['MEAN'] = np.nanmean(durations)
    stats['MEDIAN'] = np.nanmedian(durations)
    stats['STD'] = np.nanstd(durations)
    stats['NUMBER'] = len(durations)
    stats['CPU_TIME_SS'] = np.sum(durations)
    stats['TOTALTIME_SS'] = (max(end_times) - min(start_times)).to(uu.s).value
    stats['SPEED_UP'] = stats['CPU_TIME_SS'] / stats['TOTALTIME_SS']
    stats['MEANTIME_SS'] = stats['TOTALTIME_SS'] / stats['NUMBER']
    stats['CPU_TIME_HR'] = stats['CPU_TIME_SS'] / 3600
    stats['TOTALTIME_HR'] = stats['TOTALTIME_SS'] / 3600
    stats['MEANTIME_HR'] = stats['MEANTIME_SS'] / 3600
    return stats


def print_stats(recipe: str, stats: Dict[str, Any]):
    print()
    print('='*50)
    print('\t{0}'.format(recipe))
    print('='*50)
    # print stats
    statstr = ('\n\tMed Time: {MEDIAN:.3f} s +/- {STD:.3f}'
               '\n\tNruns: {NUMBER}'
               '\n\tTotal Time: {TOTALTIME_HR:.3f} hr'
               '\n\tTotal CPU Time: {CPU_TIME_HR:.3f} hr (x{SPEED_UP:.2f})'
               '\n')
    print(statstr.format(**stats))


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    # get log groups
    log_groups = glob.glob('APEROG*')
    # get log files
    files = glob.glob(log_groups[-1] + '/*.log')
    # store astropy times
    logs = dict()
    # print progress
    print('Analysing log files')
    # loop around log files
    for filename in tqdm(files):
        # construct log
        log = LogFile(filename)
        # process log
        log.process()
        # add to correct recipe
        if log.recipe in logs:
            logs[log.recipe].append(log)
        else:
            logs[log.recipe] = [log]
    # -------------------------------------------------------------------------
    # get stats and print them
    print()
    rstats = dict()
    for recipe in logs:
        # get stats
        rstats[recipe] = get_stats(logs[recipe])
        # print stats
        print_stats(recipe, rstats[recipe])

# =============================================================================
# End of code
# =============================================================================
