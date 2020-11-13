#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CODE DESCRIPTION HERE

Created on 2020-11-2020-11-04 12:07

@author: cook
"""
import numpy as np
import sys
from typing import Dict, List, Tuple


# =============================================================================
# Define variables
# =============================================================================
# default log file
LOGFILE = 'APEROL-PID-00016044429913132070-Z7JV_apero_processing.log'


# =============================================================================
# Define functions
# =============================================================================
def sort_log_rows(lines: List[str]) -> Tuple[List[str], List[float]]:
    """
    Sort a log file into valid recipes + their timings

    :param lines: list of strings from open(logfile)

    :return: list of valid recipes and list of timings for each recipe
    """
    # we want to get the lines after the line with 'ID = {0} TIME = {1}'
    valid_recipes = []
    valid_times = []
    # loop through lines
    for row, line in enumerate(lines):
        if 'ID = ' not in line:
            continue
        if 'Time = ' not in line:
            continue
        # split at \t\t
        rline = lines[row + 1]
        recipe = rline.split('\t\t')[-1].split('.py')[0]
        # add next line
        valid_recipes.append(recipe)
        # add time
        tline = line.split('Time = ')[-1].split('\n')[0]
        valid_times.append(float(tline))
    return valid_recipes, valid_times


def get_timings(valid_recipes: List[str], valid_times: List[float]
                ) -> Tuple[Dict[str, List[float]], Dict[str, float],
                           Dict[str, float]]:
    """
    Sort the valid recipes into groups (grouped by recipe name) and
    provide the average time for each recipe

    :param valid_recipes: list of recipes
    :param valid_times: list of timings for each recipe run

    :return: tuple, 1. dict - the groups recipes + times, 2. dict the mean
             timings for each recipe
    """
    # storage
    timings = dict()
    # push these into a dictionary of recipes
    for row in range(len(valid_recipes)):
        # get this iterations values
        recipe = valid_recipes[row]
        time = valid_times[row]
        # don't count ones where time is less than 5 seconds
        if time < 5:
            continue
        # add to timings
        if recipe in timings:
            timings[recipe] += [time]
        else:
            timings[recipe] = [time]
    # get average times
    mean_timings = dict()
    med_timings = dict()
    for recipe in timings:
        mean_timings[recipe] = float(np.mean(timings[recipe]))
        med_timings[recipe] = float(np.median(timings[recipe]))
    # return timings and mean timings
    return timings, mean_timings, med_timings


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    # deal with args
    if len(sys.argv) == 2:
        logfile = sys.argv[1]
    else:
        logfile = LOGFILE
    # open logf ile
    with open(logfile, 'r') as lfile:
        _lines = lfile.readlines()
    # get valid recipes and times from log file
    vrecipes, vtimes = sort_log_rows(_lines)
    # get timing stats
    times, meantimes, medtimes = get_timings(vrecipes, vtimes)
    # print out average timings (in minutes)
    print('=' * 70)
    print('Average timings')
    print('=' * 70)
    for _recipe in meantimes:
        pargs = [_recipe, meantimes[_recipe], meantimes[_recipe] / 60.0,
                 len(times[_recipe])]
        msg = '{0:30s}: {1:7.2f} seconds   ({2:5.2f} minutes)       N = {3}'

        print(msg.format(*pargs))
    # print out median timings (in minutes)
    print('=' * 70)
    print('Median timings')
    print('=' * 70)
    for _recipe in medtimes:
        pargs = [_recipe, medtimes[_recipe], medtimes[_recipe] / 60.0,
                 len(times[_recipe])]
        msg = '{0:30s}: {1:7.2f} seconds   ({2:5.2f} minutes)       N = {3}'
        print(msg.format(*pargs))


# =============================================================================
# End of code
# =============================================================================
