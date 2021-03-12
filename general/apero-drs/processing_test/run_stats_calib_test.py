#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-03-12

@author: cook
"""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Define variables
# =============================================================================
# get filename
FILENAME = ('/data/spirou/data/07000/runs/stats/'
            'APEROL-PID-00016155713574718550-A9SR_apero_processing_stats.txt')
# get night directory
NIGHTS = '/data/spirou/data/07000/tmp/'
# list of valid calibration recipes
CALIB_RECIPES = ['cal_badpix_spirou.py', 'cal_loc_spirou.py',
                 'cal_shape_spirou.py', 'cal_flat_spirou.py',
                 'cal_thermal_spirou.py', 'cal_wave_night_spirou.py',
                 'cal_shape_master_spirou.py', 'cal_wave_master_spirou.py',
                 'cal_dark_master.py']
# short names (for plot) for calib recipes
SHORTNAMES = ['BAD', 'LOC', 'SHAPEL', 'FLAT', 'THEMAL', 'WAVEL', 'SHAPEM',
              'WAVEM', 'DARKM']

# =============================================================================
# Define functions
# =============================================================================


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    # find nights
    allfiles = Path(NIGHTS).rglob('*')
    nights = []
    for filename in allfiles:
        # only consider directories
        if filename.is_dir():
            # only consider directories with files
            if len(list(filename.glob('*.fits'))) > 0:
                # just want the night
                night = str(filename).split(NIGHTS)[-1]
                nights.append(night)
    # ----------------------------------------------------------
    # now we have nights look for recipes for each night
    # ----------------------------------------------------------
    # open file
    with open(FILENAME, 'r') as wfile:
        recipe_runs = wfile.readlines()
    # storage
    count = dict()
    runs = dict()
    recipes = dict()
    for night in nights:
        # count for this night
        count[night] = dict()
        runs[night] = dict()
        # loop around recipe runs
        for recipe_run in recipe_runs:
            # skip runs not for this night
            if night not in recipe_run:
                continue
            # get recipe
            recipe = recipe_run.split(' ')[0].strip()
            # either add recipe to dict or update count
            if recipe not in count[night]:
                count[night][recipe] = 1
                runs[night][recipe] = [recipe_run.strip('\n')]
            else:
                count[night][recipe] += 1
                runs[night][recipe] += [recipe_run.strip('\n')]
            # deal with recipe
            if recipe not in recipes:
                recipes[recipe] = 1
            else:
                recipes[recipe] += 1

    # do a plot per night
    plt.close()
    fig, frames = plt.subplots(nrows=int(np.ceil(len(nights) / 2)), ncols=2)

    used = []
    # loop around nights (one plot per night)
    for it in range(len(nights)):
        # get the night name
        night = nights[it]
        # even and odd split
        if it % 2 == 0:
            frame = frames[it // 2][0]
            used += [[it // 2, 0]]
        else:
            frame = frames[it // 2][1]
            used += [[it // 2, 1]]

        # get the yvalues for this night

        xvalues, yvalues = [], []
        for cit, recipe in enumerate(CALIB_RECIPES):
            xvalues.append(SHORTNAMES[cit])
            if recipe in count[night]:
                yvalues.append(count[night][recipe])
            else:
                yvalues.append(0)

        frame.bar(xvalues, yvalues)
        frame.set(title=nights[it])

    for it in range(len(frames)):
        for jt in range(2):
            if [it, jt] not in used:
                frames[it][jt].axis('off')

    plt.show()
    plt.close()

# =============================================================================
# End of code
# =============================================================================
