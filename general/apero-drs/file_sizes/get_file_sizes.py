#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-04-30

@author: cook
"""
import numpy as np
import glob
import os

# =============================================================================
# Define variables
# =============================================================================
# search key
SEARCH_KEY = '*/DEBUG*'
# total size of groups of files to report in bytes
SIZE_TO_REPORT = 50 * (1024 ** 3)  # 50 GB


# =============================================================================
# Define functions
# =============================================================================
def human_readable_size(size, decimal_places=3):
    """
    Get human readible size

    Taken from:
        https://www.codegrepper.com/code-examples/python/
        python+get+human+readable+file+size

    :param size:
    :param decimal_places:
    :return:
    """
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            break
        size /= 1024.0
    msg = "{size:.{decimal_places}f}{unit}"
    mkwargs = dict(size=size, decimal_places=decimal_places,
                   unit=unit)
    return msg.format(**mkwargs)


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # ------------------------------------------------------------------
    # get all files
    print('\nFinding all "{0}" files'.format(SEARCH_KEY))
    files = glob.glob(SEARCH_KEY)
    print('\t{0} "{1}"" files found'.format(len(files), SEARCH_KEY))

    # get list of all prefix and suffix
    print('\nGetting all file combinations (prefix*suffix)')
    combinations = []
    for filename in files:
        prefix = os.path.basename(filename).split('2')[0]
        suffix = os.path.basename(filename).split('pp')[-1]
        # append to list
        combinations.append(prefix + '*' + suffix)
    # ------------------------------------------------------------------
    # find unique combinations of prefix + suffix
    ucombinations = np.unique(combinations)
    # ------------------------------------------------------------------
    # get the files sorted by the unique combinations
    sorted_count = dict()
    sorted_filenames = dict()
    sorted_sizes = dict()
    sorted_units = dict()
    for ucombination in ucombinations:
        print('\nDealing with combination: "{0}"'.format(ucombination))
        sorted_count[ucombination] = 0
        sorted_filenames[ucombination] = []
        sorted_sizes[ucombination] = []
        sorted_units[ucombination] = []
        # loop around files and put in groups
        for filename in files:
            prefix = ucombination.split('*')[0]
            suffix = ucombination.split('*')[-1]
            # only put in group with correct prefix / suffix
            if prefix in filename and suffix in filename:
                sorted_count[ucombination] += 1
                sorted_filenames[ucombination].append(filename)
                # get file size
                rawsize = os.path.getsize(filename)
                sorted_sizes[ucombination].append(rawsize)
    # ------------------------------------------------------------------
    # display the average size of each (and how many we have of each)
    for ucombination in ucombinations:
        # get count
        count = sorted_count[ucombination]
        # get raw sizes average / total
        average_raw = np.mean(sorted_sizes[ucombination])
        total_raw = np.sum(sorted_sizes[ucombination])
        # display these in human formats
        average_human = human_readable_size(average_raw)
        total_human = human_readable_size(total_raw)
        # only display for total size larger than the size to report
        if total_raw > SIZE_TO_REPORT:
            print('\n\n' + ucombination)
            print('\t Number of files = {0}'.format(count))
            print('\t Average size of files = {0}'.format(average_human))
            print('\t Total size of files = {0}'.format(total_human))

# =============================================================================
# End of code
# =============================================================================
