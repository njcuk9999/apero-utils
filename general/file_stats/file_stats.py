#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-12-10 at 13:14

@author: cook
"""
import os
import sys
import glob

from tqdm import tqdm

from apero.core import constants

# =============================================================================
# Define variables
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def get_dir_stats(path: str):
    nfiles, nsizes = [], []

    # get a list of all files in path
    path_files = []
    for root, dirs, files in tqdm(os.walk(path)):
        for basename in files:
            filename = os.path.join(root, basename)
            path_files.append(filename)

    # now we do the counting
    for filename in path_files:
        # get real path (including symlinks)
        filename = os.path.realpath(filename)
        # append the
        nfiles.append(filename)
        nsizes.append(os.path.getsize(filename) / (1024 ** 2))

    return nfiles, nsizes


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # get parameters
    params = constants.load()


    # lets start with raw get all directories in raw
    directories = glob.glob(params['DRS_DATA_RAW'] + '/*')

    dir_files, dir_filesizes = dict(), dict()
    dir_nums, dir_sizes = dict(), dict()
    dir_filesize_min, dir_filesize_max = dict(), dict()
    # loop around raw direcr
    for directory in directories:
        print('*' * 50)
        print(directory)
        print('*' * 50)

        print('Counting files...')
        nfile, nsize = get_dir_stats(directory)

        if len(nfile) == 0:
            continue

        dir_files[directory] = nfile
        dir_filesizes[directory] = nsize
        dir_nums[directory] = len(nfile)
        dir_sizes[directory] = sum(nsize)
        dir_filesize_min[directory] = min(nsize)
        dir_filesize_max[directory] = max(nsize)

        print('Number of files: {0}'.format(dir_nums[directory]))
        print('Total size: {0} Gb'.format(dir_sizes[directory]/1024))
        print('Average size: {0} Mb'.format(dir_sizes[directory] / dir_nums[directory]))
        print('Min size: {0} Mb'.format(dir_filesize_min[directory]))
        print('Max size: {0} Mb'.format(dir_filesize_max[directory]))



    # print stats across all directories
    print('*' * 50)
    print('Number of directories: {0}'.format(len(directories)))
    print('Total number of files: {0}'.format(sum(dir_nums.values())))
    print('Total size: {0} Tb'.format(sum(dir_sizes.values()) / (1024 ** 2)))
    print('Average files per directory: {0}'.format(sum(dir_nums.values()) / len(directories)))
    print('Average size: {0} Mb'.format(sum(dir_sizes.values()) / sum(dir_nums.values())))
    print('Average directory size {0} Gb'.format(sum(dir_sizes.values()) / (1024 * len(directories))))
    print('Min size: {0} Mb'.format(min(dir_filesize_min.values())))
    print('Max size: {0} Mb'.format(max(dir_filesize_max.values())))



# =============================================================================
# End of code
# =============================================================================
