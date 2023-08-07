#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-14 at 10:32

@author: cook
"""

import os

# =============================================================================
# Define variables
# =============================================================================
name1 = 'spirou@rali'
name2 = 'cook@jupiter'
name3 = 'cook@nb19'
name4 = 'spirou@maestria'
name5 = 'newworlds'
name6 = 'lam'
name7 = 'halau'

ref_name = 'spirou@rali'

rawpaths = dict()
rawpaths[name1] = '/spirou2/apero-data/drs-data/spirou_minidata2_07286/red/'
rawpaths[name2] = '/scratch2/spirou/drs-data/spirou_minidata2_07286_jupiter/red/'
rawpaths[name3] = '/scratch/apero-nb19/drs-data/spirou_minidata2_07286_nb19/red/'
rawpaths[name4] = '/spirou2/apero-data/drs-data/spirou_minidata2_07286_maestria/red/'

rsync_info = dict()
rsync_info[name1] = 'spirou@rali'
rsync_info[name2] = 'cook@jupiter'
rsync_info[name3] = 'cook@nb19'
rsync_info[name4] = 'spirou@maestria'

inpaths = dict()
inpaths[name1] = '/scratch2/spirou/drs-data/spirou_minidata2_07286_rali/red/'
inpaths[name2] = '/scratch2/spirou/drs-data/spirou_minidata2_07286_jupiter/reduced/'
inpaths[name3] = '/scratch2/spirou/drs-data/spirou_minidata2_07286_nb19/red/'
inpaths[name4] = '/scratch2/spirou/drs-data/spirou_minidata2_07286_maestria/red/'
# inpaths[name5] = '/scratch2/spirou/drs-data/minidata2_07275_newworld/red/'
# inpaths[name6] = '/scratch2/spirou/drs-data/minidata2_07275_lam/red/'
# inpaths[name7] = '/scratch2/spirou/drs-data/minidata2_07275_halau/red/'

# -----------------------------------------------------------------------------
rsync_cmd = 'rsync -avu {0}:{1} {2}'
TEST = False

# =============================================================================
# Define functions
# =============================================================================


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # loop around names
    for name in rawpaths:
        # get rsync info for this name
        rsyncinfo = rsync_info[name]
        # get rawpath for this name
        rawpath = rawpaths[name]
        # get inpath for this name
        inpath = inpaths[name]
        # make inpath if it doesn't exist
        if not os.path.exists(inpath):
            os.makedirs(inpath)
        # create rsync command
        command = rsync_cmd.format(rsyncinfo, rawpath, inpath)
        # run comamnd
        print('\n\nRunning command: {0}'.format(command))
        if not TEST:
            os.system(command)


# =============================================================================
# End of code
# =============================================================================
