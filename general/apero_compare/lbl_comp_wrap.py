#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run all mini data sets through lbl one after the other

Created on 2021-10-18
Last updated 2022-09-26

@author: artigau, cook
"""
from lbl import lbl_wrap
from multiprocessing import Process

# =============================================================================
# Define functions
# =============================================================================
def main(name, data_dir):
    print('=' * 50)
    print('Running LBL for {0}'.format(name))
    print('=' * 50)
    # set up parameters
    rparams = dict()
    rparams['INSTRUMENT'] = 'SPIROU'
    rparams['DATA_SOURCE'] = 'APERO'
    rparams['DATA_DIR'] = data_dir
    rparams['DATA_TYPES'] = ['SCIENCE']
    rparams['OBJECT_SCIENCE'] = ['GL699']
    rparams['OBJECT_TEMPLATE'] = ['GL699']
    rparams['OBJECT_TEFF'] = [3224]
    rparams['RUN_LBL_TELLUCLEAN'] = False
    rparams['RUN_LBL_TEMPLATE'] = True
    rparams['RUN_LBL_MASK'] = True
    rparams['RUN_LBL_COMPUTE'] = True
    rparams['RUN_LBL_COMPILE'] = True
    rparams['SKIP_LBL_TEMPLATE'] = True
    rparams['SKIP_LBL_MASK'] = True
    rparams['SKIP_LBL_COMPUTE'] = True
    rparams['SKIP_LBL_COMPILE'] = True
    # run the wrap
    lbl_wrap(rparams)

# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # define names
    name1 = 'spirou@rali'
    name2 = 'cook@jupiter'
    name3 = 'cook@nb19'
    name4 = 'spirou@maestria'
    name5 = 'newworlds'
    name6 = 'lam'
    name7 = 'halau'
    # define data directory
    outpaths = dict()
    outpaths[name1] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_rali/'
    outpaths[name2] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_jupiter/'
    outpaths[name3] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_nb19/'
    outpaths[name4] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_maestria/'
    # outpaths[name5] = '/scratch3/lbl/data/minidata_comp/minidata2_07275_newworld/science/'
    # outpaths[name6] = '/scratch3/lbl/data/minidata_comp/minidata2_07275_lam/science/'
    # outpaths[name7] = '/scratch3/lbl/data/minidata_comp/minidata2_07275_halau/rscienceed/'

    # define jobs list
    jobs = []
    # run the wrap for each mini data set
    for name in outpaths:
        job = Process(target=main, args=(name, outpaths[name]))
        jobs.append(job)
        job.start()
    # complete the jobs
    for job in jobs:
        job.join()


# =============================================================================
# End of code
# =============================================================================
