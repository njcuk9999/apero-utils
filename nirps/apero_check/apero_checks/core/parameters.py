#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 16:26

@author: cook
"""
from typing import Any, Optional

from apero_checks.core import base

# =============================================================================
# Define variables
# =============================================================================
# version, date, author
__VERSION__ = base.__VERSION__
__DATE__ = base.__DATE__
__AUTHOR__ = base.__AUTHOR__


# =============================================================================
# Define functions
# =============================================================================
class Const:
    def __init__(self, key, value, not_none=False, dtype=None,
                 report: bool = False, path: Optional[str] = None):
        self.key = key
        self.value = value
        self.not_none = not_none
        self.dtype = dtype
        self.report = report
        self.path = path

    def check(self, value: Optional[Any] = None, source: str = 'NULL'):
        # deal with internal check
        if value is None:
            value = self.value
            source = 'PARAM'
        # check if not none
        if self.not_none and value is None:
            emsg = 'Const {0} must be set in yaml'
            eargs = [self.key]
            raise base.AperoChecksError(emsg.format(*eargs))
        # force dtype if set
        if self.dtype is not None:
            try:
                return self.dtype(value), source
            except Exception as e:
                emsg = 'Const {0} must be of type {1} (error: {2})'
                eargs = [self.key, self.dtype, e]
                raise base.AperoChecksError(emsg.format(*eargs))
        # return value
        return value, source


# =============================================================================
# Define default variables to add
# =============================================================================
# dictionary to store parameters
parameters = dict()

# directory where raw data is downloaded
parameters['raw dir'] = Const('raw dir', None, not_none=True, dtype=str,
                              path='general.raw dir')

# the apero profile to use (controlled sym link raw data directory and all
#    apero reduced data directories)
parameters['apero profile'] = Const('apero profile', None, not_none=True,
                                    dtype=str, path='general.apero profile')

# the apero install directory
parameters['apero install'] = Const('apero install', None, not_none=True,
                                    dtype=str,
                                    path='general.apero install')

# the lbl data path
parameters['lbl path'] = Const('lbl path', None, not_none=True, dtype=str,
                               path='general.lbl path')

# observation directory (i.e. night name)
parameters['obsdir'] = Const('obs dir', None, dtype=str, report=True)

# the test to run (if running a single test)
parameters['test_name'] = Const('testname', None, dtype=str, report=True)

# the processing arguments
processing_dict = dict()
processing_dict['run file'] = None
parameters['processing'] = Const('processing', processing_dict, dtype=dict,
                                 path='processing')

# define the sheet id for the google sheet
parameters['raw sheet id'] = Const('raw sheet id',
                                   '1zvU_XFA1ZOJE111qZKiav7v6ptYveWpDkMHjdhiN06M',
                                   dtype=str, path='check.raw sheet id')
# define the sheet name for the google sheet
parameters['raw sheet name'] = Const('raw sheet name', None, dtype=str,
                                     not_none=True,
                                     path='check.raw sheet name')

# define the sheet id for the raw checks
parameters['raw checks sheet id'] = Const('raw checks sheet id',
                                          '1741354113', # HA: 1017212188
                                          dtype=str,
                                          path='check.raw checks sheet id')

# define the sheet id for the google sheet
parameters['red sheet id'] = Const('red sheet id',
                                   '1zvU_XFA1ZOJE111qZKiav7v6ptYveWpDkMHjdhiN06M',
                                   dtype=str, path='check.red sheet id')

# define the sheet name for the google sheet
parameters['red sheet name'] = Const('red sheet name', None, dtype=str,
                                     not_none=True,
                                     path='check.red sheet name')

# define the sheet id for the raw checks
parameters['red checks sheet id'] = Const('raw checks sheet id',
                                          '541545005', # HA: 279548254
                                          dtype=str,
                                          path='check.raw checks sheet id')

# define the sheet id for the google sheet
parameters['over sheet id'] = Const('over sheet id',
                                    '1zvU_XFA1ZOJE111qZKiav7v6ptYveWpDkMHjdhiN06M',
                                    dtype=str, path='check.over sheet id')

# define the sheet name for the google sheet
parameters['over sheet name'] = Const('over sheet name',
                                      'OVERRIDES', path='check.over sheet name',
                                      dtype=str)

# define the sheet id for the check comments
parameters['comments sheet id'] = Const('comments sheet id',
                                        '68390787', # HA: 346030560
                                        dtype=str,
                                        path='check.comments sheet id')

# allow disabling dependencies check
parameters['disable dependencies check'] = Const('disable dependencies check',
                                                 False,
                                                 dtype=bool,
                                                 path='check.disable dependencies check')

# define the critial checks csv file (set to None to skip)
parameters['critical csv file'] = Const('critical csv file',
                                        '/nirps_raw/nirps/nirps/critical-checks-output/check_status.csv',
                                        dtype=str,
                                        path='check.critical csv file')

# define the critical checks description file (set to None to skip)
parameters['critical desc file'] = Const('critical desc file',
                                         '/cosmos99/nirps/git-bin/nirpsdl/assets/critical_checks_description.csv',
                                         dtype=str,
                                         path='check.critical desc file')

# define the system disk path
parameters['system disk path'] = Const('system disk path', '/cosmos99',
                                        dtype=str,
                                        path='check.system disk path')

# define the system email server
parameters['system email server'] = Const('system email server', 'localhost',
                                          dtype=str,
                                          path='check.system email server')

# set the allocation sheet id
parameters['allocation sheet id'] = Const('allocation sheet id',
                                          '1s116aabnMH0zJ5YbXrGBWIVP17lmYXl6tYzLteoSEAc',
                                          dtype=str,
                                          path='check.allocation sheet id')

# set the allocation sheet name
parameters['allocation sheet name'] = Const('allocation sheet name',
                                            '0',
                                            dtype=str,
                                            path='check.allocation sheet name')



# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # print 'Hello World!'
    print("Hello World!")

# =============================================================================
# End of code
# =============================================================================
