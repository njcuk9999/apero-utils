#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 16:26

@author: cook
"""
from typing import Any, Optional

from apero_raw_tests.core import base

# =============================================================================
# Define variables
# =============================================================================
__VERSION__ = base.__VERSION__
__AUTHOR__ = base.__AUTHOR__

# =============================================================================
# Define functions
# =============================================================================
class Const:
    def __init__(self, key, value, not_none=False, dtype=None):
        self.key = key
        self.value = value
        self.not_none = not_none
        self.dtype = dtype

    def check(self, value: Optional[Any] = None):
        # deal with internal check
        if value is None:
            value = self.value
        # check if not none
        if self.not_none and value is None:
            emsg = 'Const {0} must be set in yaml'
            eargs = [self.key]
            raise base.AperoRawTestsError(emsg.format(*eargs))
        # force dtype if set
        if self.dtype is not None:
            try:
                return self.dtype(value)
            except Exception as e:
                emsg = 'Const {0} must be of type {1} (error: {2})'
                eargs = [self.key, self.dtype, e]
                raise base.AperoRawTestsError(emsg.format(*eargs))
        # return value
        return value


# =============================================================================
# Define default variables to add
# =============================================================================
# dictionary to store parameters
parameters = dict()

# directory where raw data is downloaded
parameters['raw dir'] = Const('raw dir', None, not_none=True, dtype=str)

# observation directory (i.e. night name)
parameters['obsdir'] = Const('obs dir', None, not_none=True, dtype=str)

# the test to run (if running a single test)
parameters['testname'] = Const('testname', None, dtype=str)



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
