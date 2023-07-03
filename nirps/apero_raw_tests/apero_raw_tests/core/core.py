#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 14:51

@author: cook
"""
from typing import Dict, Any
from apero_raw_tests.core import io


# =============================================================================
# Define variables
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def load_params(yaml_file):

    params = dict()
    # load default parameters


    # load from yaml file
    yaml_params = io.read_yaml(yaml_file)

    # push into params


    # return params
    return params


def run_tests(params: Dict[str, Any]):
    pass



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
