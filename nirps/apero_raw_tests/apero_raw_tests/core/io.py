#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 14:51

@author: cook
"""
from typing import Dict, Any

import yaml


# =============================================================================
# Define variables
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def read_yaml(yaml_filename: str) -> Dict[str, Any]:
    """
    Read the yaml file and add to settings

    :param yaml_filename: str, yaml file name
    :param settings: dict, settings dictionary

    :return: dict, updated settings dictionary
    """
    # read the yaml file
    with open(yaml_filename, 'r') as f:
        yaml_data = yaml.load(f, Loader=yaml.FullLoader)
    # add a profiles sub-dictionary
    settings = dict()
    # loop around yaml data
    for key, value in yaml_data.items():
        # if key is in settings
        settings[key] = value
    # return settings
    return settings


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
