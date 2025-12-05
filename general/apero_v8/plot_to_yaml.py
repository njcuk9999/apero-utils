#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-04-16 at 09:19

@author: cook
"""
from apero.plotting.plotter import *
import os
import yaml


# =============================================================================
# Define variables
# =============================================================================
PATH = '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/apero/plotting'


# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================



# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # Convert definitions to a serializable dictionary
    yaml_data = {}
    for name, graph in definitions.items():
        entry = {
            'kind': graph.kind,
            'func': graph.func.__name__,  # just the name, not the callable
        }
        if graph.dpi is not None:
            entry['dpi'] = graph.dpi
        if graph.description is not None:
            entry['description'] = graph.description
        if graph.figsize is not None:
            entry['figsize'] = list(graph.figsize)  # convert tuple to list for YAML
        yaml_data[name] = entry

    # Write to a YAML file
    with open(os.path.join(PATH, 'definitions.yaml'), 'w') as fdump:
        yaml.dump(yaml_data, fdump, sort_keys=False)

# =============================================================================
# End of code
# =============================================================================
