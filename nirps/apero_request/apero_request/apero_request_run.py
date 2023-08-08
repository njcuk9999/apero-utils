#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on {DATE}

@author: cook
"""
from apero_request.core import base
from apero_request.core import general
from apero_request.core import misc

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
def main():
    # load params
    params = misc.load_params()
    # -------------------------------------------------------------------------
    # get current list of requests
    # -------------------------------------------------------------------------
    current_requests = general.get_sheet(params, 'response')
    # -------------------------------------------------------------------------
    # create a list of requests
    # -------------------------------------------------------------------------
    requests, mask = general.create_requests(params, current_requests)
    # -------------------------------------------------------------------------
    # loop around requests
    for request in requests:
        # create directory for each request
        request.create_dir()
        # run apero get for each request
        request.run()
        # tar directory for each request
        request.tar()
        # copy tar to shared path
        request.copytar()
        # remove directory for each request
        request.remove_dir()
    # -------------------------------------------------------------------------
    # update google sheet
    general.update_response_sheet(params, current_requests[mask])
    # -------------------------------------------------------------------------
    # finish with an end message
    misc.end_msg()


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # run main function
    main()

# =============================================================================
# End of code
# =============================================================================
