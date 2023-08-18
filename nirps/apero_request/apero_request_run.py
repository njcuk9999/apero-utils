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
    all_requests = general.get_sheet(params, 'response')
    # -------------------------------------------------------------------------
    # create a list of requests
    # -------------------------------------------------------------------------
    requests, valid_requests = general.create_requests(params, all_requests)
    # -------------------------------------------------------------------------
    # loop around requests and run requests
    # -------------------------------------------------------------------------
    for request in requests:
        # request might already be invalid
        if request.valid:
            # run apero get for each request
            request.run_request(params)
        # request might already be invalid (run_request might invalidate)
        if request.valid:
            # copy tar to shared path
            request.copy_tar()
    # -------------------------------------------------------------------------
    # deal with informing users of results
    for request in requests:
        if request.valid:
            # email user
            request.email_success()
        else:
            # email user
            request.email_failure()
    # -------------------------------------------------------------------------
    # update google sheet
    general.update_response_sheet(params, valid_requests)
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
