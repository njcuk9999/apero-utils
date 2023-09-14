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
    """
    Wrapper around main code to catch errors (and exit nicely)
    i.e. deal with the lock file
    :return:
    """
    # lets lock (or not run if locked)
    not_running = general.lock(stop=True)
    # if running then we exit
    if not not_running:
        misc.log_msg('Request already running - exiting')
    # try to run main code
    try:
        __main__()
        # unlock the code
        general.unlock()
    except Exception as e:
        # unlock the code
        general.unlock()
        # raise exception
        raise e


def __main__():
    # load params
    params = misc.load_params()
    # get splash
    misc.splash('APERO requests')
    # -------------------------------------------------------------------------
    # get current list of requests
    # -------------------------------------------------------------------------
    misc.log_msg('Gettings requests')
    all_requests = general.get_sheet(params, 'response')
    # -------------------------------------------------------------------------
    # create a list of requests
    # -------------------------------------------------------------------------
    misc.log_msg('Creating requests')
    requests = general.create_requests(params, all_requests)
    # -------------------------------------------------------------------------
    # Remove all invalid tar files + find already processed requests
    # -------------------------------------------------------------------------
    misc.log_msg('Analysing local files')
    requests = general.remove_invalid_tars(params, requests)
    # -------------------------------------------------------------------------
    # loop around requests and run requests
    # -------------------------------------------------------------------------
    for r_it, request in enumerate(requests):
        # print where we are up to
        msg = 'Generating request {0} / {1}'
        misc.log_msg(msg.format(r_it +1, len(requests)))
        # request might already be invalid
        if request.valid:
            # run apero get for each request
            request.run_request(params)
        else:
            misc.log_msg('\tSkipping request')
    # -------------------------------------------------------------------------
    # deal with informing users of results
    for r_it, request in enumerate(requests):
        # print where we are up to
        msg = 'Emailing user for request {0} / {1}'
        misc.log_msg(msg.format(r_it +1, len(requests)))
        # email a success
        if request.valid and not request.exists:
            # email user
            misc.log_msg('\tEmailing success')
            request.email_success(params, r_it)
        elif not request.exists:
            # email user
            misc.log_msg('\tEmailing failure')
            request.email_failure(params, r_it)
        else:
            msg = 'Request {0} already exists:'
            misc.log_msg(msg.format(r_it))
            print(request)
    # -------------------------------------------------------------------------
    # copy index.html file over
    general.copy_index(params)
    # rsync data to server
    general.sync_data(params)
    # -------------------------------------------------------------------------
    # recreate the dataframe from request
    valid_requests, all_requests = general.create_dataframe(params, requests)
    # update response google sheet
    general.update_response_sheet(params, valid_requests)
    # update archive google sheet
    general.update_archive_sheet(params, all_requests)
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
