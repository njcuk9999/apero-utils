#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 17:03

@author: cook
"""
from typing import Any, Dict, List

import gspread_pandas as gspd
import pandas as pd
from astropy.table import Table

from apero_raw_tests.core import base
from apero_raw_tests.core import io
from apero_raw_tests import tests
from apero_raw_tests.core import misc

# =============================================================================
# Define variables
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def get_obs_dirs(params) -> List[str]:
    # print progress
    msg = '\n' + '*' * 50 + '\n'
    msg += 'Checking for observation directory'
    msg += '\n' + '*' * 50 + '\n'
    misc.log_msg(msg, level='info')
    # get all observation directories
    if params['obsdir'] is None:
        # ask user if this is what they want
        msg = ('\n\n\nNo observation directory set, do you want to '
               'run directories? \n\t[Y]es to continue or [N]o to exit:\t')
        # get user input
        user_input = input(msg)
        print('')
        # if user input is "n" then exit
        if 'n' in user_input.lower():
            return []
        obsdirs = io.get_obs_dirs(params)
    else:
        obsdirs = [params['obsdir']]

    # print that we found some observation directories
    if len(obsdirs) > 0:
        msg = '\tFound {0} observation directories'
        margs = [len(obsdirs)]
        misc.log_msg(msg.format(*margs), level='info')
    else:
        emsg = '\tNo observation directories found'
        misc.log_msg(emsg, level='error')
    # return the observation directories
    return obsdirs


def run_test(params: Dict[str, Any], obsdir: str, test_name: str, it: int,
             num_tests: int, log: bool = True):
    """
    Run a single test

    :param params: parameters dictionary, containing all parameters that can
                   be used
    :param obsdir: str, the observation directory to run the test on
    :param test_name: str, the name of the test to run
    :param it: int, the number of the test
    :param num_tests: int, the total number of tests
    :param log: bool, if True log messages to screen

    :return: bool, True if test passed, False otherwise
    """
    # try to run test
    try:
        # print which test we are running
        msg = '\n\n\tRunning test {0} [{1}/{2}]\n\n'
        margs = [test_name, it + 1, num_tests]
        misc.log_msg(msg.format(*margs), level='')
        output = tests.test_dict[test_name](params, obsdir, log=log)
        # print whether test passed or failed
        if output:
            misc.log_msg('\t\tPASSED', color='green')
        else:
            misc.log_msg('\t\tFAILED', color='red')
    except Exception as e:
        if log:
            raise e
        msg = '\t\tError in test {0} \n\t {1}'
        margs = [test_name, e]
        misc.log_msg(msg.format(*margs), level='warning')
        output = False
    return output


def run_tests(params: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    Run all tests in silent mode and return a dictionary of test values

    :param params: parameters dictionary, containing all parameters that can
                   be used
    :return:
    """
    # get observation directories
    obsdirs = get_obs_dirs(params)
    # storage for test values
    test_values = dict()
    # loop around observation directories
    for obsdir in obsdirs:
        # print message on which observation directory we are processing
        msg = '\n' + '*' * 50 + '\n'
        msg += f'Processing observation directory {obsdir}'
        msg += '\n' + '*' * 50 + '\n'
        misc.log_msg(msg, level='info')
        # storage for this observation directory
        test_values[obsdir] = dict()
        # loop around all tests
        for it, test_name in enumerate(tests.test_dict):
            # run the test
            output = run_test(params, obsdir, test_name, it,
                              len(tests.test_dict), False)
            # add to test values
            test_values[obsdir][test_name] = output
    # return the test values
    return test_values


def run_single_test(params: Dict[str, Any]):
    """
    Run a single test in log mode

    :param params: parameters dictionary, containing all parameters that can
                   be used
        :param obsdir: str, the observation directory to run the test on
    :param test_name:  the name of the test to run
    :return:
    """
    # get observation directories
    obsdir = params['obsdir']
    # deal with no obsdir
    if obsdir is None:
        emsg = ('Must define observation directory for single '
                'test (Either in yaml file or --obsdir')
        raise base.AperoRawTestsError(emsg)
    # -------------------------------------------------------------------------
    # get test_name
    test_name = params['test_name']
    # deal with no test_name
    if test_name is None:
        emsg = ('Must define test_name for single '
                'test (Either in yaml file or --test_name')
        raise base.AperoRawTestsError(emsg)
    # deal with bad test name
    if test_name not in tests.test_dict:
        emsg = 'test_name is not a valid test. \n\nCurrent valid tests:'
        for _test_name in tests.test_dict.keys():
            emsg += f'\n\t- "{_test_name}"'
        raise base.AperoRawTestsError(emsg)
    # -------------------------------------------------------------------------
    # print message on which observation directory we are processing
    msg = '\n' + '*' * 50 + '\n'
    msg += f'Processing single test on observatory directory {obsdir}'
    msg += '\n' + '*' * 50 + '\n'
    misc.log_msg(msg, level='info')
    # run single test
    _ = run_test(params, obsdir, test_name, 0, 1, True)


# =============================================================================
# Deal with uploading tests
# =============================================================================
def add_to_sheet(params: Dict[str, Any], dataframe: pd.DataFrame):
    """
    Add all listed astrometrics objects to the sheet

    :param params:
    :param astro_objs:
    :return:
    """
    # deal with no rows in dataframe
    if len(dataframe) == 0:
        msg = 'No rows to add to dataframe'
        misc.log_msg(msg, level='warning')
        return
    # add gspread directory and auth files
    io.gsp_setup()
    # define the sheet id and sheet name (pending)
    sheet_id = params['sheet id']
    sheet_name = params['sheet name']
    # load google sheet instance
    google_sheet = gspd.spread.Spread(sheet_id)
    # convert google sheet to pandas dataframe
    current_dataframe = google_sheet.sheet_to_df(index=0, sheet=sheet_name)
    # -------------------------------------------------------------------------
    # append empty rows to dataframe
    current_dataframe = pd.concat([current_dataframe, dataframe],
                                  ignore_index=True)
    # convert date column to date time type
    current_dataframe['date'] = pd.to_datetime(current_dataframe['date'])
    # sort the dataframe by date column in descending order
    current_dataframe = current_dataframe.sort_values(by='date',
                                                      ascending=False)
    # Group the dataframe by obsdir and get hte index of the maximum date
    # (i.e. the most recent date)
    idx = current_dataframe.groupby('obsdir')['date'].idxmax()
    # filter the dataframe using the obtained index values
    current_dataframe = current_dataframe.loc[idx]
    # reset the index of the filtered dataframe
    current_dataframe = current_dataframe.reset_index(drop=True)
    # -------------------------------------------------------------------------
    # resort values ascending in date
    current_dataframe = current_dataframe.sort_values(by='date',
                                                      ascending=True)
    # -------------------------------------------------------------------------
    # print progress
    msg = 'Pushing all rows to google-sheet'
    misc.log_msg(msg, level='info')
    # push dataframe back to server
    google_sheet.df_to_sheet(current_dataframe, index=False, replace=True)
    # print progress
    msg = 'All rows added to google-sheet'
    misc.log_msg(msg, level='info')


def upload_tests(params: Dict[str, Any], results: Dict[str, Dict[str, Any]]):
    """
    Upload test results to google sheet

    :param params: parameters dictionary, containing all parameters that can
                   be used
    :param results: a dictionary of dictionarys where each element is a
                    observation night and each sub-element is a test result
                    (key = column name, value = True or False)
    :return:
    """
    # deal with no results
    if len(results) == 0:
        msg = 'No rows to add to google sheet'
        misc.log_msg(msg, level='warning')
        return
    # print message on which observation directory we are processing
    msg = '\n' + '*' * 50 + '\n'
    msg += f'Uploading to google sheet'
    msg += '\n' + '*' * 50 + '\n'
    misc.log_msg(msg, level='info')

    # get time now
    time_now = base.AstropyTime.now().iso
    # store columns
    columns = ['obsdir', 'date']
    # turn results dictionary into a single dictionary
    results_dict = dict()
    results_dict['obsdir'] = []
    results_dict['date'] = []
    for obsdir in results:
        results_dict['obsdir'] += [obsdir]
        results_dict['date'] += [time_now]
        # loop around tests
        for test in results[obsdir]:
            if test not in results_dict:
                columns.append(test)
                results_dict[test] = []
            results_dict[test] += [results[obsdir][test]]
    # -------------------------------------------------------------------------
    # turn results dictionary into dataframe
    dataframe = pd.DataFrame(results_dict, columns=columns)
    # -------------------------------------------------------------------------
    # add to sheet
    add_to_sheet(params, dataframe)


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
