#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 17:03

@author: cook
"""
import os
import time
from typing import Any, Dict, List

import numpy as np
import gspread_pandas as gspd
import pandas as pd

from apero_checks import raw_tests
from apero_checks import red_tests
from apero_checks.core import base
from apero_checks.core import io
from apero_checks.core import misc

# =============================================================================
# Define variables
# =============================================================================
# version, date, author
__VERSION__ = base.__VERSION__
__DATE__ = base.__DATE__
__AUTHOR__ = base.__AUTHOR__
# define lock directory
LOCK_DIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
# define maximum lock out time in minutes
MAX_COUNT = 30
# -----------------------------------------------------------------------------


# =============================================================================
# Define functions
# =============================================================================
def get_obs_dirs(params) -> List[str]:
    # print progress
    msg = '*' * 50
    msg += '\nChecking for observation directory'
    msg += '\n' + '*' * 50
    misc.log_msg(msg, level='info')
    # get all observation directories
    if params['obsdir'] is None:
        # ask user if this is what they want
        msg = ('\n\nNo observation directory set, do you want to '
               'run directories? \n\t[Y]es to continue or [N]o to exit:\t')
        # get user input
        user_input = input(msg)
        print('')
        # if user input is "n" then exit
        if 'y' not in user_input.lower():
            return []
        obsdirs = io.get_obs_dirs(params['raw dir'])
    elif isinstance(params['obsdir'], list):
        obsdirs = params['obsdir']
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
             num_tests: int, log: bool = True, test_type: str = 'raw'):
    """
    Run a single test

    :param params: parameters dictionary, containing all parameters that can
                   be used
    :param obsdir: str, the observation directory to run the test on
    :param test_name: str, the name of the test to run
    :param it: int, the number of the test
    :param num_tests: int, the total number of tests
    :param log: bool, if True log messages to screen
    :param test_type: str, either "raw" or "red" (the type of tests to run)

    :return: bool, True if test passed, False otherwise
    """
    # try to run test
    try:

        if test_type == 'raw':
            # print which test we are running
            msg = '\tRunning raw test {0} [{1}/{2}]'
            margs = [test_name, it + 1, num_tests]
            misc.log_msg(msg.format(*margs), level='test')
            # run raw tests
            output = raw_tests.test_dict[test_name](params, obsdir, log=log)
        elif test_type == 'red':
            # print which test we are running
            msg = '\tRunning red test {0} [{1}/{2}]'
            margs = [test_name, it + 1, num_tests]
            misc.log_msg(msg.format(*margs), level='test')
            # run red tests
            output = red_tests.test_dict[test_name](params, obsdir, log=log)
        else:
            emsg = 'RUN_TEST error: test_type must be set to "raw" or "red"'
            raise base.AperoChecksError(emsg)
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


def run_tests(params: Dict[str, Any],
              test_type: str) -> Dict[str, Dict[str, Any]]:
    """
    Run all tests in silent mode and return a dictionary of test values

    :param params: parameters dictionary, containing all parameters that can
                   be used
    :param test_type: str, either "raw" or "red" (the type of tests to run)

    :return:
    """
    # get observation directories
    obsdirs = get_obs_dirs(params)
    # storage for test values
    test_values = dict()
    # -------------------------------------------------------------------------
    if params['test_filter'] is not None:
        current_dataframe = get_current_dataframe(params, test_type)
        # index by obsdir
        obsdir_dataframe = current_dataframe.set_index('obsdir')
    else:
        obsdir_dataframe = None
    # -------------------------------------------------------------------------
    # loop around observation directories
    for obsdir in obsdirs:
        # print message on which observation directory we are processing
        msg = '*' * 50
        msg += f'\nProcessing observation directory {obsdir}'
        msg += '\n' + '*' * 50
        misc.log_msg(msg, level='info')
        # get a list of tests
        if test_type == 'raw':
            test_list = raw_tests.test_dict
        elif test_type == 'red':
            test_list = red_tests.test_dict
        else:
            emsg = 'RUN_TESTS error: test_type must be set to "raw" or "red"'
            raise base.AperoChecksError(emsg)
        # storage for this observation directory
        test_values[obsdir] = dict()
        # loop around all tests
        for it, test_name in enumerate(test_list):
            # deal with skipping tests
            if params['test_filter'] is not None:
                if test_name not in params['test_filter']:
                    if obsdir not in obsdir_dataframe.index:
                        value = ''
                    elif test_name in obsdir_dataframe:
                        value = obsdir_dataframe.loc[obsdir, test_name]
                    else:
                        value = ''
                    # push into test_values
                    test_values[obsdir][test_name] = value
                    # skip actually running the test
                    continue
            # run the test
            output = run_test(params, obsdir, test_name, it=it,
                              num_tests=len(test_list),
                              log=False, test_type=test_type)
            # add to test values
            test_values[obsdir][test_name] = output
    # return the test values
    return test_values


def run_single_test(params: Dict[str, Any], test_type: str):
    """
    Run a single test in log mode

    :param params: parameters dictionary, containing all parameters that can
                   be used
    :param test_type: str, either "raw" or "red" (the type of tests to run)

    :return:
    """
    # get observation directories
    obsdirs = get_obs_dirs(params)
    # deal with no obsdir
    if obsdirs is None:
        emsg = ('Must define observation directory for single '
                'test (Either in yaml file or --obsdir')
        raise base.AperoChecksError(emsg)
    # -------------------------------------------------------------------------
    # get a list of tests
    if test_type == 'raw':
        test_list = raw_tests.test_dict
    elif test_type == 'red':
        test_list = red_tests.test_dict
    else:
        emsg = 'RUN_TESTS error: test_type must be set to "raw" or "red"'
        raise base.AperoChecksError(emsg)
    # get test_name
    test_name = params['test_name']
    # deal with no test_name
    if test_name is None:
        emsg = ('Must define test_name for single '
                'test (Either in yaml file or --test_name')
        raise base.AperoChecksError(emsg)
    # deal with bad test name
    if test_name not in test_list:
        emsg = 'test_name is not a valid test. \nCurrent valid tests:'
        for _test_name in test_list.keys():
            emsg += f'\n\t- "{_test_name}"'
        raise base.AperoChecksError(emsg)
    # -------------------------------------------------------------------------
    for obsdir in obsdirs:
        # print message on which observation directory we are processing
        msg = '*' * 50
        msg += f'\nProcessing single test on observatory directory {obsdir}'
        msg += '\n' + '*' * 50
        misc.log_msg(msg, level='info')
        # run single test
        _ = run_test(params, obsdir, test_name, it=0, num_tests=1, log=True,
                     test_type=test_type)
    # -------------------------------------------------------------------------
    # print a note that the single test does not update the database
    msg = ('*' * 50 + '\nPlease note\n' + '*' * 50 +
           '\n\tSingle test mode (--test) does NOT update the '
           'database. '
           '\n\tPlease run without the --test argument to update the '
           'database. '
           '\n\t--testfilter can be used to run only some tests.')
    misc.log_msg(msg, level='info', color='yellow')


# =============================================================================
# Deal with uploading tests
# =============================================================================
def get_current_dataframe(params, test_type='raw'):
    # add gspread directory and auth files
    io.gsp_setup()
    # define the sheet id and sheet name (pending)
    if test_type == 'raw':
        sheet_id = params['raw sheet id']
        sheet_name = params['raw sheet name']
    elif test_type == 'red':
        sheet_id = params['red sheet id']
        sheet_name = params['red sheet name']
    elif test_type == 'override':
        sheet_id = params['over sheet id']
        sheet_name = params['over sheet name']
    else:
        emsg = 'ADD_TO_SHEET error: test_type must be set to "raw" or "red"'
        raise base.AperoChecksError(emsg)
    # load google sheet instance
    google_sheet = gspd.spread.Spread(sheet_id)
    # convert google sheet to pandas dataframe
    return google_sheet.sheet_to_df(index=0, sheet=sheet_name)


def add_to_sheet(params: Dict[str, Any], dataframe: pd.DataFrame,
                 test_type='raw'):
    """
    Add all listed astrometrics objects to the sheet

    :param params: parameters dictionary, containing all parameters that can
                   used
    :param dataframe: pandas dataframe, the dataframe to add to the sheet
    :param test_type: str, either "raw" or "red" (the type of tests to run)

    :return: None, adds to google sheet
    """
    # deal with no rows in dataframe
    if len(dataframe) == 0:
        msg = 'No rows to add to dataframe'
        misc.log_msg(msg, level='warning')
        return
    # add gspread directory and auth files
    io.gsp_setup()
    # define the sheet id and sheet name (pending)
    if test_type == 'raw':
        sheet_id = params['raw sheet id']
        sheet_name = params['raw sheet name']
    elif test_type == 'red':
        sheet_id = params['red sheet id']
        sheet_name = params['red sheet name']
    elif test_type == 'override':
        sheet_id = params['over sheet id']
        sheet_name = params['over sheet name']
    else:
        emsg = ('ADD_TO_SHEET error: test_type must be set to "raw" or "red" '
                'or "override"')
        raise base.AperoChecksError(emsg)
    # -------------------------------------------------------------------------
    # make a local name for the file
    local_path = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    local_file = os.path.join(local_path, 'local', f'{sheet_id}_{sheet_name}.csv')
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
    # deal with grouping by date
    if test_type in ['red', 'raw']:
        # Group the dataframe by obsdir and get hte index of the maximum date
        # (i.e. the most recent date)
        idx = current_dataframe.groupby('obsdir')['date'].idxmax()
        # filter the dataframe using the obtained index values
        current_dataframe = current_dataframe.loc[idx]
        # reset the index of the filtered dataframe
        current_dataframe.reset_index(drop=True)
    elif test_type in 'override':
        # Drop dulicates based on
        dkwargs = dict(subset=['obsdir', 'test_type', 'test_name'],
                       keep='first', inplace=True)
        current_dataframe.drop_duplicates(**dkwargs)
        # reset the index of the filtered dataframe
        current_dataframe.reset_index(drop=True, inplace=True)
    # -------------------------------------------------------------------------
    # resort values ascending in date
    current_dataframe = current_dataframe.sort_values(by='obsdir',
                                                      ascending=False)
    # -------------------------------------------------------------------------
    # load last table in local directory
    if os.path.exists(local_file):
        last_dataframe = pd.read_csv(local_file)

        # check that the new dataframe isn't shorter than the old one
        if len(current_dataframe) < len(last_dataframe):
            # push dataframe back to server
            google_sheet.df_to_sheet(last_dataframe, index=False,
                                     replace=True)
            emsg = (f'Sheet {sheet_name} ({sheet_id}) has got shorter - '
                    f'something went wrong. Please delete {last_dataframe} and '
                    f'try again - note we are resetting the online version to '
                    f'this last version.')
            raise base.AperoChecksError(emsg)
    # if local file still exists remove it
    if os.path.exists(local_file):
        os.remove(local_file)
    # print progress
    msg = 'Saving local backup ({0})'.format(local_file)
    misc.log_msg(msg, level='info')
    # save table to local directory
    current_dataframe.to_csv(local_file, index=False)
    # -------------------------------------------------------------------------
    # print progress
    msg = 'Pushing all rows to google-sheet ({0})'.format(sheet_name)
    misc.log_msg(msg, level='info')
    # push dataframe back to server
    google_sheet.df_to_sheet(current_dataframe, index=False, replace=True)
    # print progress
    msg = 'All rows added to google-sheet ({0})'.format(sheet_name)
    misc.log_msg(msg, level='info')


def upload_tests(params: Dict[str, Any], results: Dict[str, Dict[str, Any]],
                 test_type: str = 'raw'):
    """
    Upload test results to google sheet

    :param params: parameters dictionary, containing all parameters that can
                   be used
    :param results: a dictionary of dictionarys where each element is a
                    observation night and each sub-element is a test result
                    (key = column name, value = True or False)
    :param test_type: str, either "raw" or "red" (the type of tests to run)

    :return:
    """
    # deal with no results
    if len(results) == 0:
        msg = 'No rows to add to google sheet'
        misc.log_msg(msg, level='warning')
        return
    # print message on which observation directory we are processing
    msg = '*' * 50
    msg += f'\nUploading to google sheet'
    msg += '\n' + '*' * 50
    misc.log_msg(msg, level='info')
    # get time now
    time_now = base.AstropyTime.now().iso
    # remove the decimal seconds
    time_now = time_now.split('.')[0]
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
    # lock codes
    lock()
    # -------------------------------------------------------------------------
    # turn results dictionary into dataframe
    dataframe = pd.DataFrame(results_dict, columns=columns)
    # -------------------------------------------------------------------------
    # add to sheet
    try:
        add_to_sheet(params, dataframe, test_type=test_type)
        # ---------------------------------------------------------------------
    finally:
        # unlock codes
        unlock()


def store_overrides(params: Dict[str, Any],
                    overrides: Dict[str, Dict[str, Any]],
                    test_type: str = 'raw'):
    # deal with no results
    if len(overrides) == 0:
        msg = 'No override rows to add to google sheet'
        misc.log_msg(msg, level='warning')
        return

    # define the sheet id and sheet name (pending)
    if test_type == 'raw':
        sheet_name = params['raw sheet name']
    elif test_type == 'red':
        sheet_name = params['red sheet name']
    else:
        return

    # print message on which observation directory we are processing
    msg = '*' * 50
    msg += f'\nUploading to override google sheet'
    msg += '\n' + '*' * 50
    misc.log_msg(msg, level='info')
    # get time now
    time_now = base.AstropyTime.now().iso
    # remove the decimal seconds
    time_now = time_now.split('.')[0]
    # store columns
    columns = ['obsdir', 'date', 'test_type', 'test_name', 'test_value',
               'who', 'comment']
    # turn results dictionary into a single dictionary
    results_dict = dict()
    for col in columns:
        results_dict[col] = []
    # loop around all overrides and populate results_dict
    for obsdir in overrides:
        # loop around tests
        for test_name in overrides[obsdir]:
            results_dict['obsdir'] += [obsdir]
            results_dict['date'] += [time_now]
            results_dict['test_type'] += [sheet_name]
            results_dict['test_name'] += [test_name]
            # get the values from the tuple
            value, name, reason = overrides[obsdir][test_name]
            # push into results
            results_dict['test_value'] += [value]
            results_dict['who'] += [name]
            results_dict['comment'] += [reason]
    # -------------------------------------------------------------------------
    # lock codes
    lock()
    # -------------------------------------------------------------------------
    # turn results dictionary into dataframe
    dataframe = pd.DataFrame(results_dict, columns=columns)
    # -------------------------------------------------------------------------
    # add to sheet
    try:
        add_to_sheet(params, dataframe, test_type='override')
        # ---------------------------------------------------------------------
        # unlock codes
        unlock()
    except Exception as e:
        # ---------------------------------------------------------------------
        # unlock codes
        unlock()
        # raise exception
        raise e


def check_override(params: Dict[str, Any],
                   results: Dict[str, Dict[str, Any]],
                   test_type: str = 'raw'):
    # get any overrides we have
    override_dataframe = get_current_dataframe(params, 'override')
    # define the sheet id and sheet name (pending)
    if test_type == 'raw':
        sheet_name = params['raw sheet name']
    elif test_type == 'red':
        sheet_name = params['red sheet name']
    else:
        emsg = 'ADD_TO_SHEET error: test_type must be set to "raw" or "red"'
        raise base.AperoChecksError(emsg)
    # filter dataframe by sheet_name
    type_mask = override_dataframe['test_type'] == sheet_name
    override_dataframe = override_dataframe[type_mask]
    # get the override obsdirs
    override_obsdirs = np.array(override_dataframe['obsdir'])
    # get the test names
    override_test_names = np.array(override_dataframe['test_name'])
    # get the test values
    override_test_values = np.array(override_dataframe['test_value'])
    # keep a counter of how many values we updated
    counter = 0
    # loop around all rows we are overriding
    for row in range(len(override_obsdirs)):
        # get the obsdir we are overriding
        override_obsdir = override_obsdirs[row]
        # get the test name
        override_test_name = override_test_names[row]
        # get the test value
        override_test_value = override_test_values[row]
        # if the obsdir is not in the results we skip
        if override_obsdir not in results:
            continue
        # if the test_name is not in the results we skip
        if override_test_name not in results[override_obsdir]:
            continue
        # update the test value - if it has changed
        if results[override_obsdir][override_test_name] != override_test_value:
            results[override_obsdir][override_test_name] = override_test_value
            counter += 1
    # log that we updated "counter" values
    msg = 'Updated {0} override values'.format(counter)
    misc.log_msg(msg, level='info')
    # return the updated results
    return results


def lock(stop: bool = False) -> bool:
    """
    Lock code while something is happening
    (avoids two codes doing this at same time)

    :return:
    """
    # get lock file
    lock_file = os.path.join(LOCK_DIR, 'lock.lock')
    # counter
    counter = 0
    # wait a small random amount of time (to avoid two codes doing this
    # at the same time)
    random_time = 0.1 + 0.05 * np.random.uniform(-1, 1)
    time.sleep(random_time)
    # if lock file exists
    while os.path.exists(lock_file):
        # if stop we assume we return straight away with a False
        if stop:
            return False
        # deal with a max count
        if counter > (MAX_COUNT * 60):
            emsg = f'Lock file timeout. Please remove: {lock_file} manually'
            raise base.AperoChecksError(emsg)
        # every 30 counts print a lock message
        if counter % 30 == 0:
            print(f'LOCKED [{counter}s]: Waiting for lock file to be removed '
                  f'by another process')
        # check again in 1 seconds - hopefully other process if fixed by then
        time.sleep(1)
        counter += 1
    # if file doesn't exist create and return
    if not os.path.exists(lock_file):
        # create lock file
        with open(lock_file, 'w') as f:
            timenow = time.time()
            f.write(f'LOCKED: {timenow}')
        return True


def unlock():
    """
    Unlock code (remove lock.lock file)
    :return:
    """
    # get lock file
    lock_file = os.path.join(LOCK_DIR, 'lock.lock')
    # if lock file exists
    if os.path.exists(lock_file):
        # delete the lock file
        os.remove(lock_file)


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
