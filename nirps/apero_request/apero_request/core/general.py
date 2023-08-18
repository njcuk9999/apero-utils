#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on {DATE}

@author: cook
"""
import os
import sys
from typing import Any, Dict, List, Tuple, Union

import gspread_pandas as gspd
import numpy as np
import pandas as pd
from astropy.time import Time

from apero_request.core import base
from apero_request.core import io
from apero_request.core import misc

# =============================================================================
# Define variables
# =============================================================================
# version, date, author
__VERSION__ = base.__VERSION__
__DATE__ = base.__DATE__
__AUTHOR__ = base.__AUTHOR__


# =============================================================================
# Define classes
# =============================================================================
class Request:
    def __init__(self, timestamp: Union[str, pd.Timestamp],
                 email: str, drsobjn: str, dprtype: str, mode: str, fibers: str,
                 drsoutid: str, passkey: str,
                 start_date: Union[str, pd.Timestamp],
                 end_date: Union[str, pd.Timestamp]):
        """
        Create a request instance

        :param timestamp:
        :param email:
        :param drsobjn:
        :param dprtype:
        :param mode:
        :param fiber:
        :param drsoutid:
        :param passkey:
        :param start_date:
        :param end_date:
        """
        self.timestamp = pd.to_datetime(timestamp)
        self.email = email
        self.drsobjn = np.char.array(drsobjn.split(',')).strip()
        self.dprtype = np.char.array(dprtype.split(',')).strip()
        self.mode = mode
        self.fibers = fibers
        self.drsoutid = np.char.array(drsoutid.split(',')).strip()
        self.passkey = passkey

        if start_date in [None, 'None', 'Null', '']:
            self.start_date = None
        else:
            self.start_date = pd.to_datetime(start_date, format="%d/%m/%Y")
        if end_date in [None, 'None', 'Null', '']:
            self.end_date = None
        else:
            self.end_date = pd.to_datetime(end_date, format="%d/%m/%Y")

        # other attributes
        self.valid: bool = False
        self.reason: str = ''
        self.allowed_runids: set = set()
        self.allowed_profiles: set = set()

    @staticmethod
    def from_pandas_row(row: Any) -> 'Request':
        """
        Create a request from a pandas row

        :param row: pandas row, the row to convert to a request

        :return: Request, the request instance
        """
        # get values from pandas row
        timestamp = row['Timestamp']
        email = row['Email Address']
        drsobjn = row['DRSOBJN']
        dprtype = row['DPRTYPE']
        mode = row['APERO_MODE']
        fiber = row['FIBER']
        drsoutid = row['DRSOUTID']
        passkey = row['Pass key']
        start_date = row['START_DATE']
        end_date = row['END_DATE']
        # make the request
        request = Request(timestamp, email, drsobjn, dprtype, mode, fiber,
                          drsoutid, passkey, start_date, end_date)
        # return the request
        return request

    def run_request(self, params: dict):
        """

        :return:
        """
        # select profile parameters
        try:
            profile_params = self._get_profile_params(params)
        except Exception as e:
            self.valid = False
            self.reason = f'Could not get profile parameters with error: {e}'
            return
        # update apero profile
        try:
            _ = update_apero_profile(profile_params)
        except Exception as e:
            self.valid = False
            self.reason = f'Could not update apero profile with error: {e}'
            return
        # ---------------------------------------------------------------------
        # sort out arguments for apero get
        # ---------------------------------------------------------------------
        try:
            objnames = ','.join(list(self.drsobjn))
            dprtypes = ','.join(list(self.dprtype))
            outtypes = ','.join(list(self.drsoutid))
            outpath = params['LOCAL_DIR']
            fibers = get_fibers(self.fibers)
            runids = get_runids(self.allowed_runids)
            test = params.get('test', False)
            if self.start_date is not None:
                start = Time(self.start_date).iso
            else:
                start = 'None'
            if self.end_date is not None:
                end = Time(self.end_date).iso
            else:
                end = 'None'
        except Exception as e:
            self.valid = False
            self.reason = f'Could not set up apero get parameters. Error {e}'
            return
        # ---------------------------------------------------------------------
        # try to run apero get
        # ---------------------------------------------------------------------
        try:
            # need to import apero_get (for this profile)
            from apero.tools.recipes.bin import apero_get
            # run apero get to make the objects dir in apero dir
            apero_get.main(objnames=objnames, dprtypes=dprtypes,
                           outtypes=outtypes, outpath=outpath,
                           fibers=fibers, runids=runids,
                           symlinks=False, tar=True,
                           test=test, since=start, latest=end)
        except Exception as e:
            self.valid = False
            self.reason = f'Apero get failed with error: {e}'
            return

    def copy_tar(self):
        pass


    def _get_profile_params(self, params: dict) -> Dict[str, str]:
        """
        Get the profile parameters from the yaml file

        :param params: parameters dictionary

        :return: dictionary, the profile parameters
        """
        # deal with mode not in params
        if self.mode not in params:
            emsg = f'Profile {self.mode} not valid.'
            raise base.AperoRequestError(emsg)
        # ---------------------------------------------------------------------
        # load yaml file from params
        yaml_file = params[self.mode]
        # read yaml file
        yaml_params = io.read_yaml(yaml_file)
        # ---------------------------------------------------------------------
        # get profile name
        profile_name = self.mode.replace('_', ' ')
        # deal with profile name not in yaml_params
        if profile_name not in yaml_params:
            emsg = f'Profile name {profile_name} not in yaml file.'
            raise base.AperoRequestError(emsg)
        # set up storage for return
        profile_params = dict()
        # ---------------------------------------------------------------------
        # deal with general not in yaml_params[profile_name]
        if 'general' not in yaml_params[profile_name]:
            emsg = f'"general" not in yaml file for profile {profile_name}.'
            raise base.AperoRequestError(emsg)
        else:
            settings = yaml_params[profile_name]['general']
        # ---------------------------------------------------------------------
        # get apero profile
        if 'apero profile' not in settings:
            emsg = f'"apero profile" not in general for profile {profile_name}.'
            raise base.AperoRequestError(emsg)
        else:
            profile_params['apero profile'] = settings['apero profile']
        # ---------------------------------------------------------------------
        # get apero install
        if 'apero install' not in settings:
            emsg = f'"apero install" not in general for profile {profile_name}.'
            raise base.AperoRequestError(emsg)
        else:
            profile_params['apero install'] = settings['apero install']
        # ---------------------------------------------------------------------
        # return profile parameters
        return profile_params

    def email_success(self):
        pass

    def email_failure(self):
        pass


# =============================================================================
# Define functions
# =============================================================================
def get_sheetkind(params: Dict[str, Any], sheetkind: str) -> Tuple[str, str]:
    """
    Get the sheet id and sheet name for a given sheet kind

    :param params: dictionary, the parameter dictionary
    :param sheetkind: str, the sheet kind (must be "response" or "pass")

    :raises: AperoRequestError if sheetkind is not "response" or "pass"

    :return: tuple, 1. the sheet id, 2. the sheet name
    """
    # get parameters from params
    if sheetkind == 'response':
        sheet_id = params['response sheet id']
        sheet_name = params['response sheet name']
    elif sheetkind == 'pass':
        sheet_id = params['pass sheet id']
        sheet_name = params['pass sheet name']
    else:
        emsg = ('Sheet kind "{0}" not recognized must be "response" '
                'or "pass"'.format(sheetkind))
        raise base.AperoRequestError(emsg)
    # return sheet id and sheet name
    return sheet_id, sheet_name


def get_sheet(params: Dict[str, Any], sheetkind: str) -> pd.DataFrame:
    """
    Get the current google sheet for a given sheet kind

    :param params: dictionary, the parameter dictionary
    :param sheetkind: str, the sheet kind (must be "response" or "pass")

    :return: DataFrame, the current google sheet for "sheetkind"
    """
    # add gspread directory and auth files
    io.gsp_setup()
    # get sheet kind parameters
    sheet_id, sheet_name = get_sheetkind(params, sheetkind)
    # ---------------------------------------------------------------------
    # print progress
    msg = f'Downloading google-sheet: {sheetkind}'
    misc.log_msg(msg, level='info')
    # load google sheet instance
    google_sheet = gspd.spread.Spread(sheet_id)
    # convert google sheet to pandas dataframe
    current_dataframe = google_sheet.sheet_to_df(index=0, sheet=sheet_name)
    # return google sheet as current dataframe
    return current_dataframe


def update_response_sheet(params: Dict[str, Any], dataframe: pd.DataFrame):
    """
    Update the google sheet with the new dataframe

    :param params: dictionary, the parameter dictionary
    :param dataframe: dataframe, the new dataframe to push to google sheet

    :return: None, updates google sheet
    """
    # add gspread directory and auth files
    io.gsp_setup()
    # get sheet kind parameters
    sheet_id, sheet_name = get_sheetkind(params, 'response')
    # print progress
    msg = 'Pushing all rows to google-sheet'
    misc.log_msg(msg, level='info')
    # load google sheet instance
    google_sheet = gspd.spread.Spread(sheet_id)
    # push dataframe back to server
    google_sheet.df_to_sheet(dataframe, index=False, replace=True,
                             sheet=sheet_name)
    # print progress
    msg = 'All rows added to google-sheet'
    misc.log_msg(msg, level='info')


def create_requests(params: Dict[str, Any],
                    dataframe: pd.DataFrame
                    ) -> Tuple[List[Request], pd.DataFrame]:
    """
    Create a list of request instances from a dataframe

    final list should reject any rows older than a certain date

    also should reject any that have the wrong password

    :param params:
    :param dataframe:

    :return: the list of requests and the dataframe with the bad requests removed
    """

    # get pass key sheet
    pass_sheet = get_sheet(params, 'pass')
    # get the valid time period in days
    valid_time = pd.Timedelta(days=-params['valid time'])

    # convert timestampe to a datetime
    dataframe['Timestamp'] = pd.to_datetime(dataframe['Timestamp'],
                                            format="mixed", dayfirst=True)
    # get the offset
    offset = dataframe['Timestamp'] - pd.Timestamp.now() > valid_time
    # valid dataframe
    valid_dataframe = dataframe[offset]
    # mask out invalid passkeys
    mask = valid_dataframe['Pass key'].isin(pass_sheet['PASSKEY'])
    # cut down the valid dataframe
    valid_dataframe = valid_dataframe[mask]
    # -------------------------------------------------------------------------
    # we need to match the passkey + email address for each row in dataframe
    #   to the passkey + email address in the passkey sheet
    #   if we find a match we add it to the list of requests
    # -------------------------------------------------------------------------
    # turn pass_sheet into a dictionary of passkeys with the emails as
    #   lists for each key
    pass_dict = dict()
    runid_dict = dict()
    profile_dict = dict()
    for row in range(len(pass_sheet)):
        # get passkey
        passkey = str(pass_sheet['PASSKEY'][row])
        # get run id
        runid = np.char.array(pass_sheet['RUN_IDS'][row].split(','))
        runid.strip()
        # get emails
        emails = np.char.array(pass_sheet['EMAIL'][row].split(','))
        emails = emails.strip()
        # get profile_dict
        profiles = np.char.array(pass_sheet['APERO PROFILES'][row].split(','))
        profiles = profiles.strip()
        # push into pass_dict using a set
        if passkey not in pass_dict:
            pass_dict[passkey] = set(emails)
            runid_dict[passkey] = set(runid)
            profile_dict[passkey] = set(profiles)
        else:
            pass_dict[passkey] = pass_dict[passkey].union(emails)
            runid_dict[passkey] = runid_dict[passkey].union(runid)
            profile_dict[passkey] = profile_dict[passkey].union(profiles)
    # -------------------------------------------------------------------------
    # add two columns to the dataframe: ALLOWED_RUNID and ALLOWED_PROFILES
    # -------------------------------------------------------------------------
    valid_dataframe['ALLOWED_RUNID'] = [None] * len(valid_dataframe)
    valid_dataframe['ALLOWED_PROFILES'] = [None] * len(valid_dataframe)
    # store a list of requests
    requests, valid_mask = [], np.zeros(len(valid_dataframe)).astype(bool)
    # loop around rows in valid_dataframe
    for it, row in enumerate(valid_dataframe.index):
        # get passkey
        passkey = valid_dataframe['Pass key'][row]
        # get email
        email = valid_dataframe['Email Address'][row]
        # get profile requested
        profile = valid_dataframe['APERO_MODE'][row]

        # deal with passkey not in pass_dict
        if passkey not in pass_dict:
            # create request
            request = Request.from_pandas_row(valid_dataframe.iloc[it])
            # set request as invalid
            request.valid = False
            request.reason = f'Passkey {passkey} incorrect.'
            # add to list of requests
            requests.append(request)
            continue

        # check whether email is valid
        if email not in pass_dict[passkey]:
            # create request
            request = Request.from_pandas_row(valid_dataframe.iloc[it])
            # set request as invalid
            request.valid = False
            request.reason = f'Email {email} not in passkey {passkey}.'
            # add to list of requests
            requests.append(request)
            continue

        # check whether profile is valid
        if profile not in profile_dict[passkey]:
            # create request
            request = Request.from_pandas_row(valid_dataframe.iloc[it])
            # set request as invalid
            request.valid = False
            request.reason = f'Profile {profile} not in passkey {passkey}.'
            # add to list of requests
            requests.append(request)
            continue
        # ---------------------------------------------------------------------
        # valid requests
        # ---------------------------------------------------------------------
        # create request
        request = Request.from_pandas_row(valid_dataframe.iloc[it])
        # update the allowed run ids and profiles
        request.allowed_runids = runid_dict[passkey]
        request.allowed_profiles = profile_dict[passkey]
        # set request as valid
        request.valid = True
        # add to list of requests
        requests.append(request)
        # set valid mask
        valid_mask[it] = True
    # -------------------------------------------------------------------------
    return requests, valid_dataframe[valid_mask]


def update_apero_profile(params: dict) -> Any:
    """
    Update the apero profile with the correct paths

    :param profile: dict, the profile to update
    :return:
    """
    # use os to add DRS_UCONFIG to the path
    os.environ['DRS_UCONFIG'] = params['apero profile']
    # allow getting apero
    sys.path.append(params['apero install'])
    # ------------------------------------------------------------------
    # cannot import until after DRS_UCONFIG set
    from apero.base import base
    from apero.core import constants
    from apero.core.constants import param_functions
    from apero.core.utils import drs_startup
    # reload DPARAMS and IPARAMS
    base.DPARAMS = base.load_database_yaml()
    base.IPARAMS = base.load_install_yaml()
    # ------------------------------------------------------------------
    # invalidate cache
    param_functions.CONFIG_CACHE = dict()
    # make sure parameters is reloaded (and not cached)
    apero_params = constants.load(cache=False)
    # set apero pid
    apero_params['PID'], apero_params['DATE_NOW'] = drs_startup.assign_pid()
    # no inputs
    apero_params['INPUTS'] = dict()
    apero_params['OBS_DIR'] = None
    # return apero parameters
    return apero_params


def get_fibers(fibers: str) -> str:
    # import constants from apero
    from apero.core import constants
    # get pseudo constants
    pconst = constants.pload()
    # get fiber combinations
    science, reference = pconst.FIBER_CCF()
    all_fibers = pconst.INDIVIDUAL_FIBERS()

    # first case: All fibers
    if 'ALL' in fibers.upper():
        fiber_list = list(np.unique((all_fibers + [science, reference])))
    elif 'REFERENCE' in fibers.upper():
        fiber_list = [reference]
    elif 'SCIENCE' in fibers.upper():
        fiber_list = [science]
    else:
        emsg = f'Fiber column = {fibers} is invalid'
        return base.AperoRequestError(emsg)
    # if we have got to here return a str comma separated list
    return ','.join(fiber_list)


def get_runids(runids: set) -> Union[str, None]:
    if 'ALL' in runids:
        return None
    else:
        return ','.join(list(runids))


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # print hello world
    print('Hello World')

# =============================================================================
# End of code
# =============================================================================
