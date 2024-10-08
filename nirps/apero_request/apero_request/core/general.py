#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on {DATE}

@author: cook
"""
import time
import os
import sys
import shutil
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
# define rsync commands
RSYNC_CMD_IN = 'rsync -avh --delete -e "{SSH}" {USER}@{HOST}:{INPATH} {OUTPATH}'
RSYNC_CMD_OUT = 'rsync -avh --delete -e "{SSH}" {INPATH} {USER}@{HOST}:{OUTPATH}'
# define lock directory
LOCK_DIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
# define maximum lock out time in minutes
MAX_COUNT = 60
# contact emails
CONTACT = ['neil.cook@umontreal.ca', 'etienne.artigau@umontreal.ca',
           'charles.cadieux.1@umontreal.ca']


# =============================================================================
# Define classes
# =============================================================================
class Request:
    def __init__(self, timestamp: Union[str, pd.Timestamp],
                 email: str, drsobjn: str, dprtype: str, mode: str, fibers: str,
                 drsoutid: str, passkey: str,
                 start_date: Union[str, pd.Timestamp, None],
                 end_date: Union[str, pd.Timestamp, None]):
        """
        Create a request instance

        :param timestamp: str or pd.Timestamp, the timestamp of the request
        :param email: str, the email address of the user
        :param drsobjn: str, the drs object names (comma separated)
        :param dprtype: str, the DPRTYPES (comma separated)
        :param mode: str, the apero profile mode
        :param fiber: str, the fiber option (All fibers, Science fiber,
                      Reference fiber)
        :param drsoutid: str, the DRSOUTIDs (comma separated) file definitions
        :param passkey: str, the passkey
        :param start_date: str or pd.Timestamp, the start date or None
        :param end_date: str or pd.Timestamp, the end date or None
        """
        self.timestampe: Union[pd.Timestamp, None] = None
        self.email: Union[str, None] = None
        self.drsobjn: Union[np.chararray, None] = None
        self.dprtype: Union[np.chararray, None] = None
        self.mode: Union[str, None] = None
        self.fibers: str = ''
        self.drsoutid: Union[np.chararray, None] = None
        self.passkey: Union[str, None] = None
        self.start_date: Union[pd.Timestamp, None] = None
        self.end_date: Union[pd.Timestamp, None] = None
        self.timestamp_str: str = ''
        self.drsobjn_str: str = ''
        self.dprtype_str: str = ''
        self.drsoutid_str: str = ''
        self.fibers_str: str = ''
        self.start_date_str: str = ''
        self.end_date_str: str = ''
        # ---------------------------------------------------------------------
        # other attributes
        self.valid: bool = False
        self.skip = False
        self.exists: bool = False
        self.error: bool = False
        self.reason: str = ''
        self.allowed_runids: set = set()
        self.allowed_profiles: set = set()
        # to fill the new table
        self.df_row: Union[Any, None] = None
        # ---------------------------------------------------------------------
        self.set_params(timestamp, email, drsobjn, dprtype, mode, fibers,
                        drsoutid, passkey, start_date, end_date)
        # generate hash key for filename based on timestamp + email + passkey
        self.hashkey: str = self._generate_hashkey()
        # set filename from hashkey
        self.tarfile: str = f'apero-request-{self.hashkey}.tar.gz'
        # set url
        self.tar_url: str = ''
        # command that would have been run
        self.cmd = ''

    def set_params(self, timestamp: Union[str, pd.Timestamp],
                 email: str, drsobjn: str, dprtype: str, mode: str, fibers: str,
                 drsoutid: str, passkey: str,
                 start_date: Union[str, pd.Timestamp],
                 end_date: Union[str, pd.Timestamp]):
        """
        Try to fill the request parameters - done so any error will just lead
        to a reported error and not a crash

        :param timestamp: str or pd.Timestamp, the timestamp of the request
        :param email: str, the email address of the user
        :param drsobjn: str, the drs object names (comma separated)
        :param dprtype: str, the DPRTYPES (comma separated)
        :param mode: str, the apero profile mode
        :param fiber: str, the fiber option (All fibers, Science fiber,
                      Reference fiber)
        :param drsoutid: str, the DRSOUTIDs (comma separated) file definitions
        :param passkey: str, the passkey
        :param start_date: str or pd.Timestamp, the start date or None
        :param end_date: str or pd.Timestamp, the end date or None

        :return: None updates attrbiutes in self
        """
        try:
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
            # -----------------------------------------------------------------
            # formatted parameters as strings (for apero get)
            self.timestamp_str = Time(self.timestamp).iso
            self.drsobjn_str =  ','.join(list(self.drsobjn))
            self.dprtype_str = ','.join(list(self.dprtype))
            self.drsoutid_str = ','.join(list(self.drsoutid))
            if self.start_date is not None:
                self.start_date_str = Time(self.start_date).iso
            else:
                self.start_date_str = 'None'
            if self.end_date is not None:
                self.end_date_str = Time(self.end_date).iso
            else:
                self.end_date_str = 'None'
        except Exception as e:
            self.error = True
            self.reason = f'\tCould not set up apero get parameters. Error {e}'

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
        Run the request using apero_get

        :return:
        """
        # deal with already existing
        if self.exists:
            return
        # select profile parameters
        try:
            profile_params = self._get_profile_params(params)
        except Exception as e:
            self.valid = False
            self.reason = f'\tCould not get profile parameters with error: {e}'
            misc.log_msg(params, self.reason, log_only=True)
            return
        # update apero profile
        try:
            _ = update_apero_profile(profile_params)
        except Exception as e:
            self.valid = False
            self.reason = f'\tCould not update apero profile with error: {e}'
            misc.log_msg(params, self.reason, log_only=True)
            return
        # ---------------------------------------------------------------------
        # sort out arguments for apero get
        # ---------------------------------------------------------------------
        try:
            outpath = params['local path']
            self.fibers_str = get_fibers(self.fibers)
            runids = get_runids(self.allowed_runids)
            test = params.get('test', False)
        except Exception as e:
            self.valid = False
            self.reason = f'\tCould not set up apero get parameters. Error {e}'
            misc.log_msg(params, self.reason, log_only=True)
            return
        # ---------------------------------------------------------------------
        # try to run apero get
        # ---------------------------------------------------------------------
        try:
            # create command (for feedback)
            cmd = 'apero_get.py'
            cmd += f' --objnames={self.drsobjn_str}'
            cmd += f' --dprtypes={self.dprtype_str}'
            cmd += f' --outtypes={self.drsoutid_str}'
            cmd += f' --outpath={outpath}'
            cmd += f' --fibers={self.fibers_str}'
            cmd += f' --runid="{runids}"'
            if test:
                cmd += ' --test'
            if self.start_date_str not in [None, 'None', '']:
                cmd += f' --since={self.start_date_str}'
            if self.start_date_str not in [None, 'None', '']:
                cmd += f' --latest={self.end_date_str}'
            cmd += f' --timekey=observed'
            cmd += f' --sizelimit={params["file size limit"]}'
            self.cmd = cmd
            # need to import apero_get (for this profile)
            from apero.tools.recipes.bin import apero_get
            # print the command we are running
            message = f'Running command: {self.cmd}'
            misc.log_msg(params, message)
            # run apero get to make the objects dir in apero dir
            llget = apero_get.main(objnames=self.drsobjn_str,
                                   dprtypes=self.dprtype_str,
                                   outtypes=self.drsoutid_str,
                                   outpath=outpath,
                                   fibers=self.fibers_str,
                                   runid=runids,
                                   symlinks=False,
                                   tar=True, tarfile=self.tarfile,
                                   test=test,
                                   since=self.start_date_str,
                                   latest=self.end_date_str,
                                   timekey='observed',
                                   sizelimit=params['file size limit'])
        except Exception as e:
            self.valid = False
            self.reason = f'\tApero get failed with error: {e}'
            misc.log_msg(params, self.reason, log_only=True)
            return
        # ---------------------------------------------------------------------
        # set url (after creation)
        ssh_url = params['ssh']['url']
        if ssh_url.endswith('/'):
            self.tar_url = params['ssh']['url'] + self.tarfile
        else:
            self.tar_url = params['ssh']['url'] + '/' + self.tarfile
        # ---------------------------------------------------------------------
        # check that tar was create
        # ---------------------------------------------------------------------
        if not os.path.exists(self.tarfile):
            # get variables from llget
            get_success = llget['success']
            get_errors = llget['params']['LOGGER_ERROR']
            get_indict = llget['params']['INDICT']

            # if apero_get ran successfully then no files were found
            if get_success:
                if len(get_indict) == 0:
                    self.valid = False
                    self.reason = f'\tNo files found.'
                    misc.log_msg(params, self.reason, log_only=True)
                    return
            # else we report the error
            elif len(get_errors) > 0:
                self.valid = False
                self.reason = '\n'.join(get_errors)
                misc.log_msg(params, self.reason, log_only=True)
                return
            else:
                # if we get here we have an unknown error
                self.valid = False
                self.reason = '\tUnknown error. Contact support with the query.'
                misc.log_msg(params, self.reason, log_only=True)
        else:
            return

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
        yaml_file = params[self.mode]['yaml']
        # read yaml file
        yaml_params = io.read_yaml(yaml_file)
        # ---------------------------------------------------------------------
        # get profile name
        profile_name = params[self.mode]['profile']
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

    def _generate_hashkey(self) -> str:
        """
        Generate a hash key for the request using all the requested parameters
        Should always be unique unless the same request was dupilcated twice

        :return: str, the unique hash key
        """
        # create a unique hash
        hash = misc.generate_arg_checksum(source=[self.timestamp_str,
                                                  self.email,
                                                  self.passkey,
                                                  self.drsobjn_str,
                                                  self.dprtype_str,
                                                  self.mode,
                                                  self.fibers,
                                                  self.drsoutid_str,
                                                  self.start_date_str,
                                                  self.end_date_str],
                                          ndigits=20)
        return hash

    def _generate_summary(self):
        summary = ''
        summary += f'tarfile: {self.tarfile}\n\n'
        summary += '='*50 + '\n'
        summary += 'Request\n'
        summary += '='*50 + '\n'
        summary += f'Timestamp: {self.timestamp_str}\n'
        summary += f'Email: {self.email}\n'
        summary += f'Passkey: {self.passkey}\n'
        summary += f'DRSOBJN: {self.drsobjn_str}\n'
        summary += f'DPRTYPE: {self.dprtype_str}\n'
        summary += f'APERO_MODE: {self.mode}\n'
        if len(self.fibers_str) == 0:
            fiber_str = self.fibers
        else:
            fiber_str = f'{self.fibers} ({self.fibers_str})'
        summary += f'FIBERS: {fiber_str}\n'
        summary += f'DRSOUTID: {self.drsoutid_str}\n'
        summary += f'START_DATE: {self.start_date_str}\n'
        summary += f'END_DATE: {self.end_date_str}\n'

        if len(self.cmd) > 0:
            summary += f''
            summary += f'APERO GET CMD:'
            summary += f''
            summary += f'\t{self.cmd}'
        return summary

    def __str__(self) -> str:
        return self._generate_summary()

    def __repr__(self) -> str:
        return self._generate_summary()

    def make_row(self) -> pd.Series:
        row = pd.Series()
        row['Timestamp'] = self.timestamp_str
        row['Email Address'] = self.email
        row['DRSOBJN'] = self.drsobjn_str
        row['DPRTYPE'] = self.dprtype_str
        row['APERO_MODE'] = self.mode
        row['FIBER'] = self.fibers
        row['DRSOUTID'] = self.drsoutid_str
        row['Pass key'] = self.passkey
        row['START_DATE'] = self.start_date_str
        row['END_DATE'] = self.end_date_str
        return row

    def email_success(self, params: Dict[str, Any], iteration: int):

        valid_time = params['valid time']

        email_string = ''
        email_string += f'Your APERO request was successful.\n\n'
        email_string += f'You can download your tar file from:\n'
        email_string += f'{self.tar_url}\n\n'
        email_string += f'Note this URL will expire in {valid_time} days\n'
        email_string += f'Your request was:\n'
        email_string += self._generate_summary()

        # if self.valid:
        #     print('Success!')
        #     print(self._generate_summary())

        email_kwargs = dict()
        email_kwargs['message'] = email_string
        email_kwargs['people'] = self.email
        email_kwargs['sender_address'] = params['email from']
        email_kwargs['subject'] = f'APERO request successful: {self.hashkey}'
        try:
            misc.send_email(**email_kwargs)
            msg = 'Succeeded to send email for request {0}'
            msg += self._generate_summary()
            margs = [iteration]
            misc.log_msg(params, msg.format(*margs), log_only=True)
        except Exception as e:
            msg = 'Failed to send email for request {0}'
            msg += '\n\t{1}: {2}'
            msg += '\n\t Request was:\n'
            msg += self._generate_summary()
            margs = [iteration, type(e), str(e)]
            misc.log_msg(params, msg.format(*margs))

    def email_failure(self, params: Dict[str, Any], iteration: int):

        email_string = ''
        email_string += f'Your APERO request failed.\n\n'
        email_string += f'The reason for the failure was:\n'
        email_string += f'{self.reason}\n\n'
        email_string += f'Your request was:\n'
        email_string += self._generate_summary()
        email_string += ('\n\nIf you believe your query was valid '
                         '(or file size was too large) please contact the '
                         'admins: {0}'.format(', '.join(CONTACT)) + '\n\n')

        # if not self.valid:
        #     print('Failure!')
        #     print(self.reason)
        #     print(self._generate_summary())

        email_kwargs = dict()
        email_kwargs['message'] = email_string
        email_kwargs['people'] = self.email
        email_kwargs['sender_address'] = params['email from']
        email_kwargs['subject'] = f'APERO request failed: {self.hashkey}'

        try:
            misc.send_email(**email_kwargs)
        except Exception as e:
            msg = 'Failed to send email for request {0}'
            msg += '\n\t{1}: {2}'
            msg += '\n\t Request was:\n'
            msg += self._generate_summary()
            margs = [iteration, type(e), str(e)]
            misc.log_msg(params, msg.format(*margs))


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
    elif sheetkind == 'archive':
        sheet_id = params['archive sheet id']
        sheet_name = params['archive sheet name']
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
    msg = f'\tDownloading google-sheet: {sheetkind}'
    misc.log_msg(params, msg, level='info')
    # load google sheet instance
    google_sheet = gspd.spread.Spread(sheet_id)
    # convert google sheet to pandas dataframe
    current_dataframe = google_sheet.sheet_to_df(index=0, sheet=sheet_name)
    # return google sheet as current dataframe
    return current_dataframe


def copy_index(params: Dict[str, Any]):
    # get path
    path = os.path.dirname(os.path.dirname(__file__))
    # get index file
    index_file = os.path.join(path, 'resources', 'index.html')
    # copy to the local path
    if not params['test']:
        shutil.copy(index_file, params['local path'])

def sync_data(params: Dict[str, Any]):
    rdict = dict()
    rdict['SSH'] = params['ssh']['options']
    rdict['USER'] = params['ssh']['user']
    rdict['HOST'] = params['ssh']['host']
    rdict['INPATH'] = params['local path']
    rdict['OUTPATH'] = os.path.join(params['ssh']['directory'])
    # print command to rsync
    misc.log_msg(params, RSYNC_CMD_OUT.format(**rdict))
    # run command (will require password)
    if not params['test']:
        os.system(RSYNC_CMD_OUT.format(**rdict))


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
    # get the current data from the google sheet (in case it has been updated)
    current_dataframe = get_sheet(params, 'response')
    # convert timestampe to a datetime
    current_dataframe['Timestamp'] = pd.to_datetime(current_dataframe['Timestamp'],
                                                    format="mixed", dayfirst=True)
    # -------------------------------------------------------------------------
    # get all entries after MOST RECENT
    new_entries = current_dataframe['Timestamp'] > params['MOST_RECENT']
    # if we have new entries add them
    if np.sum(new_entries) > 0:
        # print progress
        msg = ('New entries found since running request script.'
               '\n\tAdding {0} new entries').format(np.sum(new_entries))
        misc.log_msg(params, msg, level='info')
        # add the new entries
        dataframe = pd.concat([dataframe, current_dataframe[new_entries]],
                              ignore_index=True)
    # -------------------------------------------------------------------------
    # print progress
    msg = 'Pushing all rows to google-sheet'
    misc.log_msg(params, msg, level='info')
    # load google sheet instance
    google_sheet = gspd.spread.Spread(sheet_id, sheet=sheet_name)
    # get worksheet (we need to manually delete rows for a google form
    #   this is due to gspread_pandas trying to delete the header row
    worksheet = google_sheet.spread.worksheet('RESPONSE')
    # cannot delete row 1 as it is the header
    try:
        # hack - make sure we have at least one row
        worksheet.add_rows(1)
        # remove all rows
        worksheet.delete_rows(2, worksheet.row_count+1)
        # add a blank row
        worksheet.add_rows(1)
    except Exception as _:
        # this only happens when there are no rows
        pass
    # push dataframe back to server
    if len(dataframe) > 0:
        google_sheet.df_to_sheet(dataframe, index=False, replace=False)
    # print progress
    msg = 'Valid rows added to response google-sheet'
    misc.log_msg(params, msg, level='info')


def update_archive_sheet(params: Dict[str, Any], dataframe: pd.DataFrame):
    """
    Update the google sheet with the new dataframe

    :param params: dictionary, the parameter dictionary
    :param dataframe: dataframe, the new dataframe to push to google sheet

    :return: None, updates google sheet
    """
    # add gspread directory and auth files
    io.gsp_setup()
    # get sheet kind parameters
    sheet_id, sheet_name = get_sheetkind(params, 'archive')
    # print progress
    msg = 'Pushing all rows to google-sheet'
    misc.log_msg(params, msg, level='info')
    # make a hash key for each row in dataframe
    hashkeys = []
    for row in range(len(dataframe)):
        source = ''
        for col in dataframe.columns:
            # skip these columns (can change)
            if col in ['VALID', 'REASON', 'HASHKEY']:
                continue
            source += str(dataframe[col][row])
        # get hashkey
        hashkey = misc.generate_arg_checksum(source)
        # append to list
        hashkeys.append(hashkey)
    # push hashkeys into dataframe
    dataframe['HASHKEY'] = hashkeys

    # load google sheet instance
    google_sheet = gspd.spread.Spread(sheet_id, sheet=sheet_name)

    # get current data frame
    current_dataframe = google_sheet.sheet_to_df(index=0, sheet=sheet_name)

    # remove all rows from dataframe that are already in current_dataframe
    mask = np.in1d(dataframe['HASHKEY'], current_dataframe['HASHKEY'])
    dataframe = dataframe[~mask]

    # push dataframe back to server
    if len(dataframe) > 0:
        new_dataframe = pd.concat([current_dataframe, dataframe],
                                  ignore_index=True)

        google_sheet.df_to_sheet(new_dataframe, index=False, replace=True)
    # print progress
    msg = 'All rows added to archive google-sheet'
    misc.log_msg(params, msg, level='info')


def create_requests(params: Dict[str, Any],
                    dataframe: pd.DataFrame) -> List[Request]:
    """
    Create a list of request instances from a dataframe

    final list should reject any rows older than a certain date

    also should reject any that have the wrong password

    :param params:
    :param dataframe:

    :return: the list of requests (valid and invalid)
    """

    # get pass key sheet
    pass_sheet = get_sheet(params, 'pass')
    # get the valid time period in days
    valid_time = pd.Timedelta(days=-params['valid time'])
    # get the offset
    offset = dataframe['Timestamp'] - pd.Timestamp.now() > valid_time
    # valid dataframe
    valid_dataframe = dataframe[offset]
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
    # store a list of requests
    requests = []
    # loop around rows in valid_dataframe
    for it, row in enumerate(valid_dataframe.index):
        # get passkey
        passkey = valid_dataframe['Pass key'][row]
        # get email
        email = valid_dataframe['Email Address'][row]
        # get profile requested
        profile = valid_dataframe['APERO_MODE'][row]
        # ---------------------------------------------------------------------
        # skip profiles not in filter profiles list (if used)
        #    if this is not used then all profiles are valid
        if params['filter profiles'] is not None:
            if profile not in params['filter profiles']:
                # create request
                request = Request.from_pandas_row(valid_dataframe.iloc[it])
                # set request as invalid
                request.skip = True
                # add to list of requests
                requests.append(request)
                continue
        # ---------------------------------------------------------------------
        # deal with passkey not in pass_dict
        if passkey not in pass_dict:
            # create request
            request = Request.from_pandas_row(valid_dataframe.iloc[it])
            # set request as invalid
            request.valid = False
            request.reason = f'\tPasskey "{passkey}" incorrect.'
            # add to list of requests
            requests.append(request)
            continue

        # check whether email is valid
        if email not in pass_dict[passkey]:
            # create request
            request = Request.from_pandas_row(valid_dataframe.iloc[it])
            # set request as invalid
            request.valid = False
            request.reason = f'\tEmail {email} not valid for passkey "{passkey}".'
            # add to list of requests
            requests.append(request)
            continue

        # check whether profile is valid
        if profile not in profile_dict[passkey]:
            # create request
            request = Request.from_pandas_row(valid_dataframe.iloc[it])
            # set request as invalid
            request.valid = False
            request.reason = f'\tProfile {profile} not valid for passkey "{passkey}".'
            # add to list of requests
            requests.append(request)
            continue
        # ---------------------------------------------------------------------
        # valid requests
        # ---------------------------------------------------------------------
        # create request
        request = Request.from_pandas_row(valid_dataframe.iloc[it])
        # deal with errors
        if request.error:
            # request is not valid
            request.valid = False
            # add to list of requests
            requests.append(request)
        else:
            # update the allowed run ids and profiles
            request.allowed_runids = runid_dict[passkey]
            request.allowed_profiles = profile_dict[passkey]
            # add the dataframe row info
            request.df_row = valid_dataframe.iloc[it]
            # set request as valid
            request.valid = True
            # add to list of requests
            requests.append(request)
    # -------------------------------------------------------------------------
    return requests


def remove_invalid_tars(params: Dict[str, Any],
                        requests: List[Request]) -> List[Request]:
    # construct tar directory
    tar_dir = params['local path']
    # get a list of tar files
    tar_files_on_disk = os.listdir(tar_dir)
    # storage for valid requests
    valid_requests = dict()
    # get list of valid request files
    for it, request in enumerate(requests):
        if request.valid or request.skip:
            valid_requests[it] = request.tarfile
    # get a list of all tar files
    valid_request_tar_files = list(valid_requests.values())
    # loop around tar files on disk
    for tar_file in tar_files_on_disk:
        # case 1: tar file valid --> we do not reprocess
        if tar_file in valid_request_tar_files:
            # find the request that matches this tar file
            for it in list(valid_requests.keys()):
                if tar_file == requests[it].tarfile:
                    # set that this is valid
                    requests[it].valid = True
                    # set exists to True
                    requests[it].exists = True
                    break
        # case2: tar file is not valid --> remove
        elif tar_file.endswith('.tar.gz'):
            # remove file
            os.remove(os.path.join(tar_dir, tar_file))
    # return updated requests
    return requests


def create_dataframe(params: Dict[str, Any], requests: List[Request]
                     ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # loop around the requests and create a list of rows
    valid_rows, all_rows = [], []
    for request in requests:
        # get the row from the request
        row = request.df_row
        # deal with not having a row (i.e. it is None)
        if row is None:
            row = request.make_row()
        # copy row (for all row dataframe)
        all_row = row.copy()
        # add validity and reason for failure (if failed)
        all_row['VALID'] = request.valid
        all_row['REASON'] = request.reason
        # add to all rows
        all_rows.append(all_row)
        # only if request if valid
        if not request.valid and not request.skip:
            continue
        # add to rows
        valid_rows.append(row)
    # convert list of rows to dataframe
    valid_dataframe = pd.DataFrame(valid_rows)
    all_dataframe = pd.DataFrame(all_rows).reset_index(drop=True)
    # return dataframe
    return valid_dataframe, all_dataframe


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
        raise base.AperoRequestError(emsg)
    # if we have got to here return a str comma separated list
    return ','.join(fiber_list)


def get_runids(runids: set) -> Union[str, None]:
    if 'ALL' in runids:
        return None
    else:
        return ','.join(list(runids))




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
    # wait a very small random amount of time (to avoid two codes doing this
    # at the same time)
    random_time = 0.1 * np.random.random()
    time.sleep(random_time)
    # if lock file exists
    while os.path.exists(lock_file):
        # if stop we assume we return straight away with a False
        if stop:
            return False
        # deal with a max count
        if counter > (MAX_COUNT * 60):
            emsg = f'Lock file timeout. Please remove: {lock_file} manually'
            raise base.AperoRequestError(emsg)
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
if __name__ == "__main__":
    # print hello world
    print('Hello World')

# =============================================================================
# End of code
# =============================================================================
