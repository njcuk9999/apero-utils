#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on {DATE}

@author: cook
"""
from typing import Any, Dict, List, Tuple
import pandas as pd
import gspread_pandas as gspd

from apero_request.core import base
from apero_request.core import misc
from apero_request.core import io


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
    def __init__(self, timestamp, email, drsobjn, dprtype, online, offline,
                 drsoutid, passkey, start_date, end_date):
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
                    dataframe: pd.DataFrame) -> List[Request]:
    """
    Create a list of request instances from a dataframe

    final list should reject any rows older than a certain date

    final list should remove any rows that are currently on disk

    :param params:
    :param dataframe:
    :return:
    """

    # get pass key sheet
    pass_sheet = get_sheet(params, 'pass')
    # get the valid time period in days
    valid_time = params['valid time']

    for row in dataframe:
        pass


    return []

# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # print hello world
    print('Hello World')

# =============================================================================
# End of code
# =============================================================================
