import argparse
from typing import Any

from astropy.table import Table

from apero import lang
from apero.core import constants
from apero.core.utils import drs_startup
from apero.core.core import drs_log
from apero.tools.module.database import manage_databases


# =============================================================================
# Define variables
# =============================================================================
# get WLOG
WLOG = drs_log.wlog

# =============================================================================
# Define functions
# =============================================================================
def get_args():
    """
    Define the command line arguments

    :return: argparse namespace object containing arguments
    """
    parser = argparse.ArgumentParser(description='Rene\'s magic trigger')
    # add obs dir
    parser.add_argument('objname', type=str, default='*',
                        help='The object name')
    parser.add_argument('--date', type=str, default='*',
                        help='The date of observation')
    # load arguments with parser
    args = parser.parse_args()
    # return arguments
    return args


def find_objname(pconst: Any, objname: str, full_table: Table) -> str:
    """
    Find the object name in the object table (from OBJNAME or ORIGINAL_NAME
    or ALIASES column)

    :param objname:
    :param full_table:
    :return:
    """
    # clean the input objname
    cobjname = pconst.DRS_OBJ_NAME(objname)
    # -------------------------------------------------------------------------
    # if objname in OBJNAME column return the objname - no need to find it
    if cobjname in full_table['OBJNAME']:
        return cobjname
    # -------------------------------------------------------------------------
    # deal with a null object (should not continue)
    if cobjname == 'Null':
        WLOG(params, 'error', 'Object not found')
        return ''
    # -------------------------------------------------------------------------
    # assume we have not found our object name
    found = False
    # get aliases
    aliases = full_table['ALIASES']
    # set row to zero as a place holder
    row = 0
    # loop around each row in the table
    for row in range(len(aliases)):
        # loop around aliases until we find the alias
        for alias in aliases[row].split('|'):
            if pconst.DRS_OBJ_NAME(alias) == cobjname:
                found = True
                break
        # stop looping if we have found our object
        if found:
            break
    # get the cobjname for this target if found
    if found:
        return full_table['OBJNAME'][row]
    else:
        WLOG(params, 'error', 'Cannot find object')
        return ''


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':
    # get the parameter dictionary of constants from apero
    params = constants.load()
    pconst = constants.pload()
    # Get the text types
    textentry = lang.textentry
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # get the object database (combined with pending + user table)
    maintable = manage_databases.get_object_database(params, log=True)

    print(maintable)
