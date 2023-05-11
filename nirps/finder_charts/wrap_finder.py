from typing import Any, Dict

import pandas as pd
from astropy.time import Time

import apero_finder_chart as afc

# =============================================================================
# Define variables
# =============================================================================
# set the directory finder chart pdfs will be saved to
DIRECTORY = '/cosmos99/nirps/apero-data/misc/apero-find/'
# set the date for finder charts
DATE = Time.now()


# =============================================================================
# Define functions
# =============================================================================
def from_apero(it: int, object_table: pd.DataFrame) -> Dict[str, Any]:
    """
    Create a objdict for use in apero_finder_chart.main() from this row of the
    full object table

    :param it: int, the row of the table to create dictionary for
    :param object_table: astropy.Table, the full object database table

    :return: Dict, the objdict ofr use in apero_finder_chart.main()
    """
    # get the row from the object_table
    objdata = object_table.iloc[it]
    # set up the storage
    objdict = dict()
    # get the object name
    objdict['OBJNAME'] = objdata['OBJNAME']
    # get the object ra and dec
    objdict['RA_DEG'] = objdata['RA_DEG']
    objdict['DEC_DEG'] = objdata['DEC_DEG']
    # get the object epoch (jd)
    objdict['EPOCH'] = objdata['EPOCH']
    # get the object parallax
    objdict['PLX'] = objdata['PLX']
    # get the object proper motion
    objdict['PMRA'] = objdata['PMRA']
    objdict['PMDE'] = objdata['PMDE']
    # return the object dictionary
    return objdict


def main():
    """
    Get all objects and push into finder chart one by one
    :return:
    """
    try:
        from apero.core import constants
        from apero.core.core import drs_database
        from apero.core.core import drs_log
        from apero.core.utils import drs_startup
    except ImportError:
        emsg = 'Cannot use APERO. Please install it to use --apero'
        raise ImportError(emsg)
    # get the parameter dictionary of constants from apero
    params = constants.load()
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # print progress
    print('Getting objects from database')
    # load object database
    objdbm = drs_database.AstrometricDatabase(params)
    objdbm.load_db()
    # get all objects
    object_table = objdbm.get_entries('*')
    # loop around objnames
    for it, objname in enumerate(object_table['OBJNAME']):
        # get the objdict
        objdict = from_apero(it, object_table)
        # print progress
        print('\n' + '=' * 50)
        args = [objname, it + 1, len(object_table)]
        print('Running NIRPS finder for object: {0} [{1}/{2}]'.format(*args))
        print('=' * 50 + '\n')
        # run the main code
        afc.main(objname, DATE, objdict=objdict, directory=DIRECTORY)


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':
    # run main code
    main()

# =============================================================================
# End of code
# =============================================================================
