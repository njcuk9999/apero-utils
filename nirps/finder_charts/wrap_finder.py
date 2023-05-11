import os

from astropy.time import Time

import apero_finder_chart as afc

# =============================================================================
# Define variables
# =============================================================================
# set the directory finder chart pdfs will be saved to
DIRECTORY = '/cosmos99/nirps/apero-data/misc/apero-find/'
# set the date for finder charts
DATE = Time.now()
# Skip done
SKIP_DONE = True

# =============================================================================
# Define functions
# =============================================================================
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
        # args for printout
        args = [objname, it + 1, len(object_table)]
        # check whether file already exists
        abspath = afc.construct_savepath(DIRECTORY, DATE, objname)
        if SKIP_DONE and os.path.exists(abspath):
            msg = 'Skipping NIRPS finder for object: {0} [{1}/{2}]'
            print(msg.format(*args))
            continue
        # get the objdict
        objdict = afc.from_apero_objtable(it, object_table)
        # print progress
        print('\n' + '=' * 50)
        msg = 'Running NIRPS finder for object: {0} [{1}/{2}]'
        print(msg.format(*args))
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
