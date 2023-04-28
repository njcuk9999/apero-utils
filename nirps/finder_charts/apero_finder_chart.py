from apero.tools.module.database import manage_databases
# must import here (so that os.environ is set)
# noinspection PyPep8Naming
from apero import lang
from apero.core import constants
from apero.core.utils import drs_startup



# =============================================================================
# Define variables
# =============================================================================


# =============================================================================
# Define functions
# =============================================================================



# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':
    # get the parameter dictionary of constants from apero
    params = constants.load()
    # Get the text types
    textentry = lang.textentry
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # get the object database (combined with pending + user table)
    maintable = manage_databases.get_object_database(params, log=True)

    print(maintable)

