import os
import shutil
from typing import Any, Dict

import yaml
from apero_ri.tasks import apero_sync
from astropy.time import Time

from manual_trigger import update_apero_profile

# =============================================================================
# Define variables
# =============================================================================
# Path the the apero_ri (ARI2) profile files (relative to the apero install path)
ARI_PROFILE_REL_PATH = ['apero-ri', 'apero_ri', 'resources',
                        'aprofile_instruments']

# list of local tasks to run
LOCAL_TASKS = ['APERO_OBJECT_QUERY', 'APERO_QC_STATS']


# =============================================================================
# Define functions
# =============================================================================
def v08_settings(params):

    from aperocore import base as ac_base
    # load DPARAMS and IPARAMS
    dparams = ac_base.load_database_yaml()


    database_dict = dict()
    database_dict['DATABASE_MODE'] = dparams['TYPE']
    database_dict['DATABASE_HOST'] = dparams['HOST']
    database_dict['DATABASE_USER'] = dparams['USER']
    database_dict['DATABASE_PASSWORD'] = dparams['PASSWD']
    database_dict['DATABASE_NAME'] = dparams['DATABASE']
    database_dict['FINDEX_TABLENAME'] = dparams['FINDEX']['TABLE']
    database_dict['ASTROM_TABLENAME'] = dparams['ASTROM']['TABLE']
    database_dict['CALIB_TABLENAME'] = dparams['CALIB']['TABLE']
    database_dict['LOG_TABLENAME'] = dparams['LOG']['TABLE']
    database_dict['TELLU_TABLENAME'] = dparams['TELLU']['TABLE']
    database_dict['REJECT_TABLENAME'] = dparams['REJECT']['TABLE']

    path_dict = dict()
    path_dict['PATH_RAW'] = params['PATH.RAW']
    path_dict['PATH_PP'] = params['PATH.PP']
    path_dict['PATH_RED'] = params['PATH.RED']
    path_dict['PATH_LOG'] = params['PATH.LOG']
    path_dict['PATH_OUT'] = params['PATH.OUT']
    path_dict['PATH_LBL'] = params['PATH.LBL']
    path_dict['PATH_CALIB'] = params['PATH.CALIB']
    path_dict['PATH_TELLU'] = params['PATH.TELLU']

    return database_dict, path_dict


def v07_settings(params):

    from apero.base import base
    # reload DPARAMS and IPARAMS
    dparams = base.load_database_yaml()

    if dparams['USE_MYSQL']:
        dparams = dparams['MYSQL']
        db_mode = 'mysql+pymysql'
    else:
        dparams = dparams['SQLITE3']
        db_mode = 'sqlite'

    # get table names
    tparams = get_db_tablenames(dparams)

    database_dict = dict()
    database_dict['DATABASE_MODE'] = db_mode
    database_dict['DATABASE_HOST'] = dparams['HOST']
    database_dict['DATABASE_USER'] = dparams['USER']
    database_dict['DATABASE_PASSWORD'] = dparams['PASSWD']
    database_dict['DATABASE_NAME'] = dparams['DATABASE']
    database_dict['FINDEX_TABLENAME'] = tparams['findex']
    database_dict['ASTROM_TABLENAME'] = tparams['astrom']
    database_dict['CALIB_TABLENAME'] = tparams['calib']
    database_dict['LOG_TABLENAME'] = tparams['log']
    database_dict['TELLU_TABLENAME'] = tparams['tellu']
    database_dict['REJECT_TABLENAME'] = tparams['reject']

    path_dict = dict()
    path_dict['PATH_RAW'] = params['DRS_DATA_RAW']
    path_dict['PATH_PP'] = params['DRS_DATA_WORKING']
    path_dict['PATH_RED'] = params['DRS_DATA_REDUC']
    path_dict['PATH_LOD'] = params['DRS_DATA_MSG']
    path_dict['PATH_OUT'] = params['DRS_DATA_OUT']
    path_dict['PATH_LBL'] = params['LBL_PATH']
    path_dict['PATH_CALIB'] = params['DRS_CALIB_DB']
    path_dict['PATH_TELLU'] = params['DRS_TELLU_DB']

    return database_dict, path_dict


def get_db_tablenames(dparams) -> Dict[str, str]:
    from apero.base import base
    # storag for return
    tablenames = dict()
    # loop around database names
    for dbname in base.DATABASE_NAMES:
        # get yaml key
        ydbname = dbname.upper()
        # construct table name
        tablename = '{0}_{1}_db'.format(dbname, dparams[ydbname]['PROFILE'])
        # push into storage
        tablenames[dbname] = tablename
    # return table names dict
    return tablenames


def load_apero_ri_resource_profiles(pdict: Dict[str, Any],
                                    aprofiles: Dict[str, Any]):
    # get apero install path
    apero_install_path = pdict['general']['apero install']
    # get ari profile name
    aprofile_name = str(pdict['ari']['ari profile'])
    # construct absolute path to the yaml
    yamldir = str(os.path.join(apero_install_path, *ARI_PROFILE_REL_PATH))
    # add the file name to the yaml directory
    abspath = os.path.join(yamldir, aprofile_name)
    # deal with no file found
    if not os.path.exists(abspath):
        raise FileNotFoundError(f'ARI profile yaml not found at {abspath}')
    # read the contents of the yaml
    with open(abspath, 'r', encoding='utf-8') as yamlfile:
        raw_yaml = yaml.safe_load(yamlfile)
    # push into aprofiles
    for key in raw_yaml:
        aprofiles[key] = raw_yaml[key]
    # push these up a level
    # TODO: THESE ARE BAD HACKS from bad ARI coding
    aprofiles['SCIENCE_TYPES'] = aprofiles['general']['science_types']
    aprofiles['SCIENCE_FIBER'] = aprofiles['general']['science_fiber']

    aprofiles['general']['SCIENCE_TYPES'] = aprofiles['general']['science_types']
    aprofiles['general']['SCIENCE_FIBER'] = aprofiles['general']['science_fiber']
    # return profiles
    return aprofiles


def main(settings):

    # loop around profiles
    for profile in settings['PROFILES']:
        # get the yaml dictionary for this profile
        pdict = settings['PROFILES'][profile]
        # update the apero profile
        aparams = update_apero_profile(pdict)

        # construct the apero sync dictionary
        rparams = dict()

        # ---------------------------------------------------------------------
        # Set up local ARI directory
        # ---------------------------------------------------------------------
        local_ari_dir = pdict['ari'].get('ari path', None)
        # deal with no local ARI directory ste
        if local_ari_dir is None:
            local_ari_dir = os.path.expanduser('~/.ari')
        # deal with local ARI directory not existing
        if not os.path.exists(local_ari_dir):
            os.makedirs(local_ari_dir)
        # make log file
        logfile = f'{Time.now().fits}_apero_sync.log'
        logpath = os.path.join(local_ari_dir, 'logs', logfile)
        # set the local ARI directory in the apero sync dictionary
        rparams['LOCAL_DATA_DIR'] = os.path.join(local_ari_dir, 'sync')
        # clean out directory
        if os.path.exists(rparams['LOCAL_DATA_DIR']):
            shutil.rmtree(rparams['LOCAL_DATA_DIR'])
        # remake the sync directory
        os.makedirs(rparams['LOCAL_DATA_DIR'])

        # ---------------------------------------------------------------------
        # set up database/general/path settings
        # ---------------------------------------------------------------------
        if pdict['general']['apero version'].startswith('0.7'):
            database_dict,path_dict = v07_settings(aparams)
        else:
            database_dict, path_dict = v08_settings(aparams)

        aprofiles = dict(database=database_dict,
                         paths=path_dict)
        # Load aprofiles (from apero-ri)
        aprofiles = load_apero_ri_resource_profiles(pdict, aprofiles)
        # set instrument at global level
        rparams['INSTRUMENT'] = aprofiles['general']['instrument']
        # ---------------------------------------------------------------------
        # push into rparams
        rparams['APERO_PROFILES'] = dict()
        rparams['APERO_PROFILES'][profile] = aprofiles
        # ---------------------------------------------------------------------
        # setup task config
        rparams['TASK_CONFIG'] = dict()
        rparams['TASK_CONFIG']['ncores'] = pdict['ari']['cores']
        rparams['TASK_CONFIG']['mp_backend'] = pdict['ari']['mp_backend']
        rparams['TASK_CONFIG']['mp_start_method'] = pdict['ari']['mp_start_method']

        # ---------------------------------------------------------------------
        # run the sync code(s)
        for local_task in LOCAL_TASKS:
            apero_sync.run(local_task, rparams, verbose=True,
                           log_file=logpath)


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # fake trigger settings
    trigger_settings = dict()
    trigger_settings['PROFILES'] = dict()

    _pp = dict()
    _pp['general'] = dict()
    _pp['general']['apero profile'] = '/scratch2/spirou/drs-settings/spirou_mini2_07'
    _pp['general']['apero install'] = '/scratch2/spirou/drs-bin/apero-drs-spirou-07XXX/'
    _pp['general']['apero version'] = '0.7.296'

    _pp['ari'] = dict()
    _pp['ari']['ari path'] = '/scratch2/spirou/drs-data/spirou_mini2_07/other/ari-local'
    _pp['ari']['ari profile'] = 'spirou_v7.yaml'
    _pp['ari']['cores'] = 5
    _pp['ari']['mp_backend'] = 'process'
    _pp['ari']['mp_start_method'] = 'fork'

    trigger_settings['PROFILES']['spirou_mini2_07'] = _pp
    # ----------------------------------------------------------------------
    main(trigger_settings)

# =============================================================================
# End of code
# =============================================================================
