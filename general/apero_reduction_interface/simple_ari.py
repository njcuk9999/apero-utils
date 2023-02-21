#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-02-20 at 9:35

@author: cook
"""
import glob
import os
import re
from typing import Union

import numpy as np
import yaml
from astropy.table import Table

# =============================================================================
# Define variables
# =============================================================================
# define profile yaml file
PROFILE_FILE = 'profiles.yaml'
# define the astrometric database column names to get
ASTROMETRIC_COLUMNS = ['OBJNAME', 'RA_DEG', 'DEC_DEG', 'TEFF', 'SP_TYPE']
ASTROMETRIC_DTYPES = [str, float, float, float, str]
# define the log database column names to get
LOG_COLUMNS = ['SHORTNAME', 'RUNSTRING', 'START_TIME', 'END_TIME',
               'STARTED', 'PASSED_ALL_QC', 'ENDED']
# define which instruments have polar
HAS_POLAR = dict(SPIROU=True, NIRPS_HE=False, NIRPS_HA=False)
# define log levels to report
LOG_LEVELS = ['error', 'warning']
# define whether to reprocess stats
REPROCESS = True
# Define output path
OUTPUT_PATH = '.'
# define sphinx directories
_WORKING_DIR = '/data/spirou/cook/simple_ari/working/'
_DATA_DIR = 'data'
_RST_DIR = 'rst'
_HTML_DIR = '_build/html'
_OUT_DIR = 'output'
# ssh path
SSH_PATH = '/export/www/home/cook/www/apero-drs/'
# define path to local htpasswd file
HTPASSWD_PATH = '/home/cook/apero_ari/'
_PASS_DIR = 'pass'

# column width modifiers
COL_WIDTH_DICT = dict()
COL_WIDTH_DICT['MESSAGE'] = 100
COL_WIDTH_DICT['LOG_FILE'] = 300
COL_WIDTH_DICT['CODE'] = 20
COL_WIDTH_DICT['RUNSTRING'] = 200
COL_WIDTH_DICT['START_TIME'] = 50
COL_WIDTH_DICT['END_TIME'] = 50

# define the default column width
DEFAULT_COL_WIDTH = 30
# define table width info
DEFAULT_TABLE_WIDTH = 100
DEFAULT_TABLE_LENGTH = 6

# Currently takes ~ 1 minute for SPIROU full profile
SKIP_OBJ_TABLE = False
# Currently takes ~ 1 minute for SPIROU full profile
SKIP_RECIPE_TABLE = True
# Currently takes ~ 20 minute for SPIROU full profile
SKIP_MSG_TABLE = True


# =============================================================================
# Define functions
# =============================================================================
def get_settings(profile_name: Union[str, None] = None) -> dict:
    # storage for settings
    param_settings = dict()
    # clean profile name
    if profile_name is None:
        cpname = 'None'
    else:
        cpname = clean_profile_name(profile_name)
    # add the clean profile name to the settings
    param_settings['CPN'] = cpname
    # create paths
    param_settings['WORKING'] = _WORKING_DIR
    if profile_name is not None:
        param_settings['DATA'] = os.path.join(_WORKING_DIR, cpname, _DATA_DIR)
        param_settings['RST'] = os.path.join(_WORKING_DIR, cpname, _RST_DIR)
        param_settings['DIR'] = os.path.join(_WORKING_DIR, cpname)
        param_settings['PASS1'] = os.path.join(_WORKING_DIR, _PASS_DIR, '1',
                                               cpname)
        param_settings['PASS2'] = os.path.join(_WORKING_DIR, _PASS_DIR, '2',
                                               cpname)
    else:
        param_settings['DATA'] = None
        param_settings['RST'] = None
        param_settings['DIR'] = None
        param_settings['PASS1'] = None
        param_settings['PASS2'] = None
        # set the global paths for sphinx
    param_settings['HTML'] = os.path.join(_WORKING_DIR, _HTML_DIR)
    param_settings['OUT'] = os.path.join(_WORKING_DIR, _OUT_DIR)
    # make sure directories exist
    if profile_name is not None and not os.path.exists(param_settings['DATA']):
        os.makedirs(param_settings['DATA'])
    if profile_name is not None and not os.path.exists(param_settings['RST']):
        os.makedirs(param_settings['RST'])
    if profile_name is not None and not os.path.exists(param_settings['PASS1']):
        os.makedirs(param_settings['PASS1'])
    if profile_name is not None and not os.path.exists(param_settings['PASS2']):
        os.makedirs(param_settings['PASS2'])
    if not os.path.exists(param_settings['HTML']):
        os.makedirs(param_settings['HTML'])
    if not os.path.exists(param_settings['OUT']):
        os.makedirs(param_settings['OUT'])
    # return all parameter settings
    return param_settings


def compile_stats(settings: dict, profile: dict) -> dict:
    """
    Compile the stats for a given profile

    :param profile: dict, the profile to compile stats for
    :return: dict, the stats for the profile
    """
    # set up storage for the output dictionary
    profile_stats = dict()
    # deal with updating the path (as DRS_UCONFIG) from apero profile
    update_apero_profile(profile)
    # get paths to tables
    object_table_file = os.path.join(settings['DATA'], 'OBJECT_TABLE.fits')
    recipe_table_file = os.path.join(settings['DATA'], 'RECIPE_TABLE.fits')
    message_table_file = os.path.join(settings['DATA'], 'MESSAGE_TABLE.fits')
    # ------------------------------------------------------------------
    # deal with skipping object table
    if SKIP_OBJ_TABLE and os.path.exists(object_table_file):
        profile_stats['OBJECT_TABLE'] = Table.read(object_table_file)
    elif SKIP_OBJ_TABLE:
        profile_stats['OBJECT_TABLE'] = None
    else:
        # get the object table (astropy table)
        object_table = compile_apero_object_table()
        # add the lbl count
        object_table = add_lbl_count(profile, object_table)
        # add final object table to profile stats
        profile_stats['OBJECT_TABLE'] = object_table
    # ------------------------------------------------------------------
    # deal with skipping recipe table
    if SKIP_RECIPE_TABLE and os.path.exists(recipe_table_file):
        profile_stats['RECIPE_TABLE'] = Table.read(recipe_table_file)
    elif SKIP_RECIPE_TABLE:
        profile_stats['RECIPE_TABLE'] = None
    else:
        # get the recipe log table
        profile_stats['RECIPE_TABLE'] = compile_apero_recipe_table()
    # ------------------------------------------------------------------
    # deal with skipping message table
    if SKIP_MSG_TABLE and os.path.exists(message_table_file):
        profile_stats['MESSAGE_TABLE'] = Table.read(message_table_file)
    elif SKIP_MSG_TABLE:
        profile_stats['MESSAGE_TABLE'] = None
    else:
        # get the message log table
        profile_stats['MESSAGE_TABLE'] = compile_apero_message_table()
    # ------------------------------------------------------------------
    # return the profile stats
    return profile_stats


def update_apero_profile(profile: dict):
    """
    Update the apero profile with the correct paths
    :param profile: dict, the profile to update
    :return:
    """
    from apero.base import base
    from apero.core import constants
    from apero.core.constants import param_functions
    # use os to add DRS_UCONFIG to the path
    os.environ['DRS_UCONFIG'] = profile['apero profile']
    # reload DPARAMS and IPARAMS
    base.DPARAMS = base.load_database_yaml()
    base.IPARAMS = base.load_install_yaml()
    # ------------------------------------------------------------------
    # invalidate cache
    param_functions.CONFIG_CACHE = dict()
    # make sure parameters is reloaded (and not cached)
    _ = constants.load(cache=False)


def write_markdown(settings: dict, stats: dict):
    from apero.tools.module.documentation import drs_markdown
    # list tables to load
    tables = ['OBJECT_TABLE', 'RECIPE_TABLE', 'MESSAGE_TABLE']
    # -------------------------------------------------------------------------
    # step 1: write a page with all the different possible profiles in a
    #         sphinx table of contents
    # -------------------------------------------------------------------------
    profile_files = []
    # loop around profiles
    for profile_name in stats:
        # get settings
        settings = get_settings(profile_name)
        # reference name is just a cleaned version of the profile name
        ref_name = settings['CPN']
        # add to ref_names
        profile_files.append(f'{ref_name}/rst/profile.rst')
    # create ARI index page
    index_page = drs_markdown.MarkDownPage('ari_index')
    # add title
    index_page.add_title('APERO Reduction Interface (ARI)')
    # -------------------------------------------------------------------------
    # Add basic text
    # construct text to add
    index_page.add_text('This is the APERO Reduction Interface (ARI).')
    index_page.add_newline()
    index_page.add_text('Please select a reduction')
    index_page.add_newline()
    index_page.add_text('If you believe you should have the username/password '
                        'please contact neil.james.cook@gmail.com')
    # -------------------------------------------------------------------------
    # add table of contents
    index_page.add_table_of_contents(profile_files)
    # save index page
    index_page.write_page(os.path.join(settings['WORKING'], 'index.rst'))
    # -------------------------------------------------------------------------
    # step 2: write a page for each profile with the stats. Each page should
    #         just be a table of contents
    # -------------------------------------------------------------------------
    # loop around profiles
    for profile_name in stats:
        # get settings
        settings = get_settings(profile_name)
        # get the reference name
        cprofile_name = settings['CPN']
        # create a page
        profile_page = drs_markdown.MarkDownPage(cprofile_name)
        # add title
        profile_page.add_title(profile_name)
        # store the reference name for profile page table of contents
        table_files = []
        # loop around tables
        for table_name in tables:
            # get table
            table = stats[profile_name][table_name]
            # create a table page for this table
            table_ref_name = cprofile_name.lower() + '_' + table_name.lower()
            # make a markdown page for the table
            table_page = drs_markdown.MarkDownPage(table_ref_name)
            # add object table
            table_filename = f'{table_name}.csv'
            table_title = table_name.lower().replace('_', ' ')
            table_page.add_title(table_title)
            # deal with column widths for this file type
            if table is not None:
                # get column widths
                cwidth, cwidths = _get_column_widths(table)
                # add the csv version of this table
                table_page.add_csv_table('', '../data/' + table_filename,
                                         width=cwidth, widths=cwidths)
            else:
                # if we have no table then add a message
                table_page.add_text('No table created.')
            # save table page
            table_markdown_file = os.path.basename(table_filename).lower()
            table_markdown_file = table_markdown_file.replace('.csv', '.rst')
            # write table page
            print('Writing table page: {0}'.format(table_markdown_file))
            table_page.write_page(os.path.join(settings['RST'],
                                               table_markdown_file))
            # store the reference name for profile page table of contents
            #  these are in the same dir as profile_page so do not add the
            #  rst sub-dir
            table_files.append(os.path.basename(table_markdown_file))
        # add table of contents to profile page
        profile_page.add_table_of_contents(table_files)
        # save profile page
        profile_page.write_page(os.path.join(settings['RST'], 'profile.rst'))


def compile_docs(settings: dict):
    """
    Compile the documentation

    :param settings: dict, the settings for the documentation
    :return:
    """
    # must import here (so that os.environ is set)
    # noinspection PyPep8Naming
    from apero import lang
    from apero.core import constants
    from apero.core.core import drs_log
    from apero.core.utils import drs_startup
    from apero.io import drs_path
    # ------------------------------------------------------------------
    # get the parameter dictionary of constants from apero
    params = constants.load()
    # Get the text types
    textentry = lang.textentry
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # get WLOG
    wlog = drs_log.wlog
    # ------------------------------------------------------------------
    # get current directory
    cwd = os.getcwd()
    # change to docs directory
    os.chdir(settings['WORKING'])
    # ------------------------------------------------------------------
    # clean build
    # ------------------------------------------------------------------
    # Cleaning build directories
    wlog(params, '', textentry('40-506-00001'))
    os.system('make clean')
    # ------------------------------------------------------------------
    # build html
    # ------------------------------------------------------------------
    # Compiling HTML
    wlog(params, '', textentry('40-506-00002'))
    # Running make html
    wlog(params, '', textentry('40-506-00003'))
    # make html using sphinx
    os.system('make html')
    # ------------------------------------------------------------------
    # copy files to output directory
    # ------------------------------------------------------------------
    # clear out the output directory
    # Removing content of {0}
    wlog(params, '', textentry('40-506-00007', args=[settings['HTML']]))
    os.system(f'rm -rf {settings["OUT"]}/*')
    # copy
    drs_path.copytree(settings['HTML'], settings["OUT"])
    # ------------------------------------------------------------------
    # change back to current directory
    os.chdir(cwd)


def upload_docs(settings: dict):
    """
    Upload the documentation to the web server

    :param settings: dict, the settings for the documentation

    :return:
    """
    # must import here (so that os.environ is set)
    # noinspection PyPep8Naming
    from apero.core import constants
    from apero.tools.module.documentation import drs_documentation
    from apero.core.utils import drs_startup
    from apero.core.core import drs_log
    # ------------------------------------------------------------------
    # get the parameter dictionary of constants from apero
    params = constants.load()
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # get WLOG
    wlog = drs_log.wlog
    # get output directory
    out_dir = settings['OUT']
    # ------------------------------------------------------------------
    # change permission of all files and directories
    os.system(f'chmod 777 -R {out_dir}')
    # make sure we copy contents not directory
    if not out_dir.endswith(os.sep):
        out_dir += os.sep
    # get rsync dict
    rdict = dict()
    rdict['SSH'] = drs_documentation.SSH_OPTIONS
    rdict['USER'] = drs_documentation.SSH_USER
    rdict['HOST'] = drs_documentation.SSH_HOST
    rdict['INPATH'] = out_dir
    rdict['OUTPATH'] = os.path.join(SSH_PATH, 'ari/')
    # print command to rsync
    wlog(params, '', drs_documentation.RSYNC_CMD.format(**rdict))
    # run command (will require password)
    os.system(drs_documentation.RSYNC_CMD.format(**rdict))


def protect(profiles: dict):
    """
    Protect the documentation with .htaccess and .htpasswd files

    :param profiles: dict, the profiles to protect

    :return:
    """
    import getpass
    from apero.core import constants
    from apero.tools.module.documentation import drs_documentation
    from apero.core.utils import drs_startup
    from apero.core.core import drs_log
    # ------------------------------------------------------------------
    # get the parameter dictionary of constants from apero
    params = constants.load()
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # get WLOG
    wlog = drs_log.wlog
    # ------------------------------------------------------------------
    # password store
    password_store = dict()
    # create a .htpasswrd file for each profile
    for profile_name in profiles:
        # get the profile settings
        settings = get_settings(profile_name)
        # get the output directory
        cpn = settings['CPN']
        # get username from profile
        username = profiles[profile_name]['username']
        # get the password file
        local_pass_dir = os.path.join(settings['PASS1'])
        remote_pass_dir = HTPASSWD_PATH

        # check if account in profiles
        if 'account' in profiles[profile_name]:
            if len(profiles[profile_name]['account']) > 0:
                password_store[username] = profiles[profile_name]['account']

        # store password for same usernames
        if username in password_store:
            password = password_store[username]
        else:
            # try to get the password from the user and make sure it matches
            while True:
                # get the password
                password = getpass.getpass(f'Enter password for {profile_name}: ')
                # repeat password
                password2 = getpass.getpass(f'Repeat password for {profile_name}: ')
                if password == password2:
                    break
                else:
                    # Passwords do not match
                    print('Password must match')
        # create the password file
        os.system(f'htpasswd -b -c {local_pass_dir}/.htpasswd {username} {password}')

        # get the .htaccess file
        local_htaccess_file = os.path.join(settings['PASS2'], '.htaccess')
        remote_htaccess_file = os.path.join(SSH_PATH, 'ari/', cpn, '.htaccess')
        # construct lines to add to .htaccess file
        lines = []
        lines.append('AuthType Basic\n')
        lines.append('AuthName "Protected Site"\n')
        lines.append(f'AuthUserFile {remote_pass_dir}/{cpn}/.htpasswd\n')
        lines.append('require valid-user\n')
        # create the .htaccess file
        with open(local_htaccess_file, 'w') as htafile:
            htafile.writelines(lines)

        # copy local files to remote
        localfiles = [local_pass_dir, local_htaccess_file]
        remotefiles = [remote_pass_dir, remote_htaccess_file]
        # set up rsync dict
        rdict = dict()
        rdict['SSH'] = drs_documentation.SSH_OPTIONS
        rdict['USER'] = drs_documentation.SSH_USER
        rdict['HOST'] = drs_documentation.SSH_HOST
        # loop around files
        for it in range(len(localfiles)):
            # modify input/output rsync dict
            rdict['INPATH'] = localfiles[it]
            rdict['OUTPATH'] = remotefiles[it]
            # print command to rsync
            wlog(params, '', drs_documentation.RSYNC_CMD.format(**rdict))
            # run command (will require password)
            os.system(drs_documentation.RSYNC_CMD.format(**rdict))


# =============================================================================
# Functions for the apero stats
# =============================================================================
def compile_apero_object_table() -> Table:
    """
    Compile the apero object table

    :return: Table, the apero object table
    """
    # must import here (so that os.environ is set)
    # noinspection PyPep8Naming
    from apero.base.base import TQDM as tqdm
    from apero.core import constants
    from apero.core.core import drs_database
    from apero.core.core import drs_log
    from apero.core.utils import drs_startup
    # ------------------------------------------------------------------
    # get the parameter dictionary of constants from apero
    params = constants.load()
    # load pseudo-constants from apero
    pconst = constants.pload()
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # get WLOG
    wlog = drs_log.wlog
    # get polar condition
    has_polar = HAS_POLAR[params['INSTRUMENT']]
    # get the science fiber from pconst
    science_fibers, ref_fiber = pconst.FIBER_KINDS()
    # we assume the first science fiber is the primary science fiber
    science_fiber = science_fibers[0]
    # ------------------------------------------------------------------
    # get the astrometric database from apero
    astrodbm = drs_database.AstrometricDatabase(params)
    astrodbm.load_db()
    # ------------------------------------------------------------------
    # log that we are loading
    # get the object table from the astrometric database
    object_table = astrodbm.get_entries(columns=','.join(ASTROMETRIC_COLUMNS))
    # force columns to data types
    for itcol, col in enumerate(object_table):
        object_table[col] = np.array(object_table[col],
                                     dtype=ASTROMETRIC_DTYPES[itcol])
    # ------------------------------------------------------------------
    # get the index database from file index database
    indexdbm = drs_database.FileIndexDatabase(params)
    indexdbm.load_db()
    # ------------------------------------------------------------------
    # convert pandas dataframe to astropy table
    object_table = Table.from_pandas(object_table)
    # ------------------------------------------------------------------
    # add counting columns to the object table
    object_table['RAW_FILES'] = [0] * len(object_table)
    object_table['PP_FILES'] = [0] * len(object_table)
    object_table['EXT_FILES'] = [0] * len(object_table)
    object_table['TCORR_FILES'] = [0] * len(object_table)
    object_table['e.fits'] = [0] * len(object_table)
    object_table['t.fits'] = [0] * len(object_table)
    # deal with instruments that have polarimetry
    if has_polar:
        object_table['POLAR_FILES'] = [0] * len(object_table)
        object_table['p.fits'] = [0] * len(object_table)
    # ------------------------------------------------------------------
    # log progress
    wlog(params, '', 'Compiling object table (this may take a while)')
    # ------------------------------------------------------------------
    # for each object we run several counts
    # loop around objects in the object table
    for pos in tqdm(range(len(object_table))):
        # get the object name
        objname = object_table['OBJNAME'][pos]
        # set up where condition for raw files
        raw_cond = f'KW_OBJNAME="{objname}" AND BLOCK_KIND="raw"'
        # setup up where condition for pp files
        pp_cond = f'KW_OBJNAME="{objname}" AND BLOCK_KIND="tmp"'
        # setup up where condition for ext files
        # TODO: Get KW_OUTPUT from file definitions + deal with QC
        ext_cond = (f'KW_OBJNAME="{objname}" AND BLOCK_KIND="red" '
                    f'AND KW_OUTPUT="EXT_E2DS_FF" '
                    f'AND KW_FIBER="{science_fiber}"')
        # setup up where condition for tcorr files
        # TODO: Get KW_OUTPUT from file definitions + deal with QC
        tcorr_cond = (f'KW_OBJNAME="{objname}" AND BLOCK_KIND="red" '
                      f'AND KW_OUTPUT="TELLU_OBJ" '
                      f'AND KW_FIBER="{science_fiber}"')
        # setup up where condition for e.fits files
        e_cond = (f'KW_OBJNAME="{objname}" AND BLOCK_KIND="out" '
                  f'AND KW_OUTPUT="DRS_POST_E"')
        # setup up where condition for t.fits files
        t_cond = (f'KW_OBJNAME="{objname}" AND BLOCK_KIND="out" '
                  f'AND KW_OUTPUT="DRS_POST_T"')
        # deal with instruments that have polarimetry
        if has_polar:
            # setup up where condition for polar files
            polar_cond = (f'KW_OBJNAME="{objname}" AND BLOCK_KIND="red" '
                          f'AND KW_OUTPUT="POL_DEG"')
            # setup up where condition for p.fits files
            p_cond = (f'KW_OBJNAME="{objname}" AND BLOCK_KIND="out" '
                      f'AND KW_OUTPUT="DRS_POST_P"')
        else:
            polar_cond = None
            p_cond = None
        # ------------------------------------------------------------------
        # run counting conditions using indexdbm
        # ------------------------------------------------------------------
        conditions = [raw_cond, pp_cond, ext_cond, tcorr_cond, polar_cond,
                      e_cond, t_cond, p_cond]
        colnames = ['RAW_FILES', 'PP_FILES', 'EXT_FILES', 'TCORR_FILES',
                    'POLAR_FILES', 'e.fits', 't.fits', 'p.fits']
        # storage of count (for chains
        counts = dict()
        # define chains (if this number is zero do not count)
        chains = [None, 'RAW_FILES', 'PP_FILES', 'EXT_FILES',
                  'TCORR_FILES', 'PP_FILES', 'e.fits', 't.fits']
        # loop around conditions
        for it, condition in enumerate(conditions):
            # deal with polar conditions
            if condition is None:
                continue
            # deal with chains
            if chains[it] is not None:
                # get count from counts
                count = counts[chains[it]]
                # if the count is zero then continue
                if count == 0:
                    counts[colnames[it]] = 0
                    continue
            # run the count
            count = indexdbm.database.count(condition=condition)
            # add to object table
            object_table[colnames[it]][pos] = count
            # set counts
            counts[colnames[it]] = count
    # ------------------------------------------------------------------
    # remove rows with no raw entries
    mask = object_table['RAW_FILES'] > 0
    # apply mask
    object_table = object_table[mask]
    # ------------------------------------------------------------------
    # return object table
    return object_table


def compile_apero_recipe_table() -> Table:
    """
    Compile the apero recipe log table

    :return: apero recipe log table
    """
    # must import here (so that os.environ is set)
    from apero.core import constants
    from apero.core.core import drs_database
    from apero.core.core import drs_log
    from apero.core.utils import drs_startup
    # ------------------------------------------------------------------
    # get the parameter dictionary of constants from apero
    params = constants.load()
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # get WLOG
    wlog = drs_log.wlog
    # log progress
    wlog(params, '', 'Compiling apero log table (this may take a while)')
    # ------------------------------------------------------------------
    # get the log database from apero
    logdbm = drs_database.LogDatabase(params)
    logdbm.load_db()
    # get the log table using the LOG_OCLUMNS
    log_table = logdbm.database.get(columns=','.join(LOG_COLUMNS),
                                    sort_by='START_TIME',
                                    sort_descending=False,
                                    return_table=True)
    # ------------------------------------------------------------------
    return log_table


def compile_apero_message_table() -> Table:
    """
    Compile the apero message log table

    :return:
    """
    # TODO: Do this via database in future?
    # must import here (so that os.environ is set)
    # noinspection PyPep8Naming
    from apero.base.base import TQDM as tqdm
    from apero.core import constants
    from apero.core.core import drs_log
    from apero.core.utils import drs_startup
    # ------------------------------------------------------------------
    # get the parameter dictionary of constants from apero
    params = constants.load()
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # load pseudo-constants from apero
    pconst = constants.pload()
    # get WLOG
    wlog = drs_log.wlog
    # ------------------------------------------------------------------
    # get the log trigger characters
    log_trig_keys = pconst.LOG_TRIG_KEYS()
    # get the required log trig patterns
    log_trig_patterns = []
    for log_key in log_trig_keys:
        if log_key in LOG_LEVELS:
            log_trig_patterns.append(log_trig_keys[log_key])
    # get the pattern and the groups expected in a log message
    # Split format string into parts, consisting of literal text and
    # replacement fields
    # noinspection RegExpRedundantEscape
    log_parts = re.split(r'(\{[^}]*\})', params['DRS_LOG_FORMAT'])
    # ------------------------------------------------------------------
    # get the path to the log file directory
    log_dir = params['DRS_DATA_MSG']
    # get the list of log files
    # use os.walk to find all .log files in log_dir
    log_files = []
    for root, dirs, files in os.walk(log_dir):
        for file in files:
            if file.endswith(".log"):
                log_files.append(os.path.join(root, file))
    # ------------------------------------------------------------------
    # dictionary storage for appending to table
    table_dict = dict()
    # set table_dict columns
    table_dict['TIME'] = []
    table_dict['LEVEL'] = []
    table_dict['SHORTNAME'] = []
    table_dict['CODE'] = []
    table_dict['MESSAGE'] = []
    table_dict['LOG_FILE'] = []
    # log progress
    wlog(params, '', 'Compiling apero message table (this may take a while)')
    # loop around log files
    for log_file in tqdm(log_files):
        # read log file
        with open(log_file, 'r') as logfile:
            # read all lines
            lines = logfile.readlines()
            # -----------------------------------------------------------------
            # loop around lines
            for line in lines:
                # -------------------------------------------------------------
                # use the fmt to extract data from line
                # noinspection PyBroadException
                try:
                    variables = split_line(log_parts, line)
                except Exception as _:
                    continue
                # -------------------------------------------------------------
                # deal with no variables from this line
                if variables is None:
                    continue
                # now we can extracted the variables
                raw_time, raw_level, raw_id, raw_msg = variables
                # -------------------------------------------------------------
                # only continue if raw_level in required levels
                if raw_level not in log_trig_patterns:
                    continue
                # -------------------------------------------------------------
                # get time
                # TODO: update when we have real time (use Time class)
                time = raw_time
                # -------------------------------------------------------------
                # get level
                level = 'all'
                for key in log_trig_keys:
                    if raw_level == log_trig_keys[key]:
                        # add the key to LEVEL
                        level = key
                        # don't need to keep looking for the key
                        break
                # -------------------------------------------------------------
                # get the shortname
                shortname = raw_id
                # -------------------------------------------------------------
                # check to see whether '[' and ']' are in raw_msg
                try:
                    # get the code (we assume code is before the first "]")
                    code = raw_msg.split(']')[0].split('[')[1]
                    # get message (assume it is after "[code]:")
                    message = raw_msg.split(f'[{code}]:')[1]
                except Exception as _:
                    code = ''
                    message = raw_msg
                # -------------------------------------------------------------
                # clean message of escape characters
                message = message.replace('\\n', ' ')
                message = message.replace('\\t', ' ')
                message = message.replace('\\r', ' ')
                # -------------------------------------------------------------
                # finally append all data to table dict
                table_dict['TIME'].append(time)
                table_dict['LEVEL'].append(level)
                table_dict['SHORTNAME'].append(shortname)
                table_dict['CODE'].append(code)
                table_dict['MESSAGE'].append(message)
                table_dict['LOG_FILE'].append(log_file)
    # ------------------------------------------------------------------
    # convert table_dict to table
    message_table = Table(table_dict)
    # return table
    return message_table


# =============================================================================
# Functions for compiling the lbl stats
# =============================================================================
def add_lbl_count(profile: dict, object_table: Table) -> Table:
    # TODO: Fill out this function
    # must import here (so that os.environ is set)
    # noinspection PyPep8Naming
    from apero.base.base import TQDM as tqdm
    from apero.core import constants
    from apero.core.core import drs_log
    from apero.core.utils import drs_startup
    # get the parameter dictionary of constants from apero
    params = constants.load()
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # get WLOG
    wlog = drs_log.wlog
    # -------------------------------------------------------------------------
    # get the lbl path
    lbl_path = profile['lbl path']
    # -------------------------------------------------------------------------
    # we are adding three columns to the object table
    # first column: templates used for lbl
    # second column: template selected for lbl count
    # third column: number of lbl files for template with most lbl files
    object_table['LBL_TEMPLATES'] = np.zeros(len(object_table), dtype=str)
    object_table['LBL_SELECT'] = np.zeros(len(object_table), dtype=str)
    object_table['LBL_COUNT'] = np.zeros(len(object_table), dtype=int)
    lbl_templates = []
    lbl_select = []
    lbl_count = []
    # -------------------------------------------------------------------------
    # deal with no valid lbl path
    if lbl_path is None:
        return object_table
    # deal with lbl path not existing
    if not os.path.exists(lbl_path):
        return object_table
    # print that we are analysing lbl outputs
    wlog(params, '', 'Analysing LBL files')
    # -------------------------------------------------------------------------
    # get object name
    objnames = np.array(object_table['OBJNAME'])
    # loop around objects
    for pos, objname in tqdm(enumerate(objnames)):
        # for each object find all directories in lbl path that match this
        #   object name
        directories = glob.glob(os.path.join(lbl_path, 'lblrv', f'{objname}_*'))
        # deal with no directories --> skip
        if len(directories) == 0:
            lbl_templates.append('')
            lbl_select.append('')
            lbl_count.append(0)
            continue
        # store a list of templates
        templates = []
        # store a list of counts
        counts = []
        # loop around each directory
        for directory in directories:
            # get the template name for each directory
            template = os.path.basename(directory).split('_')[1]
            # get the number of lbl files in each directory
            count = len(glob.glob(os.path.join(directory, '*lbl.fits')))
            # append storage
            templates.append(template)
            counts.append(count)
        # decide which template to use (using max of counts)
        select = np.argmax(counts)
        lbl_templates.append(','.join(templates))
        lbl_select.append(templates[select])
        lbl_count.append(int(counts[select]))
    # -------------------------------------------------------------------------
    # add to object table
    object_table['LBL_TEMPLATES'] = lbl_templates
    object_table['LBL_SELECT'] = lbl_select
    object_table['LBL_COUNT'] = lbl_count
    # -------------------------------------------------------------------------
    # return the object table
    return object_table


# =============================================================================
# Worker functions
# =============================================================================
def split_line(parts, rawstring):
    # store list of variables
    variables = []
    number = 0
    # Split each replacement field into its parts
    for i, part in enumerate(parts):
        # skip empty parts
        if len(part) == 0:
            continue
        # skip variables
        if '{' in part and '}' in part:
            number += 1
            continue

        variable = rawstring.split(part)[0]
        # add to variables
        variables.append(variable)
        # update rawstring (remove variable and part)
        rawstring = rawstring.replace(variable + part, '', 1)
    # add the last variable as everything left in rawstring (unless empty)
    if len(rawstring) > 0:
        variables.append(rawstring)
    # return variables if equal to the number of {variables}
    if len(variables) == number:
        return variables
    # otherwise things went wrong
    else:
        return None


def save_stats(settings: dict, stats: dict):
    """
    Save the stats for a given profile

    :param settings: dict, the settings for the profile:
    :param stats: dict, the stats to save

    :return:
    """
    # list tables to load
    tables = ['OBJECT_TABLE', 'RECIPE_TABLE', 'MESSAGE_TABLE']
    # loop around tables and load them
    for table_name in tables:
        # deal with no table
        if stats[table_name] is None:
            continue
        # otherwise we assume it is an astropy table
        else:
            table = stats[table_name]
        # construct filename for table
        table_filename = f'{table_name}.fits'
        # construct full path
        table_path = os.path.join(settings['DATA'], table_filename)
        # write the table to the file
        table.write(table_path, format='fits',
                    overwrite=True)
        table.write(table_path.replace('.fits', '.csv'),
                    format='csv', overwrite=True)


def load_stats(settings: dict) -> dict:
    """
    Load the stats for a given profile

    :param settings: dict, the settings for the profile:

    :return: dict, the stats for the profile
    """
    # set up storage for the output dictionary
    profile_stats = dict()
    # list tables to load
    tables = ['OBJECT_TABLE', 'RECIPE_TABLE', 'MESSAGE_TABLE']
    # loop around tables and load them
    for table_name in tables:
        # construct filename for table
        table_filename = f'{table_name}.fits'
        # construct full path
        table_path = os.path.join(settings['DATA'], table_filename)
        # might not be on disk
        if not os.path.exists(table_path):
            profile_stats[table_name] = None
            continue
        # load the table
        profile_stats[table_name] = Table.read(table_path)
    # return the profile stats
    return profile_stats


def read_yaml_file() -> dict:
    """
    Read the yaml file PROFILE_FILE for the list of profiles using pyyaml

    :return:
    """
    # read the yaml file "profiles.yaml" for the list of profiles using pyyaml
    with open(PROFILE_FILE, 'r') as stream:
        # try to load the yaml file
        try:
            profiles = yaml.safe_load(stream)
        # deal with a yaml error (for now just raise the error)
        except yaml.YAMLError as exc:
            raise exc
    # return the profiles
    return profiles


def clean_profile_name(profile_name: str) -> str:
    """
    Clean up the profile name

    :param profile_name: str, the profile name
    :return: str, the cleaned up profile name
    """
    # clean up profile name
    profile_name = profile_name.replace(' ', '_')
    profile_name = profile_name.replace('-', '_')
    profile_name = profile_name.replace('(', '_')
    profile_name = profile_name.replace(')', '_')
    # remove ugly double underscores
    while '__' in profile_name:
        profile_name = profile_name.replace('__', '_')
    # remove starting and ending underscores
    profile_name = profile_name.strip('_')
    # return cleaned up profile name
    return profile_name


def _get_column_widths(table: Table):
    """
    Take a table and get columns widths from lookup table
    (or assign default value)
    """
    cwidths = []
    # loop around column names and look up the widths in the lookup table
    for colname in table.colnames:
        # if they are in the look up table use this width
        if colname in COL_WIDTH_DICT:
            cwidths.append(COL_WIDTH_DICT[colname])
        # otherwise use the default width
        else:
            cwidths.append(DEFAULT_COL_WIDTH)
    # widths must be percentages (100% total)
    cwidths = np.array(cwidths)
    cwidths = np.floor(100 * cwidths / np.sum(cwidths)).astype(int) - 1
    # widths must be strings
    cwidths = list(cwidths.astype(str))
    # -------------------------------------------------------------------------
    # deal with table width
    if len(table.colnames) <= DEFAULT_TABLE_LENGTH:
        cwidth = None
    else:
        cfrac = DEFAULT_TABLE_WIDTH / DEFAULT_TABLE_LENGTH
        cwidth = '{0}%'.format(int(np.floor(len(table.colnames) * cfrac)))
    # -------------------------------------------------------------------------
    # return a list of the columns
    return cwidth, cwidths


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # step 1: read the yaml file
    apero_profiles = read_yaml_file()
    # ----------------------------------------------------------------------
    # step 2: for each profile compile all stats
    all_apero_stats = dict()
    # loop around profiles from yaml file
    for apero_profile_name in apero_profiles:
        # sort out settings
        ari_settings = get_settings(apero_profile_name)
        # we reprocess if the file does not exist or if REPROCESS is True
        if REPROCESS:
            # print progress
            print('=' * 50)
            print('Compiling stats for profile: {0}'.format(apero_profile_name))
            print('=' * 50)
            # get profile
            apero_profile = apero_profiles[apero_profile_name]
            # compile stats
            apero_stats = compile_stats(ari_settings, apero_profile)
            # add to all_apero_stats
            all_apero_stats[apero_profile_name] = apero_stats
            # -----------------------------------------------------------------
            # Save stats to disk
            save_stats(ari_settings, apero_stats)
        # otherwise we load the stats from disk
        else:
            # print progress
            print('=' * 50)
            print('Loading stats for profile: {0}'.format(apero_profile_name))
            print('=' * 50)
            # load stats
            apero_stats = load_stats(ari_settings)
            # add to all_apero_stats
            all_apero_stats[apero_profile_name] = apero_stats
    # ----------------------------------------------------------------------
    # sort out settings
    ari_settings = get_settings()
    # ----------------------------------------------------------------------
    # step 3: write markdown files
    write_markdown(ari_settings, all_apero_stats)
    # ----------------------------------------------------------------------
    # step 4: compile sphinx files
    compile_docs(ari_settings)
    # ----------------------------------------------------------------------
    # step 5: upload to hosting
    upload_docs(ari_settings)
    # ----------------------------------------------------------------------
    # step 6: protect profiles
    protect(apero_profiles)

# =============================================================================
