#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-02-20 at 9:35

@author: cook
"""
import os
import re

import numpy as np
import pandas as pd
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
WORKING_DIR = 'working'
DATA_DIR = 'working/data/'
RST_DIR = 'working/rst/'
HTML_DIR = 'working/_build/html/'
OUT_DIR = 'working/output/'
SSH_PATH = '/export/www/home/cook/www/apero-drs/'

# =============================================================================
# Define functions
# =============================================================================
def compile_stats(profile: dict) -> dict:
    """
    Compile the stats for a given profile

    :param profile: dict, the profile to compile stats for
    :return: dict, the stats for the profile
    """
    # set up storage for the output dictionary
    profile_stats = dict()
    # deal with updating the path (as DRS_UCONFIG) from apero profile
    update_apero_profile(profile)
    # get the object table (astropy table)
    object_table = compile_apero_object_table()
    # add the lbl count
    object_table = add_lbl_count(profile, object_table)
    # add final object table to profile stats
    profile_stats['OBJECT_TABLE'] = object_table
    # get the recipe log table
    profile_stats['RECIPE_TABLE'] = compile_apero_recipe_table()
    # get the message log table
    profile_stats['MESSAGE_TABLE'] = compile_apero_message_table()
    # return the profile stats
    return profile_stats


def update_apero_profile(profile: dict):
    """
    Update the apero profile with the correct paths
    :param profile: dict, the profile to update
    :return:
    """
    # use os to add DRS_UCONFIG to the path
    os.environ['DRS_UCONFIG'] = profile['apero profile']


def write_markdown(stats: dict):
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
        # reference name is just a cleaned version of the profile name
        ref_name = clean_profile_name(profile_name)
        # add to ref_names
        profile_files.append('rst/profile_' + ref_name + '.rst')
    # create ARI index page
    index_page = drs_markdown.MarkDownPage('ari_index')
    # add title
    index_page.add_title('APERO Reduction Interface (ARI)')
    # -------------------------------------------------------------------------
    # add table of contents
    index_page.add_table_of_contents(profile_files)
    # save index page
    index_page.write_page(os.path.join(WORKING_DIR, 'index.rst'))
    # -------------------------------------------------------------------------
    # step 2: write a page for each profile with the stats. Each page should
    #         just be a table of contents
    # -------------------------------------------------------------------------
    # loop around profiles
    for profile_name in stats:
        # get the reference name
        ref_name = clean_profile_name(profile_name)
        # create a page
        profile_page = drs_markdown.MarkDownPage(ref_name)
        # add title
        profile_page.add_title(profile_name)
        # store the reference name for profile page table of contents
        table_files = []
        # loop around tables
        for table_name in tables:
            # create a table page for this table
            table_ref_name = ref_name.lower() + '_' + table_name.lower()
            # make a markdown page for the table
            table_page = drs_markdown.MarkDownPage(table_ref_name)
            # add object table
            table_filename = construct_table_filename(profile_name, table_name,
                                                      fmt='csv')
            table_title = table_name.lower().replace('_', ' ')
            table_page.add_csv_table(table_title, table_filename)
            # save table page
            table_markdown_file = os.path.basename(table_filename).lower()
            table_markdown_file = table_markdown_file.replace('.csv', '.rst')
            # write table page
            print('Writing table page: {0}'.format(table_markdown_file))
            table_page.write_page(os.path.join(RST_DIR, table_markdown_file))
            # store the reference name for profile page table of contents
            table_files.append('rst/' + os.path.basename(table_markdown_file))
        # add table of contents to profile page
        profile_page.add_table_of_contents(table_files)
        # save profile page
        profile_file = 'profile_{0}.rst'.format(ref_name)
        profile_page.write_page(os.path.join(RST_DIR, profile_file))


def compile_docs():
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
    os.chdir(os.path.realpath(WORKING_DIR))
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
    wlog(params, '', textentry('40-506-00007', args=[HTML_DIR]))
    os.system(f'rm -rf {os.path.realpath(OUT_DIR)}/*')
    # copy
    drs_path.copytree(os.path.realpath(HTML_DIR), os.path.realpath(OUT_DIR))
    # ------------------------------------------------------------------
    # change back to current directory
    os.chdir(cwd)


def upload_docs():
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
    out_dir = os.path.realpath(OUT_DIR)
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


# =============================================================================
# Functions for the apero stats
# =============================================================================
def compile_apero_object_table() -> pd.DataFrame:
    """
    Compile the apero object table

    :return: pd.DataFrame, the object table
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
    # add counting columns to the object table
    object_table['RAW_FILES'] = [0] * len(object_table)
    object_table['PP_FILES'] = [0] * len(object_table)
    object_table['EXT_FILES'] = [0] * len(object_table)
    object_table['TCORR_FILES'] = [0] * len(object_table)
    object_table['e.fits'] = [0] * len(object_table)
    object_table['t.fits'] = [0] * len(object_table)
    # deal with instruments that have polarimetry
    if HAS_POLAR[params['INSTRUMENT']]:
        object_table['POLAR_FILES'] = [0] * len(object_table)
        object_table['p.fits'] = [0] * len(object_table)
    # ------------------------------------------------------------------
    # log progress
    wlog(params, '', 'Compiling object table (this may take a while)')
    # ------------------------------------------------------------------
    # for each object we run several counts
    # loop around objects in the object table
    for pos in tqdm(object_table.index):
        # get the object name
        objname = object_table.loc[pos, 'OBJNAME']
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
        if HAS_POLAR:
            # setup up where condition for polar files
            polar_cond = (f'KW_OBJNAME = "{objname}" AND BLOCK_KIND="red" '
                          f'AND KW_OUTPUT="POL_DEG"')
            # setup up where condition for p.fits files
            p_cond = (f'KW_OBJNAME = "{objname}" AND BLOCK_KIND="out" '
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
            object_table.at[pos, colnames[it]] = count
            # set counts
            counts[colnames[it]] = count
    # ------------------------------------------------------------------
    # return object table
    return object_table


def compile_apero_recipe_table() -> pd.DataFrame:
    """
    Compile the apero recipe log table

    :return: pd.DataFrame, the recipe log table
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
    wlog(params, '', 'Compiling apero log table')
    # ------------------------------------------------------------------
    # get the log database from apero
    logdbm = drs_database.LogDatabase(params)
    logdbm.load_db()
    # get the log table using the LOG_OCLUMNS
    log_table = logdbm.database.get(columns=','.join(LOG_COLUMNS),
                                    sort_by='START_TIME',
                                    sort_descending=False,
                                    return_pandas=True)
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
                if '[' in raw_msg and ']' in raw_msg:
                    # get the code (we assume code is before the first "]")
                    code = raw_msg.split(']')[0].split('[')[1]
                    # get message (assume it is after "[code]:")
                    message = raw_msg.split(f'[{code}]:')[1]
                else:
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
def add_lbl_count(profile: dict, object_table: pd.DataFrame) -> pd.DataFrame:
    # TODO: Fill out this function
    _ = profile
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


def save_stats(stats: dict, profile_name: str):
    """
    Save the stats for a given profile

    :param stats: dict, the stats to save
    :param stats_filename: str, the filename to save the stats to

    :return:
    """
    # list tables to load
    tables = ['OBJECT_TABLE', 'RECIPE_TABLE', 'MESSAGE_TABLE']
    # loop around tables and load them
    for table_name in tables:
        # deal with pandas dataframes
        if isinstance(stats[table_name], pd.DataFrame):
            table = Table.from_pandas(stats[table_name])
        # otherwise we assume it is an astropy table
        else:
            table = stats[table_name]
        # construct filename for table
        table_filename = construct_table_filename(profile_name, table_name)
        # construct full path
        table_path = os.path.join(OUTPUT_PATH, DATA_DIR, 'data', table_filename)
        # write the table to the file
        table.write(table_path, format='fits',
                    overwrite=True)
        table.write(table_path.replace('.fits', '.csv'),
                    format='csv', overwrite=True)


def load_stats(profile_name: str) -> dict:
    """
    Load the stats for a given profile

    :param stats_filename: str, the filename to load the stats from
    :return: dict, the stats for the profile
    """
    # set up storage for the output dictionary
    profile_stats = dict()
    # list tables to load
    tables = ['OBJECT_TABLE', 'RECIPE_TABLE', 'MESSAGE_TABLE']
    # loop around tables and load them
    for table_name in tables:
        # construct filename for table
        table_filename = construct_table_filename(profile_name, table_name)
        # construct full path
        table_path = os.path.join(OUTPUT_PATH, DATA_DIR, table_filename)
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
    # return cleaned up profile name
    return profile_name


def construct_filename(profile_name: str) -> str:
    """
    Construct the filename for the stats file

    :param profile_name: str, the name of the profile
    :return: str, the apero stats filename
    """
    # clean up profile name
    profile_name = clean_profile_name(profile_name)
    # return cleaned up profile name
    return 'apero_ari_{0}.fits'.format(profile_name)


def construct_table_filename(profile_name: str, table_name: str,
                             fmt='fits') -> str:
    """
    Construct the filename for the stats file

    :param profile_name: str, the name of the profile
    :param table_name: str, the name of the table
    :return: str, the apero stats filename
    """
    # clean up profile name
    profile_name = clean_profile_name(profile_name)
    # construct filename
    filename = 'apero_ari_{0}_{1}.{2}'.format(profile_name, table_name, fmt)
    # return filename
    return filename


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # step 1: read the yaml file
    apero_profiles = read_yaml_file()
    # make output path
    outpath = os.path.join(OUTPUT_PATH, WORKING_DIR)
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)
    if not os.path.exists(RST_DIR):
        os.makedirs(RST_DIR)
    # make out direcory
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    # ----------------------------------------------------------------------
    # step 2: for each profile compile all stats
    all_apero_stats = dict()
    # loop around profiles from yaml file
    for apero_profile_name in apero_profiles:
        # construct filename for apero stats
        ari_filename = construct_filename(apero_profile_name)
        apero_stats_filename = os.path.join(outpath, ari_filename)
        # we reprocess if the file does not exist or if REPROCESS is True
        if REPROCESS:
            # print progress
            print('=' * 50)
            print('Compiling stats for profile: {0}'.format(apero_profile_name))
            print('=' * 50)
            # get profile
            apero_profile = apero_profiles[apero_profile_name]
            # compile stats
            apero_stats = compile_stats(apero_profile)
            # add to all_apero_stats
            all_apero_stats[apero_profile_name] = apero_stats
            # -----------------------------------------------------------------
            # Save stats to disk
            save_stats(apero_stats, apero_profile_name)
        # otherwise we load the stats from disk
        else:
            # print progress
            print('=' * 50)
            print('Loading stats for profile: {0}'.format(apero_profile_name))
            print('=' * 50)
            # load stats
            apero_stats = load_stats(apero_profile_name)
            # add to all_apero_stats
            all_apero_stats[apero_profile_name] = apero_stats
    # ----------------------------------------------------------------------
    # step 3: write markdown files
    write_markdown(all_apero_stats)
    # ----------------------------------------------------------------------
    # step 4: compile sphinx files
    compile_docs()
    # ----------------------------------------------------------------------
    # step 5: upload to hosting
    upload_docs()

# =============================================================================
