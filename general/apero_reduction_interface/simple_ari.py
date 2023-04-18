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
import shutil
import sys
from typing import Any, Dict, List, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import yaml
from astropy.table import Table
from astropy.time import Time
from astropy.io import fits
from astropy import units as uu

# =============================================================================
# Define variables
# =============================================================================
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
_DATA_DIR = 'data'
_RST_DIR = 'rst'
_HTML_DIR = '_build/html'
_OUT_DIR = 'output'
_ITEM_DIR = 'items'
_DOWN_DIR = 'downloads'
# Define output object data dir
_OBJ_OUT_DIR = 'objects'
# define path to local htpasswd file
_PASS_DIR = 'pass'
# define the column which is the object name
OBJECT_COLUMN = 'OBJNAME'
# define the count / file columns
COUNT_COLS = ['RAW', 'PP', 'EXT', 'TCORR', 'CCF', 'POL', 'e', 't', 'v', 'p']
# define chains (if this number is zero do not count)
COUNT_CHAINS = [None, 'RAW', 'PP', 'EXT', 'TCORR', 'TCORR', 'PP', 'EXT',
                'TCORR', 'TCORR']
# define the lbl rdb suffix (lbl or lbl2)
LBL_SUFFIX = 'lbl'
# object page styling
DIVIDER_COLOR = '#FFA500'
DIVIDER_HEIGHT = 6
# column width modifiers
COL_WIDTH_DICT = dict()
COL_WIDTH_DICT['MESSAGE'] = 100
COL_WIDTH_DICT['LOG_FILE'] = 300
COL_WIDTH_DICT['CODE'] = 20
COL_WIDTH_DICT['RUNSTRING'] = 200
COL_WIDTH_DICT['START_TIME'] = 60
COL_WIDTH_DICT['END_TIME'] = 60
COL_WIDTH_DICT['RA_DEG'] = 3
COL_WIDTH_DICT['DEC_DEG'] = 3
COL_WIDTH_DICT['TEFF'] = 3
COL_WIDTH_DICT['SP_TYPE'] = 3
COL_WIDTH_DICT['RAW'] = 3
COL_WIDTH_DICT['PP'] = 3
COL_WIDTH_DICT['EXT'] = 3
COL_WIDTH_DICT['TCORR'] = 3
COL_WIDTH_DICT['e'] = 3
COL_WIDTH_DICT['t'] = 3
COL_WIDTH_DICT['POL'] = 3
COL_WIDTH_DICT['p'] = 3
COL_WIDTH_DICT['LBL_COUNT'] = 3
COL_WIDTH_DICT['LBL_SELECT'] = 3
# define the default column width
DEFAULT_COL_WIDTH = 10
# define table width info
DEFAULT_TABLE_WIDTH = 100
DEFAULT_TABLE_LENGTH = 6


# =============================================================================
# Define functions
# =============================================================================
def get_settings(settings: Dict[str, Any],
                 profile_name: Union[str, None] = None) -> dict:
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
    working_dir = settings['working directory']
    param_settings['WORKING'] = working_dir
    if profile_name is not None:
        param_settings['DATA'] = os.path.join(working_dir, cpname, _DATA_DIR)
        param_settings['RST'] = os.path.join(working_dir, cpname, _RST_DIR)
        param_settings['DIR'] = os.path.join(working_dir, cpname)
        param_settings['PASS1'] = os.path.join(working_dir, _PASS_DIR, '1',
                                               cpname)
        param_settings['PASS2'] = os.path.join(working_dir, _PASS_DIR, '2',
                                               cpname)
        param_settings['ITEMS'] = os.path.join(working_dir, cpname, _ITEM_DIR)
        param_settings['DOWNS'] = os.path.join(working_dir, cpname, _DOWN_DIR)
        param_settings['OBJ_OUT'] = os.path.join(working_dir, cpname,
                                                 _OBJ_OUT_DIR)
    else:
        param_settings['DATA'] = None
        param_settings['RST'] = None
        param_settings['DIR'] = None
        param_settings['PASS1'] = None
        param_settings['PASS2'] = None
        param_settings['ITEMS'] = None
        param_settings['DOWNS'] = None
        param_settings['OBJ_OUT'] = None
    # set the global paths for sphinx
    param_settings['HTML'] = os.path.join(working_dir, _HTML_DIR)
    param_settings['OUT'] = os.path.join(working_dir, _OUT_DIR)
    # get params from apero
    param_settings['INSTRUMENT'] = settings['instrument']
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
    if not os.path.exists(param_settings['OBJ_OUT']):
        os.makedirs(param_settings['OBJ_OUT'])
    if not os.path.exists(param_settings['ITEMS']):
        os.makedirs(param_settings['ITEMS'])
    if not os.path.exists(param_settings['DOWNS']):
        os.makedirs(param_settings['DOWNS'])
    # return all parameter settings
    return param_settings


def compile_stats(settings: dict, profile: dict, headers: dict) -> dict:
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
    # get skip criteria
    skip_obj_table = profile['skip obj table']
    skip_recipe_table = profile['skip recipe table']
    skip_msg_table = profile['skip msg table']
    # ------------------------------------------------------------------
    # deal with skipping object table
    if skip_obj_table and os.path.exists(object_table_file):
        profile_stats['OBJECT_TABLE'] = Table.read(object_table_file)
    elif skip_obj_table:
        profile_stats['OBJECT_TABLE'] = None
    else:
        # get the object table (astropy table)
        object_table, filedict = compile_apero_object_table()
        # add the lbl count
        object_table, filedict = add_lbl_count(profile, object_table, filedict)
        # add the object pages
        add_obj_pages(settings, profile, headers, object_table, filedict)
        # remove the temporary object column
        del object_table['_OBJECT']
        # add final object table to profile stats
        profile_stats['OBJECT_TABLE'] = object_table
    # ------------------------------------------------------------------
    # deal with skipping recipe table
    if skip_recipe_table and os.path.exists(recipe_table_file):
        profile_stats['RECIPE_TABLE'] = Table.read(recipe_table_file)
    elif skip_recipe_table:
        profile_stats['RECIPE_TABLE'] = None
    else:
        # get the recipe log table
        profile_stats['RECIPE_TABLE'] = compile_apero_recipe_table()
    # ------------------------------------------------------------------
    # deal with skipping message table
    if skip_msg_table and os.path.exists(message_table_file):
        profile_stats['MESSAGE_TABLE'] = Table.read(message_table_file)
    elif skip_msg_table:
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


def write_markdown(gsettings: dict, settings: dict, stats: dict):
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
        settings = get_settings(gsettings, profile_name)
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
    index_page.add_newline()
    index_page.add_text('Last updated: {0} [UTC]'.format(Time.now()))
    index_page.add_newline()
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
        settings = get_settings(gsettings, profile_name)
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
                # add the csv version of this table
                table_page.add_csv_table('', f'../{_DATA_DIR}/' +
                                         table_filename)
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
    # get list of files to copy
    copy_files = ['conf.py', 'make.bat', 'Makefile']
    # loop around files and copy
    for copy_file in copy_files:
        # copy conf.py make.bat and Makefile to the working directory
        shutil.copy(__file__.replace('simple_ari.py', copy_file),
                    os.path.join(settings['WORKING'], copy_file))
    # get _static directory
    static_outdir = os.path.join(settings['WORKING'], '_static')
    # deal with static_outdir existing
    if os.path.exists(static_outdir):
        shutil.rmtree(static_outdir)
    # copy the static directory as well
    shutil.copytree(__file__.replace('simple_ari.py', '_static'),
                    os.path.join(settings['WORKING'], '_static'))
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


def upload_docs(gsettings: dict, settings: dict):
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
    # get the instrument
    instrument = settings['INSTRUMENT'].lower()
    # ------------------------------------------------------------------
    # change permission of all files and directories
    os.system(f'chmod 777 -R {out_dir}')
    # make sure we copy contents not directory
    if not out_dir.endswith(os.sep):
        out_dir += os.sep
    # get rsync dict
    rdict = dict()
    rdict['SSH'] = gsettings['ssh']['options']
    rdict['USER'] = gsettings['ssh']['user']
    rdict['HOST'] = gsettings['ssh']['host']
    rdict['INPATH'] = out_dir
    rdict['OUTPATH'] = os.path.join(gsettings['ssh']['directory'],
                                    f'ari/{instrument}/')
    # print command to rsync
    wlog(params, '', drs_documentation.RSYNC_CMD.format(**rdict))
    # run command (will require password)
    os.system(drs_documentation.RSYNC_CMD.format(**rdict))


# =============================================================================
# Functions for the apero stats
# =============================================================================
FileDictReturn = Dict[str, Dict[str, List[str]]]


def compile_apero_object_table() -> Tuple[Table, FileDictReturn]:
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
    # Deal with object pages
    # ------------------------------------------------------------------
    # add a temporary column
    object_table['_OBJECT'] = np.array(object_table[OBJECT_COLUMN])
    object_table[OBJECT_COLUMN] = [' ' * 255] * len(object_table)
    # ------------------------------------------------------------------
    # add counting columns to the object table
    object_table['RAW'] = [0] * len(object_table)
    object_table['LATEST_RAW'] = [' ' * 25] * len(object_table)
    object_table['PP'] = [0] * len(object_table)
    object_table['EXT'] = [0] * len(object_table)
    object_table['TCORR'] = [0] * len(object_table)
    # deal with instruments that have polarimetry
    if has_polar:
        object_table['POL'] = [0] * len(object_table)
    object_table['e'] = [0] * len(object_table)
    object_table['t'] = [0] * len(object_table)
    # deal with instruments that have polarimetry
    if has_polar:
        object_table['p'] = [0] * len(object_table)
    # ------------------------------------------------------------------
    # storage for files for each type
    file_dict = dict()
    # loop around column names and add a empty list for each
    for objname in object_table[OBJECT_COLUMN]:
        file_dict[objname] = dict()
        for colname in COUNT_COLS:
            file_dict[objname][colname] = []
    # ------------------------------------------------------------------
    # log progress
    wlog(params, '', 'Compiling object table (this may take a while)')
    # ------------------------------------------------------------------
    # for each object we run several counts
    # loop around objects in the object table
    for pos in tqdm(range(len(object_table))):
        # get the object name
        objname = object_table['_OBJECT'][pos]
        # set up where condition for raw files
        raw_cond = f'KW_OBJNAME="{objname}" AND BLOCK_KIND="raw"'
        # setup up where condition for pp files
        pp_cond = f'KW_OBJNAME="{objname}" AND BLOCK_KIND="tmp"'
        # setup where condition for ext files
        # TODO: Get KW_OUTPUT from file definitions + deal with QC
        ext_cond = (f'KW_OBJNAME="{objname}" AND BLOCK_KIND="red" '
                    f'AND KW_OUTPUT="EXT_E2DS_FF" '
                    f'AND KW_FIBER="{science_fiber}"')
        # setup where condition for tcorr files
        # TODO: Get KW_OUTPUT from file definitions + deal with QC
        tcorr_cond = (f'KW_OBJNAME="{objname}" AND BLOCK_KIND="red" '
                      f'AND KW_OUTPUT="TELLU_OBJ" '
                      f'AND KW_FIBER="{science_fiber}"')
        # setup where condition for ccf files
        # TODO: Get KW_OUTPUT from file definitions + deal with QC
        ccf_cond = (f'KW_OBJNAME="{objname}" AND BLOCK_KIND="red" '
                      f'AND KW_OUTPUT="CCF_RV" '
                      f'AND KW_FIBER="{science_fiber}"')
        # setup where condition for e.fits files
        e_cond = (f'KW_OBJNAME="{objname}" AND BLOCK_KIND="out" '
                  f'AND KW_OUTPUT="DRS_POST_E"')
        # setup where condition for t.fits files
        t_cond = (f'KW_OBJNAME="{objname}" AND BLOCK_KIND="out" '
                  f'AND KW_OUTPUT="DRS_POST_T"')
        # setup where condition for v.fits files
        v_cond = (f'KW_OBJNAME="{objname}" AND BLOCK_KIND="out" '
                  f'AND KW_OUTPUT="DRS_POST_V"')
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
        # deal with finding the latest raw file
        # ------------------------------------------------------------------
        latest_cond = f'KW_OBJNAME="{objname}" AND BLOCK_KIND="raw"'
        times = indexdbm.get_entries('KW_MID_OBS_TIME', condition=latest_cond)
        # if there are no entries we have no raw files for this object
        if len(times) == 0:
            continue
        # find maximum time value and convert to human time
        latest_time = Time(np.max(times), format='mjd')
        object_table['LATEST_RAW'][pos] = latest_time.iso
        # ------------------------------------------------------------------
        # run counting conditions using indexdbm
        # ------------------------------------------------------------------
        conditions = [raw_cond, pp_cond, ext_cond, tcorr_cond, ccf_cond,
                      polar_cond, e_cond, t_cond, v_cond, p_cond]
        # storage of count (for chains
        counts = dict()
        # loop around conditions
        for it, condition in enumerate(conditions):
            # deal with polar conditions
            if condition is None:
                continue
            # deal with chains
            if COUNT_CHAINS[it] is not None:
                # get count from counts
                count = counts[COUNT_CHAINS[it]]
                # if the count is zero then continue
                if count == 0:
                    counts[COUNT_COLS[it]] = 0
                    continue
            # get the files
            files = indexdbm.get_entries('ABSPATH', condition=condition)
            # get the count
            count = len(files)
            # append to files
            file_dict[objname][COUNT_COLS[it]] = files
            # add to object table
            object_table[COUNT_COLS[it]][pos] = count
            # set counts
            counts[COUNT_COLS[it]] = count
    # ------------------------------------------------------------------
    # remove rows with no raw entries
    mask = object_table['RAW'] > 0
    # apply mask
    object_table = object_table[mask]
    # ------------------------------------------------------------------
    # remove objects for file_dict with zero raw entries
    for objname in object_table[OBJECT_COLUMN]:
        # deal with objects not in the file dict
        if objname not in file_dict:
            continue
        # deal with objects with no raw files
        if len(file_dict[objname]['RAW']) == 0:
            del file_dict[objname]
    # ------------------------------------------------------------------
    # return object table
    return object_table, file_dict


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
def add_lbl_count(profile: dict, object_table: Table,
                  file_dict: FileDictReturn) -> Tuple[Table, FileDictReturn]:
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
    lbl_templates = []
    lbl_select = []
    lbl_count = []
    # -------------------------------------------------------------------------
    # deal with no valid lbl path
    if lbl_path is None:
        return object_table, file_dict
    # deal with lbl path not existing
    if not os.path.exists(lbl_path):
        return object_table, file_dict
    # print that we are analysing lbl outputs
    wlog(params, '', 'Analysing LBL files')
    # -------------------------------------------------------------------------
    # get object name
    objnames = np.array(object_table['_OBJECT'])
    # loop around objects
    for pos, objname in tqdm(enumerate(objnames)):
        # for each object find all directories in lbl path that match this
        #   object name
        lblrv_dir = glob.glob(os.path.join(lbl_path, 'lblrv', f'{objname}_*'))
        # get all the lbl rdb files for this object name
        rdb_glob = f'{LBL_SUFFIX}_{objname}_*.rdb'
        all_lblrdb_files = glob.glob(os.path.join(lbl_path, 'lblrdb', rdb_glob))
        # remove drift files from the lbl rdb files
        lblrdb_files = []
        for lblrdb_file in all_lblrdb_files:
            if 'drift' not in lblrdb_file:
                lblrdb_files.append(lblrdb_file)
        # ---------------------------------------------------------------------
        # deal with no directories --> skip
        if len(lblrv_dir) == 0:
            lbl_templates.append('')
            lbl_select.append('')
            lbl_count.append(0)
            # add an empty list to the LBLRDB file dict for this object
            file_dict[objname]['LBLRDB'] = []
            continue
        # store a list of templates
        templates = []
        # store a list of counts
        counts = []
        # loop around each directory
        for directory in lblrv_dir:
            # get the template name for each directory
            basename = os.path.basename(directory)
            template = basename.split(f'{objname}_')[-1]
            # get the number of lbl files in each directory
            count = len(glob.glob(os.path.join(directory, '*lbl.fits')))
            # append storage
            templates.append(template)
            counts.append(count)
        # decide which template to use (using max of counts)
        select = np.argmax(counts)
        # get strings to add to storage
        _template = ','.join(templates)
        _select = templates[select]
        _count = int(counts[select])
        # append to lists
        lbl_templates.append(_template)
        lbl_select.append(_select)
        lbl_count.append(_count)
        # add an empty list to the LBLRDB file dict for this object
        file_dict[objname]['LBLRDB'] = lblrdb_files
    # -------------------------------------------------------------------------
    # add to object table
    object_table['LBL_TEMPLATES'] = lbl_templates
    object_table['LBL_SELECT'] = lbl_select
    object_table['LBL_COUNT'] = lbl_count
    # -------------------------------------------------------------------------
    # return the object table
    return object_table, file_dict


# =============================================================================
# Functions for making object pages
# =============================================================================
class ObjectData:
    def __init__(self, profile, settings, headers,
                 objname, file_dict, object_table):
        # get the object name
        self.objname = objname
        # get the file_dictionary for this object
        self.raw_files = file_dict[objname]['RAW']
        self.pp_files = file_dict[objname]['PP']
        self.ext_files = file_dict[objname]['EXT']
        self.tcorr_files = file_dict[objname]['TCORR']
        self.ccf_files = file_dict[objname]['CCF']
        self.lbl_rdb_files = file_dict[objname]['LBLRDB']
        # get the object_table row for this object
        self.object_table = object_table[object_table['_OBJECT'] == objname]
        # ---------------------------------------------------------------------
        # get spectrum output parameters (for page integration)
        self.spec_snr_plot_path = None
        self.spec_stats_table = None
        self.spec_download_table = None
        # get lbl output parameters (for page integration)
        self.lbl_combinations = []
        self.lbl_plot_path = dict()
        self.lbl_stats_table = dict()
        self.lbl_dwn_table = dict()
        # get ccf output parameters (for page integration)
        self.ccf_rv_plot_path = None
        self.ccf_med_plot_path = None
        self.ccf_stats_table = None
        self.ccf_download_table = None
        # get time series output parameters (for page integration)
        self.time_series_stats_table = None
        self.time_series_download_table = None
        # ---------------------------------------------------------------------
        # misc (save for use throughout)
        self.settings = settings
        self.profile = profile
        self.headers = headers
        # ---------------------------------------------------------------------
        # store the required header info
        self.header_dict = dict()
        # file lists to match COUNT COLS
        self.file_lists = [self.raw_files, self.pp_files, self.ext_files,
                           self.tcorr_files, self.ccf_files, None, None,
                           None, None, None]

    def populate_header_dict(self):
        """
        Populate the header dictionary with the required header keys
        :return:
        """
        # loop around COUNT COLS and populate header dict
        for it, file_kind in enumerate(COUNT_COLS):
            if file_kind in self.headers and self.file_lists[it] is not None:
                self.get_header_keys(self.headers[file_kind],
                                     self.file_lists[it])

    def get_header_keys(self, keys: Dict[str, Dict[str]], files: List[str]):
        """
        Get the header keys from the files
        :param keys: dictionary of keys to get (with their properties)
        :param files: list of files (for error reporting)
        :return:
        """

        # deal with no keys
        if len(keys) == 0:
            return
        # get an array of length of the files for each key
        for keydict in keys:
            # get key
            key = keys[keydict]['key']
            # deal with no files
            if len(files) == 0:
                self.header_dict[key] = None
            else:
                self.header_dict[key] = np.full(len(files), np.nan)
        # loop around files and populate the header_dict
        for pos, filename in enumerate(files):
            header = fits.getheader(filename)
            for keydict in keys:
                # get key (for header dict)
                key = keys[keydict]['key']
                # get value (for header dict)
                value = _header_value(keys[keydict], header, filename)
                # set value in header dict
                self.header_dict[key][pos] = value

    def get_spec_parameters(self):
        pass

    def get_lbl_parameters(self):
        """
        Get the LBL parameters for this object

        :return:
        """
        # don't go here is lbl rdb files are not present
        if len(self.lbl_rdb_files) == 0:
            return
        # ---------------------------------------------------------------------
        # generate place to save figures
        item_save_path = self.settings['ITEMS']
        item_rel_path = f'../{_ITEM_DIR}/'
        down_save_path = self.settings['DOWNS']
        down_rel_path = f'../{_DOWN_DIR}/'
        # ---------------------------------------------------------------------
        # get the ext h-band key
        ext_h_key = self.headers['LBL']['EXT_H']['key']
        # store the object+ template combinations
        lbl_objtmps = dict()
        # get the lbl objname+templates combinations
        for lbl_rdb_file in self.lbl_rdb_files:
            # get the basename
            basename = os.path.basename(lbl_rdb_file)
            # get the objname+template
            lbl_objtmp = basename.split(f'{LBL_SUFFIX}_')[-1].split('.rdb')[0]
            # append to list
            lbl_objtmps[lbl_objtmp] = lbl_rdb_file
        # loop around the objname+template combinations
        for lbl_objtmp in lbl_objtmps:
            # load rdb file
            rdb_table = Table.read(lbl_objtmps[lbl_objtmp], format='ascii.rdb')
            # get the values required
            lbl_rjd = np.array(rdb_table['rjd'])
            lbl_vrad = np.array(rdb_table['vrad'])
            lbl_svrad = np.array(rdb_table['svrad'])
            lbl_plot_date = np.array(rdb_table['plot_date'])
            lbl_snr_h = np.array(rdb_table[ext_h_key])
            # -----------------------------------------------------------------
            # plot the figure
            # -----------------------------------------------------------------
            # get the plot base name
            plot_base_name = 'lbl_plot_' + lbl_objtmp + '.png'
            # get the plot path
            plot_path = os.path.join(item_save_path, plot_base_name)
            # plot the lbl figure
            lbl_props = lbl_plot(lbl_plot_date, lbl_vrad, lbl_svrad, lbl_snr_h,
                                 plot_path, plot_title=lbl_objtmp)
            # -----------------------------------------------------------------
            # construct the stats
            # -----------------------------------------------------------------
            # get the stats base name
            stat_base_name = 'lbl_stat_' + lbl_objtmp + '.txt'
            # get the stat path
            stat_path = os.path.join(item_save_path, stat_base_name)
            # compute the stats
            lbl_stats_table(lbl_rjd, lbl_vrad, lbl_svrad, lbl_props,
                            stat_path)
            # -----------------------------------------------------------------
            # construct the download table
            # -----------------------------------------------------------------
            # get the download base name
            dwn_base_name = 'lbl_download_' + lbl_objtmp + '.txt'
            # get the download path
            download_path = os.path.join(down_save_path, dwn_base_name)
            # define the download files
            down_files = [lbl_objtmps[lbl_objtmp]]
            # define the download descriptions
            down_descs = ['RDB file']
            # compute the download table
            download_table(down_files, down_descs, down_rel_path,
                           download_path)
            # -----------------------------------------------------------------
            # update the paths
            self.lbl_plot_path[lbl_objtmp] = item_rel_path + plot_base_name
            self.lbl_stats_table[lbl_objtmp] = item_rel_path + stat_base_name
            self.lbl_dwn_table[lbl_objtmp] = item_rel_path + dwn_base_name
        # ---------------------------------------------------------------------
        # set the lbl combinations
        self.lbl_combinations = lbl_objtmps.keys()

    def get_ccf_parameters(self):
        pass

    def get_time_series_parameters(self):
        pass


def add_obj_pages(settings: dict, profile: dict, headers: dict,
                  object_table: Table, file_dict: FileDictReturn):
    # must import here (so that os.environ is set)
    # noinspection PyPep8Naming
    from apero.base.base import TQDM as tqdm
    from apero.core import constants
    from apero.core.core import drs_log
    from apero.tools.module.documentation import drs_markdown
    # get the parameter dictionary of constants from apero
    params = constants.load()
    # get WLOG
    wlog = drs_log.wlog
    # ------------------------------------------------------------------
    # get the parameters for lbl website url
    instrument = settings['INSTRUMENT']
    outdir = os.path.basename(settings['OBJ_OUT'])
    # ------------------------------------------------------------------
    # print progress
    wlog(params, '', 'Creating object pages')
    # change the object column to a url
    for it, row in tqdm(enumerate(object_table)):
        # get the object name for this row
        objname = row['_OBJECT']
        # create the object class
        object_instance = ObjectData(settings, profile, headers,
                                     objname, file_dict, object_table)
        # ---------------------------------------------------------------------
        # populate the header dictionary for this object instance
        wlog(params, '', f'\tPopulating header dictionary for {objname}')
        object_instance.populate_header_dict()
        # ---------------------------------------------------------------------
        # generate url for object
        object_url = f'{instrument}_{outdir}_{objname}_index'
        # replace  object name with the object name + object url
        object_table[OBJECT_COLUMN][it] = _make_url(objname, object_url)
        # create ARI object page
        object_page = drs_markdown.MarkDownPage(object_url)
        # get the profile name
        name = profile['profile name']
        # add title
        object_page.add_title(f'{objname} ({name})')
        # ---------------------------------------------------------------------
        # Add basic text
        # construct text to add
        object_page.add_text(f'This page was last modified: {Time.now()}')
        object_page.add_newline()
        # ---------------------------------------------------------------------
        # table of contents
        # ---------------------------------------------------------------------
        # Add the names of the sections
        names = ['Spectrum', 'LBL', 'CCF', 'Time series']
        # add the links to the pages
        items = [f'spectrum_objpage_{objname}',
                 f'lbl_objpage_{objname}',
                 f'ccf_objpage_{objname}', f'timeseries_objpage_{objname}']
        # add table of contents
        object_page.add_table_of_contents(items=items, names=names)
        # ---------------------------------------------------------------------
        # Spectrum section
        # ---------------------------------------------------------------------
        # print progress
        wlog(params, '', f'\tCreating spectrum section for {objname}')
        # add spectrum section
        objpage_spectrum(object_page, names[0], items[0], object_instance)
        # ---------------------------------------------------------------------
        # LBL section
        # ---------------------------------------------------------------------
        # print progress
        wlog(params, '', f'\tCreating LBL section for {objname}')
        # add LBL section
        objpage_lbl(object_page, names[1], items[1], object_instance)
        # ---------------------------------------------------------------------
        # CCF section
        # ---------------------------------------------------------------------
        # print progress
        wlog(params, '', f'\tCreating CCF section for {objname}')
        # add CCF section
        objpage_ccf(object_page, names[2], items[2], object_instance)
        # ---------------------------------------------------------------------
        # Time series section
        # ---------------------------------------------------------------------
        # print progress
        wlog(params, '', f'\tCreating time series section for {objname}')
        # add time series section
        objpage_timeseries(object_page, names[3], items[3], object_instance)
        # ---------------------------------------------------------------------
        # construct a path for the object name
        object_page_path = settings['OBJ_OUT']
        # construct the rst filename
        rst_filename = f'{objname}.rst'
        # save index page
        object_page.write_page(os.path.join(object_page_path, rst_filename))


def objpage_spectrum(page: Any, name: str, ref: str,
                     object_instance: ObjectData):
    # add divider
    page.add_divider(color=DIVIDER_COLOR, height=DIVIDER_HEIGHT)
    # add a reference to this section
    page.add_reference(ref)
    # add the section heading
    page.add_section(name)
    # ------------------------------------------------------------------
    # deal with no spectrum found
    if len(object_instance.ext_files) == 0:
        page.add_text('No spectrum found')
        return
    # ------------------------------------------------------------------
    # get the spectrum parameters
    object_instance.get_spec_parameters()
    # ------------------------------------------------------------------
    # add the snr plot
    if object_instance.spec_snr_plot_path is not None:
        # add the snr plot to the page
        page.add_image(object_instance.spec_snr_plot_path)
        # add a new line
        page.add_newline()
    # ------------------------------------------------------------------
    # add stats
    if object_instance.spec_stats_table is not None:
        # add some stats
        page.add_sub_section('Spectrum Stats')
        # add the stats table
        page.add_csv_table('', object_instance.spec_stats_table)
    # ------------------------------------------------------------------
    # add download links
    if object_instance.spec_download_table is not None:
        page.add_sub_section('Spectrum Downloads')
        # add the stats table
        page.add_csv_table('', object_instance.spec_download_table)


def objpage_lbl(page: Any, name: str, ref: str,
                object_instance: ObjectData):
    # add divider
    page.add_divider(color=DIVIDER_COLOR, height=DIVIDER_HEIGHT)
    # add a reference to this section
    page.add_reference(ref)
    # add the section heading
    page.add_section(name)
    # ------------------------------------------------------------------
    # deal with no spectrum found
    if len(object_instance.lbl_combinations) == 0:
        page.add_text('No LBL reduction found')
        return
    # ------------------------------------------------------------------
    # get the spectrum parameters
    object_instance.get_lbl_parameters()
    # ------------------------------------------------------------------
    # loop around the object+template combinations
    for objcomb in object_instance.lbl_combinations:
        # add the lbl plot
        if object_instance.lbl_plot_path[objcomb] is not None:
            # add the snr plot to the page
            page.add_image(object_instance.lbl_plot_path[objcomb])
            # add a new line
            page.add_newline()
        # ------------------------------------------------------------------
        # add stats
        if object_instance.lbl_stats_table[objcomb] is not None:
            # add some stats
            page.add_sub_section('LBL Stats')
            # add the stats table
            page.add_csv_table('', object_instance.lbl_stats_table[objcomb])
        # ------------------------------------------------------------------
        # add download links
        if object_instance.spec_download_table is not None:
            page.add_sub_section('LBL Downloads')
            # add the stats table
            page.add_csv_table('', object_instance.lbl_dwn_table[objcomb])


def objpage_ccf(page: Any, name: str, ref: str, object_instance: ObjectData):
    # add divider
    page.add_divider(color=DIVIDER_COLOR, height=DIVIDER_HEIGHT)
    # add a reference to this section
    page.add_reference(ref)
    # add the section heading
    page.add_section(name)
    # ------------------------------------------------------------------
    # get the spectrum parameters
    object_instance.get_ccf_parameters()
    # ------------------------------------------------------------------
    # add the ccf rv plot
    if object_instance.ccf_rv_plot_path is not None:
        # add the snr plot to the page
        page.add_image(object_instance.ccf_rv_plot_path)
        # add a new line
        page.add_newline()
    # add the ccf med plot
    if object_instance.ccf_rv_plot_path is not None:
        # add the snr plot to the page
        page.add_image(object_instance.ccf_med_plot_path)
        # add a new line
        page.add_newline()
    # ------------------------------------------------------------------
    # add stats
    if object_instance.ccf_stats_table is not None:
        # add some stats
        page.add_sub_section('CCF Stats')
        # add the stats table
        page.add_csv_table('', object_instance.ccf_stats_table)
    # ------------------------------------------------------------------
    # add download links
    if object_instance.ccf_download_table is not None:
        page.add_sub_section('CCF Downloads')
        # add the stats table
        page.add_csv_table('', object_instance.ccf_download_table)


def objpage_timeseries(page: Any, name: str, ref: str,
                       object_instance: ObjectData):
    # add divider
    page.add_divider(color=DIVIDER_COLOR, height=DIVIDER_HEIGHT)
    # add a reference to this section
    page.add_reference(ref)
    # add the section heading
    page.add_section(name)
    # ------------------------------------------------------------------
    # get the spectrum parameters
    object_instance.get_time_series_parameters()
    # ------------------------------------------------------------------
    # add stats
    if object_instance.time_series_stats_table is not None:
        # add some stats
        page.add_sub_section('Time Series Stats')
        # add the stats table
        page.add_csv_table('', object_instance.time_series_stats_table)
    # ------------------------------------------------------------------
    # add download links
    if object_instance.ccf_download_table is not None:
        page.add_sub_section('Time Series Downloads')
        # add the stats table
        page.add_csv_table('', object_instance.time_series_download_table)


# =============================================================================
# Plots, stat tables and download tables
# =============================================================================
def download_table(files: List[str], descriptions: List[str],
                   down_rel_path: str, down_path: str):
    """
    Generic download table saving to the item relative path

    :param files:
    :param descriptions:
    :param down_rel_path:
    :param down_path:
    :return:
    """
    # --------------------------------------------------------------------------
    # storage for ref file
    ref_paths = dict()
    # storage for outpath
    out_paths = dict()
    # storage for inpath
    in_paths = dict()
    # storage for the descriptions
    descs = dict()
    # get the out dir (item directory)
    out_dir = os.path.dirname(down_path)

    # loop around files and get the paths
    for it, filename in enumerate(files):
        # get the basename
        basename = os.path.basename(filename)
        # get the in path
        in_paths[basename] = filename
        # get the reference path (for the link)
        ref_paths[basename] = down_rel_path + basename
        # get the outpath
        out_paths[basename] = os.path.join(out_dir, basename)
        # get the descriptions
        descs[basename] = descriptions[it]
    # --------------------------------------------------------------------------
    # construct the stats table
    # --------------------------------------------------------------------------
    # start with a download dictionary
    down_dict = dict(Description=[], Value=[])
    # loop around files
    for basename in in_paths:
        # add the rdb file
        down_dict['Description'].append(descs[basename])
        down_dict['Value'].append(ref_paths[basename])
        # copy the file from in path to out path
        shutil.copy(in_paths[basename], out_paths[basename])
    # --------------------------------------------------------------------------
    # convert to table
    down_table = Table(down_dict)
    # write to file as csv file
    down_table.write(down_path, format='ascii.csv', overwrite=True)


def spec_plot():
    pass


def spec_stats_table():
    pass


def lbl_plot(lbl_plot_date: np.ndarray, lbl_vrad: np.ndarray,
             lbl_svrad: np.ndarray, lbl_snr_h: np.ndarray,
             plot_path: str, plot_title: str) -> Dict[str, Any]:
    # setup the figure
    fig, frame = plt.subplots(2, 1, figsize=(12, 6), sharex='all')
    # plot the points
    frame[0].plot_date(lbl_plot_date, lbl_vrad, marker='.', alpha=0.5,
                       color='green', ls='None')
    # plot the error bars
    frame[0].errorbar(lbl_plot_date, lbl_vrad, yerr=lbl_svrad,
                      marker='o', alpha=0.5, color='green', ls='None')
    # find percentile cuts that will be expanded by 150% for the ylim
    pp = np.nanpercentile(lbl_vrad, [10, 90])
    diff = pp[1] - pp[0]
    central_val = np.nanmean(pp)
    # used for plotting but also for the flagging of outliers
    ylim = [central_val - 1.5 * diff, central_val + 1.5 * diff]
    # length of the arrow flagging outliers
    l_arrow = (ylim[1] - ylim[0]) / 10.0
    # flag the low outliers
    low = lbl_vrad < ylim[0]
    # get the x and y values of the outliers to be looped over within
    # the arrow plotting
    xpoints = np.array(lbl_plot_date[low], dtype=float)
    ypoints = np.array(lbl_vrad[low], dtype=float)
    for ix in range(len(xpoints)):
        frame[0].arrow(xpoints[ix], ylim[0] + l_arrow * 2, 0, -l_arrow,
                       color='red',  head_width=0.5, head_length=1.5, alpha=0.5)
    # same as above for the high outliers
    high = lbl_vrad > ylim[1]
    xpoints = np.array(lbl_plot_date[high], dtype=float)
    ypoints = np.array(lbl_vrad[high], dtype=float)
    for ix in range(len(xpoints)):
        frame[0].arrow(xpoints[ix], ylim[1] - l_arrow * 2, 0, l_arrow,
                       color='red', head_width=0.5, head_length=1.5, alpha=0.5)
    # setting the plot
    frame[0].set(ylim=ylim)
    frame[0].set(title=plot_title)
    frame[0].grid(which='both', color='lightgray', linestyle='--')
    frame[0].set(ylabel='Velocity [m/s]')
    # simple plot of the SNR in a sample order. You need to
    # update the relevant ketword for SPIRou
    frame[1].plot_date(lbl_plot_date, lbl_snr_h, marker='.',
                       alpha=0.5, color='green', ls='None')
    frame[1].grid(which='both', color='lightgray', linestyle='--')
    frame[1].set(xlabel='Date')
    frame[1].set(ylabel='EXTSN060')
    plt.tight_layout()
    # save figure and close the plot
    plt.savefig(plot_path)
    plt.close()
    # some parameters are required later save them in a dictionary
    props = dict()
    props['low'] = low
    props['high'] = high
    props['ylim'] = ylim
    # return the props
    return props


def lbl_stats_table(lbl_rjd: np.ndarray, lbl_vrad: np.ndarray,
                    lbl_svrad: np.ndarray,
                    props: Dict[str, Any], stat_path: str):
    # --------------------------------------------------------------------------
    # compute the stats
    # --------------------------------------------------------------------------
    # get the 25, 50 and 75 percentile of the velocity uncertainty
    p_sigma = np.nanpercentile(lbl_svrad, [25,50,75])
    # get the 25, 50 and 75 percentile of the velocity
    v_sigma = np.nanpercentile(lbl_vrad - np.nanmedian(lbl_vrad), [25, 50, 75])
    # get the low outliers
    low = props['low']
    # get the high outliers
    high = props['high']
    # calculate the number of nights
    n_nights = len(np.unique(np.floor(lbl_rjd)))
    # calculate the systemetic velocity
    sys_vel = np.nanmedian(lbl_vrad)
    # calculate the velocity domain considered valid
    vel_domain = props['ylim']
    # --------------------------------------------------------------------------
    # construct the stats table
    # --------------------------------------------------------------------------
    # start with a stats dictionary
    stat_dict = dict(Description=[], Value=[])
    # add rv uncertainty
    stat_dict['Description'].append('RV Uncertainty (25, 50, 75 percentile)')
    stat_dict['Value'].append('{:.2f}, {:.2f}, {:.2f} m/s'.format(*p_sigma))
    # add the absolute deviation
    stat_dict['Description'].append('RV Absolute Deviation (25, 50, 75 '
                                    'percentile)')
    stat_dict['Value'].append('{:.2f}, {:.2f}, {:.2f} m/s'.format(*v_sigma))
    # add the number of measurements
    stat_dict['Description'].append('Number of Measurements')
    stat_dict['Value'].append(len(lbl_vrad))
    # add the spurious low points
    stat_dict['Description'].append('Number of Spurious Low Points')
    stat_dict['Value'].append(np.sum(low))
    # add the spurious high points
    stat_dict['Description'].append('Number of Spurious High Points')
    stat_dict['Value'].append(np.sum(high))
    # add the number of nights
    stat_dict['Description'].append('Number of Nights')
    stat_dict['Value'].append(n_nights)
    # add the systemic velocity
    stat_dict['Description'].append('Systemic Velocity')
    stat_dict['Value'].append('{:.2f} m/s'.format(sys_vel))
    # add the Velocity domain considered
    stat_dict['Description'].append('Velocity Domain considered valid')
    stat_dict['Value'].append('{:.2f} - {:.2f} m/s'.format(*vel_domain))
    # --------------------------------------------------------------------------
    # convert to table
    stat_table = Table(stat_dict)
    # write to file as csv file
    stat_table.write(stat_path, format='ascii.csv', overwrite=True)


def ccf_plot():
    pass


def ccf_stats_table():
    pass


def time_series_stats_table():
    pass

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


def read_yaml_file(profile_filename: str) -> dict:
    """
    Read the yaml file  for the list of profiles using pyyaml

    :param profile_filename: str, the filename of the yaml file

    :return:
    """
    # read the yaml file "profiles.yaml" for the list of profiles using pyyaml
    with open(profile_filename, 'r') as stream:
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
    # any cwidths that are 0 are set to 1
    cwidths[cwidths == 0] = 1
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


def _make_url(value: str, url: str) -> str:
    """
    Make a url from a value and a url
    """
    return f':ref:`{value} <{url}>`'


def _header_value(keydict: Dict[str, str], header: fits.Header,
                  filename: str):
    # get properties of header key
    key = keydict['key']
    unit = keydict.get('unit', None)
    dtype = keydict.get('dtype', None)
    timefmt = keydict.get('timefmt', None)
    # -------------------------------------------------------------------------
    # get raw value from header
    rawvalue = header.get(key, None)
    # deal with no value
    if rawvalue is None:
        raise ValueError(f'HeaderKey: {key} not found in header'
                         f'\n\tFile: {filename}')
    # -------------------------------------------------------------------------
    # deal with dtype
    if dtype is not None:
        try:
            if dtype == 'int':
                rawvalue = int(rawvalue)
            elif dtype == 'float':
                rawvalue = float(rawvalue)
            elif dtype == 'bool':
                rawvalue = bool(rawvalue)
        except Exception as _:
            raise ValueError(f'HeaderDtype: {dtype} not valid for '
                             f'{key}={rawvalue}\n\tFile: {filename}')
    # -------------------------------------------------------------------------
    # deal with time
    if timefmt is not None:
        try:
            return Time(rawvalue, format=timefmt)
        except Exception as _:
            raise ValueError(f'HeaderTime: {timefmt} not valid for '
                             f'{key}={rawvalue}\n\tFile: {filename}')
    # -------------------------------------------------------------------------
    # deal with units
    if unit is not None:
        try:
            rawvalue = uu.Quantity(rawvalue, unit)
        except ValueError:
            raise ValueError(f'HeaderUnit: {unit} not valid for '
                             f'{key}={rawvalue}\n\tFile: {filename}')
    # -------------------------------------------------------------------------
    # return the raw value
    return rawvalue


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # assume the first argument is the profile filename
    if len(sys.argv) != 2:
        raise IOError('Please supply the profile filename')
    profile_filename = sys.argv[1]
    # step 1: read the yaml file
    apero_profiles = read_yaml_file(profile_filename)
    # get global settings
    gsettings = dict(apero_profiles['settings'])
    header_settings = dict(apero_profiles['headers'])
    del apero_profiles['settings']
    del apero_profiles['headers']
    # deal with a reset
    if gsettings['reset']:
        # remove working directory
        shutil.rmtree(gsettings['working directory'], ignore_errors=True)
    # ----------------------------------------------------------------------
    # step 2: for each profile compile all stats
    all_apero_stats = dict()
    # loop around profiles from yaml file
    for apero_profile_name in apero_profiles:
        # sort out settings
        ari_settings = get_settings(gsettings, apero_profile_name)
        # add profile name to settings
        apero_profiles[apero_profile_name]['profile name'] = apero_profile_name
        # we reprocess if the file does not exist or if REPROCESS is True
        if REPROCESS:
            # print progress
            print('=' * 50)
            print('Compiling stats for profile: {0}'.format(apero_profile_name))
            print('=' * 50)
            # get profile
            apero_profile = apero_profiles[apero_profile_name]
            # compile stats
            apero_stats = compile_stats(ari_settings, apero_profile,
                                        header_settings)
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
    ari_settings = get_settings(gsettings)
    # ----------------------------------------------------------------------
    # step 3: write markdown files
    write_markdown(gsettings, ari_settings, all_apero_stats)
    # ----------------------------------------------------------------------
    # step 4: compile sphinx files
    compile_docs(ari_settings)
    # ----------------------------------------------------------------------
    # step 5: upload to hosting
    upload_docs(gsettings, ari_settings)

# =============================================================================
