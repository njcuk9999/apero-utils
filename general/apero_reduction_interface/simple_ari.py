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
from scipy.optimize import curve_fit

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
# list tables to load
TABLE_NAMES = ['OBJECT_TABLE', 'RECIPE_TABLE', 'MESSAGE_TABLE']
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
COUNT_COLS = ['RAW FILE  ',
              'PP FILE   ',
              'EXT FILE  ',
              'TCORR FILE',
              'CCF FILE  ',
              'POL FILE  ',
              'e.fits    ',
              't.fits    ',
              'v.fits    ',
              'p.fits    ']
# header cols
HEADER_COL = ['RAW', 'PP', 'EXT', 'TCORR', 'CCF', 'POL', 'e', 't', 'v', 'p']
# time series column names
TIME_SERIES_COLS = ['Obs Dir', 'First obs mid',
                    'Last obs mid', 'Number of obs', 'Seeing', 'Airmass',
                    'Mean Exptime', 'Total Exptime', 'DPRTYPEs']
# define chains (if this number is zero do not count)
COUNT_CHAINS = [None, COUNT_COLS[0], COUNT_COLS[1], COUNT_COLS[2],
                COUNT_COLS[3], COUNT_COLS[3], COUNT_COLS[1], COUNT_COLS[2],
                COUNT_COLS[3], COUNT_COLS[3]]
# define the lbl rdb suffix (lbl or lbl2)
LBL_SUFFIX = 'lbl'
# define the LBL stat dir (inside the lbl directory)
LBL_STAT_DIR = 'lblstats'
# define the LBL stat files {0} is for the object name + template name
#  i.e. {objname}_{template}
LBL_STAT_FILES = dict()
LBL_STAT_FILES['LBL Diagnostic Plots'] = 'lbl_{0}_plots.pdf'
# define how many ccf files to use
MAX_NUM_CCF = 100
# object page styling
DIVIDER_COLOR = '#FFA500'
DIVIDER_HEIGHT = 6
PLOT_BACKGROUND_COLOR = '#FEFDE1'
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
    if profile_name is not None and not os.path.exists(param_settings['OBJ_OUT']):
        os.makedirs(param_settings['OBJ_OUT'])
    if profile_name is not None and not os.path.exists(param_settings['ITEMS']):
        os.makedirs(param_settings['ITEMS'])
    if profile_name is not None and not os.path.exists(param_settings['DOWNS']):
        os.makedirs(param_settings['DOWNS'])
    # return all parameter settings
    return param_settings


def compile_stats(gsettings: dict, settings: dict, profile: dict,
                  headers: dict) -> dict:
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
        profile_stats[TABLE_NAMES[0]] = Table.read(object_table_file)
    elif skip_obj_table:
        profile_stats[TABLE_NAMES[0]] = None
    else:
        # get the object table (astropy table)
        object_table, filedict = compile_apero_object_table(gsettings)
        # add the lbl count
        object_table, filedict = add_lbl_count(profile, object_table, filedict)
        # add the object pages
        object_table = add_obj_pages(gsettings, settings, profile, headers,
                                     object_table, filedict)
        # add final object table to profile stats
        profile_stats[TABLE_NAMES[0]] = object_table
    # ------------------------------------------------------------------
    # deal with skipping recipe table
    if skip_recipe_table and os.path.exists(recipe_table_file):
        profile_stats[TABLE_NAMES[1]] = Table.read(recipe_table_file)
    elif skip_recipe_table:
        profile_stats[TABLE_NAMES[1]] = None
    else:
        # get the recipe log table
        profile_stats[TABLE_NAMES[1]] = compile_apero_recipe_table()
    # ------------------------------------------------------------------
    # deal with skipping message table
    if skip_msg_table and os.path.exists(message_table_file):
        profile_stats[TABLE_NAMES[2]] = Table.read(message_table_file)
    elif skip_msg_table:
        profile_stats[TABLE_NAMES[2]] = None
    else:
        # get the message log table
        profile_stats[TABLE_NAMES[2]] = compile_apero_message_table()
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
    index_page.add_text('Please note: Your object may be under another '
                        'name. Please check `here <https://docs.google.com/'
                        'spreadsheets/d/'
                        '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/'
                        'edit?usp=sharing>`_, the name displayed in ARI will '
                        'be the first column [OBJNAME]')
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
        settings['TABLE_REFS'] = dict()
        # create a page
        profile_page = drs_markdown.MarkDownPage(cprofile_name)
        # add title
        profile_page.add_title(profile_name)
        # store the reference name for profile page table of contents
        table_files = []
        # loop around tables
        for table_name in TABLE_NAMES:
            # get table
            table = stats[profile_name][table_name]
            # create a table page for this table
            table_ref_name = cprofile_name.lower() + '_' + table_name.lower()
            # make a markdown page for the table
            table_page = drs_markdown.MarkDownPage(table_ref_name)
            # add object table
            table_filename = f'{table_name}.csv'
            table_title = table_name.lower().replace('_', ' ')
            title = f'APERO reduction {table_title} ({cprofile_name})'
            table_page.add_title(title)
            # deal with column widths for this file type
            if table is not None:
                # add the csv version of this table
                table_page.add_csv_table('', f'../{_DATA_DIR}/' +
                                         table_filename, cssclass='csvtable2')
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
FileDictReturn = Dict[str, Dict[str, Any]]


def compile_apero_object_table(gsettings) -> Tuple[Table, FileDictReturn]:
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
    # deal with filtering by object
    if gsettings['filter objects']:
        subconditions = []
        for objname in gsettings['objects']:
            subconditions.append(f'OBJNAME="{objname}"')
        condition = ' OR '.join(subconditions)
        condition = f'({condition})'
    else:
        condition = None
    # ------------------------------------------------------------------
    # log that we are loading
    # get the object table from the astrometric database
    object_table = astrodbm.get_entries(columns=','.join(ASTROMETRIC_COLUMNS),
                                        condition=condition)
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
    # add counting columns to the object table
    object_table['DPRTYPES'] = [' '*255] * len(object_table)
    object_table[COUNT_COLS[0]] = [0] * len(object_table)
    object_table['FIRST_RAW'] = [None] * len(object_table)
    object_table['LAST_RAW'] = [None] * len(object_table)
    object_table[COUNT_COLS[1]] = [0] * len(object_table)
    object_table[COUNT_COLS[2]] = [0] * len(object_table)
    object_table[COUNT_COLS[3]] = [0] * len(object_table)
    object_table[COUNT_COLS[4]] = [0] * len(object_table)
    # deal with instruments that have polarimetry
    if has_polar:
        object_table[COUNT_COLS[5]] = [0] * len(object_table)
    object_table[COUNT_COLS[6]] = [0] * len(object_table)
    object_table[COUNT_COLS[7]] = [0] * len(object_table)
    object_table[COUNT_COLS[8]] = [0] * len(object_table)
    # deal with instruments that have polarimetry
    if has_polar:
        object_table[COUNT_COLS[9]] = [0] * len(object_table)
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
        objname = object_table[OBJECT_COLUMN][pos]
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
        # Add a dpr type column
        dprtypes = indexdbm.get_entries('KW_DPRTYPE', condition=ext_cond)
        object_table['DPRTYPES'][pos] = ','.join(list(np.unique(dprtypes)))
        # ------------------------------------------------------------------
        # deal with finding the first and last raw file
        # ------------------------------------------------------------------
        times = indexdbm.get_entries('KW_MID_OBS_TIME', condition=raw_cond)
        # if there are no entries we have no raw files for this object
        if len(times) == 0:
            continue
        # find minimum time value and convert to human time
        first_time = Time(np.min(times), format='mjd')
        object_table['FIRST_RAW'][pos] = first_time
        # find maximum time value and convert to human time
        last_time = Time(np.max(times), format='mjd')
        object_table['LAST_RAW'][pos] = last_time
        # ------------------------------------------------------------------
        # deal with getting the observation directories (For EXT files)
        # ------------------------------------------------------------------
        obs_dirs = indexdbm.get_entries('OBS_DIR', condition=ext_cond)
        # get the observation directories
        file_dict[objname]['OBS_DIRS'] = obs_dirs
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
    mask = object_table[COUNT_COLS[0]] > 0
    # apply mask
    object_table = object_table[mask]
    # ------------------------------------------------------------------
    # remove objects for file_dict with zero raw entries
    for objname in object_table[OBJECT_COLUMN]:
        # deal with objects not in the file dict
        if objname not in file_dict:
            continue
        # deal with objects with no raw files
        if len(file_dict[objname][COUNT_COLS[0]]) == 0:
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
# Functions for compiling the lbl/ccf stats
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
    objnames = np.array(object_table[OBJECT_COLUMN])
    # loop around objects
    for pos, objname in tqdm(enumerate(objnames)):
        # ---------------------------------------------------------------------
        # LBL RV files
        # ---------------------------------------------------------------------
        # for each object find all directories in lbl path that match this
        #   object name
        lblrv_dir = glob.glob(os.path.join(lbl_path, 'lblrv', f'{objname}_*'))
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
        # ---------------------------------------------------------------------
        # LBL RDB files
        # ---------------------------------------------------------------------
        # get all the lbl rdb files for this object name
        rdb_glob = f'{LBL_SUFFIX}_{objname}_*.rdb'
        all_lblrdb_files = glob.glob(os.path.join(lbl_path, 'lblrdb', rdb_glob))
        # remove drift files from the lbl rdb files
        lblrdb_files = []
        for lblrdb_file in all_lblrdb_files:
            if 'drift' not in lblrdb_file:
                lblrdb_files.append(lblrdb_file)
        # add list to the LBLRDB file dict for this object
        file_dict[objname]['LBLRDB'] = lblrdb_files
        # ---------------------------------------------------------------------
        # LBL Stats (generated by Charles)
        # ---------------------------------------------------------------------
        # define lbl stat dir
        lbl_stat_dir = os.path.join(lbl_path, LBL_STAT_DIR)
        # loop around the stat files
        for lblfilekey in LBL_STAT_FILES:
            # get obj+template glob
            objtmp_glob = f'{objname}_*'
            # get all obj+template directories
            lblsdirs = glob.glob(os.path.join(lbl_stat_dir, objtmp_glob))
            # make the file_dict
            file_dict[objname][lblfilekey] = dict()
            # loop around obj_templates
            for objtmpdir in lblsdirs:
                # get the obj+tmp name
                objtmp = os.path.basename(objtmpdir)
                # get the expected stats file name
                lblsfile = LBL_STAT_FILES[lblfilekey].format(objtmp)
                # get the expected stats file path
                lblspath = os.path.join(objtmpdir, lblsfile)
                # check that expected file exists
                if os.path.exists(lblspath):
                    # add a list to file dict for this object
                    file_dict[objname][lblfilekey][objtmp] = lblspath
    # -------------------------------------------------------------------------
    # add to object table
    object_table['LBL'] = lbl_count
    # -------------------------------------------------------------------------
    # return the object table
    return object_table, file_dict


def choose_ccf_files(ccf_props: Dict[str, Any]) -> List[str]:
    """
    Choose CCF files based on the most numerious of a single mask
    and then the MAX_NUM_CCF selected uniformly in time
    """
    from apero.core.utils import drs_utils
    # get parameters from props
    masks = ccf_props['masks']
    files = np.array(ccf_props['files'])
    mjd = ccf_props['mjd']
    # storage of the ccf mask count
    mask_names, mask_counts = [], []
    # loop around masks
    for mask in masks:
        # count how many files of each mask
        mask_names = np.append(mask_names, mask)
        mask_counts = np.append(mask_counts, np.sum(masks == mask))
    # choose the mask with the most counts
    chosen_mask = mask_names[np.argmax(mask_counts)]
    # filter files and mjd by chosen mask
    mask_mask = masks == chosen_mask
    sfiles = files[mask_mask]
    smjd = mjd[mask_mask]
    # now choose files distributed equally in time
    time_mask = drs_utils.uniform_time_list(smjd, MAX_NUM_CCF)
    # filter files by the time mask
    sfiles = sfiles[time_mask]
    # return ccf props
    return sfiles


def fit_ccf(ccf_props: Dict[str, Any]) -> Dict[str, Any]:
    """
    Fit the median CCF in order to plot the graph

    :param ccf_props:
    :return:
    """
    from apero.core.math.gauss import gauss_function
    # get parameters from props
    rv_vec = ccf_props['rv_vec']
    med_ccf = ccf_props['med_ccf']
    y1_1sig = ccf_props['y1_1sig']
    y2_1sig = ccf_props['y2_1sig']
    # set up guess
    amp0 = 1 - med_ccf[np.argmin(med_ccf)]
    pos0 = rv_vec[np.argmin(med_ccf)]
    sig0 = 4.0
    dc0 = 1.0
    guess = [-amp0, pos0, sig0, dc0]
    # try to fit the median ccf
    try:
        coeffs, _ = curve_fit(gauss_function, rv_vec, med_ccf, p0=guess)
        fit = gauss_function(rv_vec, *coeffs)
        xlim = [coeffs[1] - coeffs[2] * 20, coeffs[1] + coeffs[2] * 20]
        ylim = [np.min(y1_1sig - fit), np.max(y2_1sig - fit)]
        has_fit = True
    except Exception as _:
        fit = np.full(len(rv_vec), np.nan)
        xlim = [np.min(rv_vec), np.max(rv_vec)]
        ylim = [np.min(med_ccf), np.max(med_ccf)]
        has_fit = False
    # adjust ylim
    ylim = [ylim[0] - 0.1 * (ylim[1] - ylim[0]),
            ylim[1] + 0.1 * (ylim[1] - ylim[1])]
    # add to ccf_props
    ccf_props['has_fit'] = has_fit
    ccf_props['fit'] = fit
    ccf_props['xlim'] = xlim
    ccf_props['ylim'] = ylim
    # return ccf props
    return ccf_props


# =============================================================================
# Functions for making object pages
# =============================================================================
class ObjectData:
    def __init__(self, profile, settings, headers,
                 objname, file_dict, object_table):
        # get the object name
        self.objname = objname
        # get the file_dictionary for this object
        self.raw_files = file_dict[objname].get(COUNT_COLS[0], [])
        self.pp_files = file_dict[objname].get(COUNT_COLS[1], [])
        self.ext_files = file_dict[objname].get(COUNT_COLS[2], [])
        self.tcorr_files = file_dict[objname].get(COUNT_COLS[3], [])
        self.ccf_files = file_dict[objname].get(COUNT_COLS[4], [])
        self.lbl_rdb_files = file_dict[objname].get('LBLRDB', [])
        # add the lbl stat files into a dictionary
        self.lbl_stat_files = dict()
        # loop around lbl stats files and load lblfilekey dict in
        #   each dict should contain obj+temp combinations
        for lblfilekey in LBL_STAT_FILES:
            lblsfiles = file_dict[objname].get(lblfilekey, dict())
            self.lbl_stat_files[lblfilekey] = lblsfiles
        # mask for this object
        objmask = object_table[OBJECT_COLUMN] == objname
        # get the object_table row for this object
        self.object_table = object_table[objmask]
        # get the observation directories for this object
        self.obs_dirs = file_dict[objname].get('OBS_DIRS', [])
        # ---------------------------------------------------------------------
        # get spectrum output parameters (for page integration)
        self.spec_plot_path = None
        self.spec_stats_table = None
        self.spec_dwn_table = None
        # get lbl output parameters (for page integration)
        self.lbl_combinations = []
        self.lbl_plot_path = dict()
        self.lbl_stats_table = dict()
        self.lbl_dwn_table = dict()
        # get ccf output parameters (for page integration)
        self.ccf_plot_path = None
        self.ccf_stats_table = None
        self.ccf_dwn_table = None
        # get time series output parameters (for page integration)
        self.time_series_stats_table = None
        self.time_series_dwn_table = None
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
        for it, file_kind in enumerate(HEADER_COL):
            if file_kind in self.headers and self.file_lists[it] is not None:
                self.get_header_keys(self.headers[file_kind],
                                     self.file_lists[it])

    def get_header_keys(self, keys: Dict[str, Dict[str, str]], files: List[str]):
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
            kstore = keys[keydict]
            dtype = kstore.get('dtype', None)
            timefmt = kstore.get('timefmt', None)
            unit = kstore.get('unit', None)
            # deal with no files
            if len(files) == 0:
                self.header_dict[keydict] = None
            elif dtype == 'float' and timefmt is None and unit is None:
                self.header_dict[keydict] = np.full(len(files), np.nan)
            elif unit is not None:
                null_list = [np.nan] * len(files)
                self.header_dict[keydict] = uu.Quantity(null_list, unit)
            else:
                self.header_dict[keydict] = np.array([None] * len(files))
        # loop around files and populate the header_dict
        for pos, filename in enumerate(files):
            header = fits.getheader(filename)
            for keydict in keys:
                # get value (for header dict)
                value = _header_value(keys[keydict], header, filename)
                # set value in header dict
                self.header_dict[keydict][pos] = value

    def get_spec_parameters(self):
        # don't go here is lbl rdb files are not present
        if len(self.ext_files) == 0:
            return
        # ---------------------------------------------------------------------
        # generate place to save figures
        item_save_path = self.settings['ITEMS']
        item_rel_path = f'../{_ITEM_DIR}/'
        down_save_path = self.settings['DOWNS']
        down_rel_path = f'../{_DOWN_DIR}/'
        # -----------------------------------------------------------------
        # storage for ccf values
        spec_props = dict()
        # get values for use in plot
        spec_props['mjd'] = Time(np.array(self.header_dict['EXT_MJDMID']))
        spec_props['EXT_Y'] = np.array(self.header_dict['EXT_Y'])
        spec_props['EXT_H'] = np.array(self.header_dict['EXT_H'])
        spec_props['EXT_Y_LABEL'] = self.headers['EXT']['EXT_Y']['label']
        spec_props['EXT_H_LABEL'] = self.headers['EXT']['EXT_H']['label']
        spec_props['NUM_RAW_FILES'] = len(self.raw_files)
        spec_props['NUM_PP_FILES'] = len(self.pp_files)
        spec_props['NUM_EXT_FILES'] = len(self.ext_files)
        spec_props['NUM_TCORR_FILES'] = len(self.tcorr_files)
        spec_props['FIRST_RAW'] = Time(self.object_table['FIRST_RAW'][0])
        spec_props['LAST_RAW'] = Time(self.object_table['LAST_RAW'][0])
        spec_props['FIRST_PP'] = Time(np.min(self.header_dict['PP_MJDMID']))
        spec_props['LAST_PP'] = Time(np.max(self.header_dict['PP_MJDMID']))
        spec_props['FIRST_EXT'] = Time(np.min(self.header_dict['EXT_MJDMID']))
        spec_props['LAST_EXT'] = Time(np.max(self.header_dict['EXT_MJDMID']))
        spec_props['FIRST_TCORR'] = Time(np.min(self.header_dict['TCORR_MJDMID']))
        spec_props['LAST_TCORR'] = Time(np.max(self.header_dict['TCORR_MJDMID']))
        # -----------------------------------------------------------------
        # plot the figure
        # -----------------------------------------------------------------
        # get the plot base name
        plot_base_name = 'spec_plot_' + self.objname + '.png'
        # get the plot path
        plot_path = os.path.join(item_save_path, plot_base_name)
        # plot the lbl figure
        spec_plot(spec_props, plot_path, plot_title=f'{self.objname}')
        # -----------------------------------------------------------------
        # construct the stats
        # -----------------------------------------------------------------
        # get the stats base name
        stat_base_name = 'spec_stat_' + self.objname + '.txt'
        # get the stat path
        stat_path = os.path.join(item_save_path, stat_base_name)
        # compute the stats
        spec_stats_table(spec_props, stat_path, title='Spectrum stats')
        # -----------------------------------------------------------------
        # Create the file lists for this object
        # -----------------------------------------------------------------
        # construct the save path for ext files
        ext_file = os.path.join(down_save_path, 'ext_file_list.txt')
        create_file_list(self.ext_files, ext_file)
        # construct the save path for the tcorr files
        tcorr_file = os.path.join(down_save_path, 'tcorr_file_list.txt')
        create_file_list(self.tcorr_files, tcorr_file)
        # -----------------------------------------------------------------
        # construct the download table
        # -----------------------------------------------------------------
        # get the download base name
        dwn_base_name = 'spec_download_' + self.objname + '.txt'
        # get the download table path
        item_path = os.path.join(item_save_path, dwn_base_name)
        # define the download files
        down_files = [ext_file, tcorr_file]
        # define the download descriptions
        down_descs = ['Extracted 2D spectra', 'Telluric corrected 2D spectra']
        # compute the download table
        download_table(down_files, down_descs, item_path, down_rel_path,
                       down_save_path, title='Spectrum Downloads')
        # -----------------------------------------------------------------
        # update the paths
        self.spec_plot_path = item_rel_path + plot_base_name
        self.spec_stats_table = item_rel_path + stat_base_name
        self.spec_dwn_table = item_rel_path + dwn_base_name

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
        # storage of properties
        lbl_props = dict()
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
            lbl_props['rjd'] = np.array(rdb_table['rjd'])
            lbl_props['vrad'] = np.array(rdb_table['vrad'])
            lbl_props['svrad'] = np.array(rdb_table['svrad'])
            lbl_props['plot_date'] = np.array(rdb_table['plot_date'])
            lbl_props['snr_h'] = np.array(rdb_table[ext_h_key])
            lbl_props['SNR_H_LABEL'] = self.headers['LBL']['EXT_H']['label']
            # -----------------------------------------------------------------
            # plot the figure
            # -----------------------------------------------------------------
            # get the plot base name
            plot_base_name = 'lbl_plot_' + lbl_objtmp + '.png'
            # get the plot path
            plot_path = os.path.join(item_save_path, plot_base_name)
            # plot the lbl figure
            lbl_props = lbl_plot(lbl_props, plot_path,
                                 plot_title=f'LBL {lbl_objtmp}')
            # -----------------------------------------------------------------
            # construct the stats
            # -----------------------------------------------------------------
            # get the stats base name
            stat_base_name = 'lbl_stat_' + lbl_objtmp + '.txt'
            # get the stat path
            stat_path = os.path.join(item_save_path, stat_base_name)
            # compute the stats
            lbl_stats_table(lbl_props, stat_path, title='LBL stats')
            # -----------------------------------------------------------------
            # construct the download table
            # -----------------------------------------------------------------
            # get the download base name
            dwn_base_name = 'lbl_download_' + lbl_objtmp + '.txt'
            # get the download table path
            item_path = os.path.join(item_save_path, dwn_base_name)
            # define the download files
            down_files = [lbl_objtmps[lbl_objtmp]]
            # define the download descriptions
            down_descs = ['RDB file']
            # Add the lbl stat files
            for lblfilekey in LBL_STAT_FILES:
                # deal with no lbl file key for this object
                if lblfilekey not in self.lbl_stat_files:
                    continue
                # deal with this obj+templ being present
                if lbl_objtmp in self.lbl_stat_files[lblfilekey]:
                    # get lblsfile
                    lblsfile = self.lbl_stat_files[lblfilekey][lbl_objtmp]
                    # add to the list
                    down_files.append(lblsfile)
                    # add to the list
                    down_descs.append(lblfilekey)
            # compute the download table
            download_table(down_files, down_descs, item_path, down_rel_path,
                           down_save_path, title='LBL Downloads')
            # -----------------------------------------------------------------
            # update the paths
            self.lbl_plot_path[lbl_objtmp] = item_rel_path + plot_base_name
            self.lbl_stats_table[lbl_objtmp] = item_rel_path + stat_base_name
            self.lbl_dwn_table[lbl_objtmp] = item_rel_path + dwn_base_name
        # ---------------------------------------------------------------------
        # set the lbl combinations
        self.lbl_combinations = lbl_objtmps.keys()

    def get_ccf_parameters(self):

        from apero.core.math import normal_fraction
        # don't go here is lbl rdb files are not present
        if len(self.ccf_files) == 0:
            return
        # ---------------------------------------------------------------------
        # generate place to save figures
        item_save_path = self.settings['ITEMS']
        item_rel_path = f'../{_ITEM_DIR}/'
        down_save_path = self.settings['DOWNS']
        down_rel_path = f'../{_DOWN_DIR}/'
        # -----------------------------------------------------------------
        # storage for ccf values
        ccf_props = dict()
        # get values for use in plot
        ccf_props['mjd'] = Time(np.array(self.header_dict['CCF_MJDMID']))
        dv_vec = self.header_dict['CCF_DV'].to(uu.m / uu.s).value
        ccf_props['dv'] = np.array(dv_vec)
        sdv_vec = self.header_dict['CCF_SDV'].to(uu.m / uu.s).value
        ccf_props['sdv'] = np.array(sdv_vec)
        fwhm_vec = self.header_dict['CCF_FWHM'].to(uu.m / uu.s).value
        ccf_props['fwhm'] = np.array(fwhm_vec)
        ccf_props['masks'] = np.array(self.header_dict['CCF_MASK'])
        ccf_props['files'] = np.array(self.ccf_files)
        # -----------------------------------------------------------------
        # select ccf files to use
        select_files = choose_ccf_files(ccf_props)
        # load the first file to get the rv vector
        ccf_table0 = Table.read(select_files[0], format='fits', hdu=1)
        # get the rv vector
        ccf_props['rv_vec'] = ccf_table0['RV']
        # storage for the CCF vectors
        all_ccf = np.zeros((len(select_files), len(ccf_table0)))
        # loop around all other files, load them and load into all_ccf
        for row, select_file in enumerate(select_files):
            table_row = Table.read(select_file, format='fits', hdu=1)
            # get the combined CCF for this file
            ccf_row = table_row['Combined']
            # normalize ccf
            ccf_row = ccf_row / np.nanmedian(ccf_row)
            # push into vector
            all_ccf[row] = ccf_row
        # -----------------------------------------------------------------
        # get the 1 and 2 sigma limits
        lower_sig1 = 100 * (0.5 - normal_fraction(1) / 2)
        upper_sig1 = 100 * (0.5 + normal_fraction(1) / 2)
        lower_sig2 = 100 * (0.5 - normal_fraction(2) / 2)
        upper_sig2 = 100 * (0.5 + normal_fraction(2) / 2)
        # -----------------------------------------------------------------
        # y1 1sig is the 15th percentile of all ccfs
        ccf_props['y1_1sig'] = np.nanpercentile(all_ccf, lower_sig1, axis=0)
        # y2 1sig is the 84th percentile of all ccfs
        ccf_props['y2_1sig'] = np.nanpercentile(all_ccf, upper_sig1, axis=0)
        # y1 1sig is the 15th percentile of all ccfs
        ccf_props['y1_2sig'] = np.nanpercentile(all_ccf, lower_sig2, axis=0)
        # y2 1sig is the 84th percentile of all ccfs
        ccf_props['y2_2sig'] = np.nanpercentile(all_ccf, upper_sig2, axis=0)
        # med ccf is the median ccf (50th percentile)
        ccf_props['med_ccf'] = np.nanmedian(all_ccf, axis=0)
        # delete all_ccf to save memeory
        del all_ccf
        # fit the median ccf
        ccf_props = fit_ccf(ccf_props)
        # -----------------------------------------------------------------
        # plot the figure
        # -----------------------------------------------------------------
        # get the plot base name
        plot_base_name = 'ccf_plot_' + self.objname + '.png'
        # get the plot path
        plot_path = os.path.join(item_save_path, plot_base_name)
        # plot the lbl figure
        ccf_plot(ccf_props, plot_path, plot_title=f'CCF {self.objname}')
        # -----------------------------------------------------------------
        # construct the stats
        # -----------------------------------------------------------------
        # get the stats base name
        stat_base_name = 'ccf_stat_' + self.objname + '.txt'
        # get the stat path
        stat_path = os.path.join(item_save_path, stat_base_name)
        # compute the stats
        ccf_stats_table(ccf_props, stat_path, title='CCF stats')
        # -----------------------------------------------------------------
        # Create the file lists for this object
        # -----------------------------------------------------------------
        # construct the save path for ccf files
        ccf_file = os.path.join(down_save_path, 'ccf_file_list.txt')
        create_file_list(self.ccf_files, ccf_file)
        # -----------------------------------------------------------------
        # construct the download table
        # -----------------------------------------------------------------
        # get the download base name
        dwn_base_name = 'ccf_download_' + self.objname + '.txt'
        # get the download table path
        item_path = os.path.join(item_save_path, dwn_base_name)
        # define the download files
        down_files = [ccf_file]
        # define the download descriptions
        down_descs = ['CCF Table']
        # compute the download table
        download_table(down_files, down_descs, item_path, down_rel_path,
                       down_save_path, title='CCF Downloads')
        # -----------------------------------------------------------------
        # update the paths
        self.ccf_plot_path = item_rel_path + plot_base_name
        self.ccf_stats_table = item_rel_path + stat_base_name
        self.ccf_dwn_table = item_rel_path + dwn_base_name

    def get_time_series_parameters(self):
        # don't go here is lbl rdb files are not present
        if len(self.ext_files) == 0:
            return
        # ---------------------------------------------------------------------
        # generate place to save figures
        item_save_path = self.settings['ITEMS']
        item_rel_path = f'../{_ITEM_DIR}/'
        down_save_path = self.settings['DOWNS']
        down_rel_path = f'../{_DOWN_DIR}/'
        # -----------------------------------------------------------------
        # storage for ccf values
        time_series_props = dict()
        # get labels
        snr_y_label = self.headers['EXT']['EXT_Y']['label']
        snr_y_label = snr_y_label.replace('$\mu$', 'u')
        snr_h_label = self.headers['EXT']['EXT_H']['label']
        snr_h_label = snr_h_label.replace('$\mu$', 'u')
        ext_col = 'ext_files'
        tcorr_col = 'tcorr_files'
        # ---------------------------------------------------------------------
        # construct the stats table
        # ---------------------------------------------------------------------
        # columns
        time_series_props['columns'] = TIME_SERIES_COLS[0:8]
        time_series_props['columns'] += [snr_y_label, snr_h_label]
        time_series_props['columns'] += [TIME_SERIES_COLS[8]]
        time_series_props['columns'] += [ext_col, tcorr_col]
        # get values for use in time series table
        for time_series_col in TIME_SERIES_COLS:
            time_series_props[time_series_col] = []
        time_series_props[snr_y_label] = []
        time_series_props[snr_h_label] = []
        time_series_props[ext_col] = []
        time_series_props[tcorr_col] = []
        # get values from self.header_dict
        mjd_vec = np.array(self.header_dict['EXT_MJDMID'])
        seeing_vec = np.array(self.header_dict['EXT_SEEING'])
        airmass_vec = np.array(self.header_dict['EXT_AIRMASS'])
        exptime_vec = np.array(self.header_dict['EXT_EXPTIME'])
        snry_vec = np.array(self.header_dict['EXT_Y'])
        snyh_vec = np.array(self.header_dict['EXT_H'])
        dprtype_vec = np.array(self.header_dict['EXT_DPRTYPE'])
        # get unique object directories (for this object)
        u_obs_dirs = np.unique(self.obs_dirs)
        # loop around observation directories
        for obs_dir in u_obs_dirs:
            # create a mask for this observation directory
            obs_mask = self.obs_dirs == obs_dir
            # get the first and last mjd for this observation directory
            first_mjd = Time(np.min(mjd_vec[obs_mask])).iso
            last_mjd = Time(np.max(mjd_vec[obs_mask])).iso
            # get the number of observations for this observation
            num_obs = str(np.sum(obs_mask))
            # get the seeing for this observation directory
            seeing = np.mean(seeing_vec[obs_mask])
            seeing = '{:.3f}'.format(seeing)
            # get the airmass for this observation directory
            airmass = np.mean(airmass_vec[obs_mask])
            airmass = '{:.3f}'.format(airmass)
            # get the mean exposure time
            exptime = np.mean(exptime_vec[obs_mask])
            exptime =  '{:.3f}'.format(exptime)
            # get the total exposure time
            texptime = np.sum(exptime_vec[obs_mask])
            texptime = '{:.3f}'.format(texptime)
            # get the mean snr_y
            snry = np.mean(snry_vec[obs_mask])
            snry = '{:.3f}'.format(snry)
            # get the mean snr_h
            snyh = np.mean(snyh_vec[obs_mask])
            snyh = '{:.3f}'.format(snyh)
            # get the dprtypes
            dprtype = ','.join(list(np.unique(dprtype_vec[obs_mask])))
            # -----------------------------------------------------------------
            # Create the ext and tellu for this object
            # -----------------------------------------------------------------
            ext_files = np.array(self.ext_files)[obs_mask]
            tcorr_files = np.array(self.tcorr_files)[obs_mask]
            # -----------------------------------------------------------------
            # Create the file lists for this object
            # -----------------------------------------------------------------
            # construct the save path for ext files
            ext_file = f'ext_file_list_{obs_dir}_{self.objname}.txt'
            ext_path = os.path.join(down_save_path, ext_file)
            ext_rel_path = os.path.join(down_rel_path, ext_file)
            create_file_list(ext_files, ext_path)
            ext_download = _make_download('[download]', ext_rel_path)
            ext_value = f'{len(ext_files)} {ext_download}'
            # construct the save path for the tcorr files
            tcorr_file = f'tcorr_file_list_{obs_dir}_{self.objname}.txt'
            tcorr_path = os.path.join(down_save_path, tcorr_file)
            tcorr_rel_path = os.path.join(down_rel_path, tcorr_file)
            create_file_list(tcorr_files, tcorr_path)
            tcorr_download = _make_download('[download]', tcorr_rel_path)
            tcorr_value = f'{len(tcorr_files)} {tcorr_download}'
            # -----------------------------------------------------------------
            # append to the time series properties
            time_series_props[TIME_SERIES_COLS[0]].append(obs_dir)
            time_series_props[TIME_SERIES_COLS[1]].append(first_mjd)
            time_series_props[TIME_SERIES_COLS[2]].append(last_mjd)
            time_series_props[TIME_SERIES_COLS[3]].append(num_obs)
            time_series_props[TIME_SERIES_COLS[4]].append(seeing)
            time_series_props[TIME_SERIES_COLS[5]].append(airmass)
            time_series_props[TIME_SERIES_COLS[6]].append(exptime)
            time_series_props[TIME_SERIES_COLS[7]].append(texptime)
            time_series_props[snr_y_label].append(snry)
            time_series_props[snr_h_label].append(snyh)
            time_series_props[TIME_SERIES_COLS[8]].append(dprtype)
            time_series_props[ext_col].append(ext_value)
            time_series_props[tcorr_col].append(tcorr_value)
        # -----------------------------------------------------------------
        # construct the stats
        # -----------------------------------------------------------------
        # get the stats base name
        time_series_base_name = 'time_series_stat_' + self.objname + '.txt'
        # get the stat path
        stat_path = os.path.join(item_save_path, time_series_base_name)
        # compute the stats
        time_series_stats_table(time_series_props, stat_path)
        # -----------------------------------------------------------------
        # update the paths
        self.time_series_plot_path = None
        self.time_series_stats_table = item_rel_path + time_series_base_name
        self.time_series_dwn_table = None


def add_obj_page(it: int, profile: dict, settings: dict,
                 headers, object_table: Table,
                 file_dict: FileDictReturn) -> Dict[str, Any]:
    # get the object name for this row
    objname = object_table[OBJECT_COLUMN][it]
    # must import here (so that os.environ is set)
    # noinspection PyPep8Naming
    from apero.core import constants
    from apero.core.core import drs_log
    from apero.core.utils import drs_startup
    from apero.tools.module.documentation import drs_markdown
    # get the parameter dictionary of constants from apero
    params = constants.load()
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # get WLOG
    wlog = drs_log.wlog
    # ------------------------------------------------------------------
    # get the parameters for lbl website url
    outdir = os.path.basename(settings['OBJ_OUT'])
    # get the profile name
    name = profile['profile name']
    clean_name = settings['CPN']
    # ------------------------------------------------------------------
    # print progress
    msg = f'\tCreating page for {objname} [{it + 1} of {len(object_table)}]'
    margs = [objname, it + 1, len(object_table)]
    wlog(params, '', msg.format(*margs))
    # create the object class
    object_instance = ObjectData(profile, settings, headers,
                                 objname, file_dict, object_table)
    # ---------------------------------------------------------------------
    # populate the header dictionary for this object instance
    # wlog(params, '', f'\t\tPopulating header dictionary')
    object_instance.populate_header_dict()
    # ---------------------------------------------------------------------
    # generate url for object
    object_url = f'{clean_name}_{outdir}_{objname}_index'
    obj_url = _make_url(objname, object_url)
    # create ARI object page
    object_page = drs_markdown.MarkDownPage(object_url)
    # add title
    object_page.add_title(f'{objname} ({name})')
    # ---------------------------------------------------------------------
    # Add basic text
    # construct text to add
    object_page.add_text(f'This page was last modified: {Time.now()} (UTC)')
    object_page.add_newline()
    # link back to the table
    # create a table page for this table
    table_ref_name = clean_name.lower() + '_' + TABLE_NAMES[0].lower()
    table_ref_url = _make_url('object table', table_ref_name)
    object_page.add_text(f'Back to the {table_ref_url}')
    object_page.add_newline()
    # ---------------------------------------------------------------------
    # table of contents
    # ---------------------------------------------------------------------
    # Add the names of the sections
    names = ['Spectrum', 'LBL', 'CCF', 'Time series']
    # add the links to the pages
    items = [f'spectrum_{clean_name}_objpage_{objname}',
             f'lbl_{clean_name}_objpage_{objname}',
             f'ccf_{clean_name}_objpage_{objname}',
             f'timeseries_{clean_name}_objpage_{objname}']
    # add table of contents
    object_page.add_table_of_contents(items=items, names=names)
    # ---------------------------------------------------------------------
    # Spectrum section
    # ---------------------------------------------------------------------
    # print progress
    # wlog(params, '', f'\t\tCreating spectrum section')
    # add spectrum section
    objpage_spectrum(object_page, names[0], items[0], object_instance)
    # ---------------------------------------------------------------------
    # LBL section
    # ---------------------------------------------------------------------
    # print progress
    # wlog(params, '', f'\t\tCreating LBL section')
    # add LBL section
    objpage_lbl(object_page, names[1], items[1], object_instance)
    # ---------------------------------------------------------------------
    # CCF section
    # ---------------------------------------------------------------------
    # print progress
    # wlog(params, '', f'\t\tCreating CCF section')
    # add CCF section
    objpage_ccf(object_page, names[2], items[2], object_instance)
    # ---------------------------------------------------------------------
    # Time series section
    # ---------------------------------------------------------------------
    # print progress
    # wlog(params, '', f'\t\tCreating time series section')
    # add time series section
    objpage_timeseries(object_page, names[3], items[3], object_instance)
    # ---------------------------------------------------------------------
    # print progress
    # wlog(params, '', f'\t\tWriting to disk')
    # construct a path for the object name
    object_page_path = settings['OBJ_OUT']
    # construct the rst filename
    rst_filename = f'{objname}.rst'
    # save index page
    object_page.write_page(os.path.join(object_page_path, rst_filename))

    # return dictioanry of properties
    rprops = dict()
    rprops['OBJURL'] = obj_url
    # things to return
    return rprops


def add_obj_pages(gsettings: dict, settings: dict, profile: dict,
                  headers: dict, object_table: Table,
                  file_dict: FileDictReturn) -> Table:
    # must import here (so that os.environ is set)
    # noinspection PyPep8Naming
    from apero.core import constants
    from apero.core.core import drs_log
    from apero.core.utils import drs_startup
    # get the parameter dictionary of constants from apero
    params = constants.load()
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # get WLOG
    wlog = drs_log.wlog
    # ------------------------------------------------------------------
    # print progress
    wlog(params, 'info', 'Creating object pages')
    # set up the arguments for the multiprocessing
    args = [0, profile, settings, headers, object_table, file_dict]

    import multiprocessing
    pool = multiprocessing.Pool(processes=gsettings['N_CORES'])
    # use a Manager to create a shared dictionary to collect results
    manager = multiprocessing.Manager()
    results_dict = manager.dict()
    # change the object column to a url
    for it, row in enumerate(object_table):
        # combine arguments
        itargs = [it] + args[1:]
        # run the pool
        results = pool.apply_async(add_obj_page, args=itargs)
        # push result to results storage
        results_dict[it] = results
    # Wait for all jobs to finish
    pool.close()
    pool.join()
    # -------------------------------------------------------------------------
    # Use a dictionary comprehension to extract trhe results from each
    # async results object
    nprops = {key: result.get() for key, result in results_dict.items()}
    # -------------------------------------------------------------------------
    # replace object name with the object name + object url
    object_table[OBJECT_COLUMN] = nprops['OBJURL']
    # remove the FIRST and LAST RAW column
    object_table.remove_column('FIRST_RAW')
    object_table.remove_column('LAST_RAW')
    # return the object table
    return object_table


def objpage_spectrum(page: Any, name: str, ref: str,
                     object_instance: ObjectData):
    # add divider
    # page.add_divider(color=DIVIDER_COLOR, height=DIVIDER_HEIGHT)
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
    if object_instance.spec_plot_path is not None:
        # add the snr plot to the page
        page.add_image(object_instance.spec_plot_path, align='left')
        # add a new line
        page.add_newline()
    # ------------------------------------------------------------------
    # add stats
    if object_instance.spec_stats_table is not None:
        # add the stats table
        page.add_csv_table('', object_instance.spec_stats_table,
                           cssclass='csvtable2')
    # ------------------------------------------------------------------
    # add download links
    if object_instance.spec_dwn_table is not None:
        # add the stats table
        page.add_csv_table('', object_instance.spec_dwn_table,
                           cssclass='csvtable2')


def objpage_lbl(page: Any, name: str, ref: str,
                object_instance: ObjectData):
    # add divider
    # page.add_divider(color=DIVIDER_COLOR, height=DIVIDER_HEIGHT)
    # add a reference to this section
    page.add_reference(ref)
    # ------------------------------------------------------------------
    # deal with no spectrum found
    if len(object_instance.lbl_rdb_files) == 0:
        # add the section heading
        page.add_section(name)
        # print that there is no LBL reduction found
        page.add_text('No LBL reduction found')
        return
    # ------------------------------------------------------------------
    # get the spectrum parameters
    object_instance.get_lbl_parameters()
    # ------------------------------------------------------------------
    # loop around the object+template combinations
    for objcomb in object_instance.lbl_combinations:
        # add subsection for the object+template combination
        page.add_section(f'LBL ({objcomb})')
        # add the lbl plot
        if object_instance.lbl_plot_path[objcomb] is not None:
            # add the snr plot to the page
            page.add_image(object_instance.lbl_plot_path[objcomb], align='left')
            # add a new line
            page.add_newline()
        # ------------------------------------------------------------------
        # add stats
        if object_instance.lbl_stats_table[objcomb] is not None:
            # add the stats table
            page.add_csv_table('', object_instance.lbl_stats_table[objcomb],
                               cssclass='csvtable2')
        # ------------------------------------------------------------------
        # add download links
        if object_instance.lbl_dwn_table is not None:
            # add the stats table
            page.add_csv_table('', object_instance.lbl_dwn_table[objcomb],
                               cssclass='csvtable2')


def objpage_ccf(page: Any, name: str, ref: str, object_instance: ObjectData):
    # add divider
    # page.add_divider(color=DIVIDER_COLOR, height=DIVIDER_HEIGHT)
    # add a reference to this section
    page.add_reference(ref)
    # ------------------------------------------------------------------
    # deal with no spectrum found
    if len(object_instance.ccf_files) == 0:
        # add the section heading
        page.add_section(name)
        # print that there is no LBL reduction found
        page.add_text('No CCF files found')
        return
    # ------------------------------------------------------------------
    # add the section heading
    page.add_section(name)
    # ------------------------------------------------------------------
    # get the spectrum parameters
    object_instance.get_ccf_parameters()
    # ------------------------------------------------------------------
    # add the ccf plot
    if object_instance.ccf_plot_path is not None:
        # add the snr plot to the page
        page.add_image(object_instance.ccf_plot_path, align='left')
        # add a new line
        page.add_newline()
    # ------------------------------------------------------------------
    # add stats
    if object_instance.ccf_stats_table is not None:
        # add the stats table
        page.add_csv_table('', object_instance.ccf_stats_table,
                           cssclass='csvtable2')
    # ------------------------------------------------------------------
    # add download links
    if object_instance.ccf_dwn_table is not None:
        # add the stats table
        page.add_csv_table('', object_instance.ccf_dwn_table,
                           cssclass='csvtable2')


def objpage_timeseries(page: Any, name: str, ref: str,
                       object_instance: ObjectData):
    # add divider
    # page.add_divider(color=DIVIDER_COLOR, height=DIVIDER_HEIGHT)
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
        # add the stats table
        page.add_csv_table('', object_instance.time_series_stats_table,
                           cssclass='csvtable2')
    # ------------------------------------------------------------------
    # add download links
    if object_instance.ccf_dwn_table is not None:
        # add the stats table
        page.add_csv_table('', object_instance.time_series_dwn_table,
                           cssclass='csvtable2')


# =============================================================================
# Plots, stat tables and download tables
# =============================================================================
def download_table(files: List[str], descriptions: List[str],
                   item_path: str, dwn_rel_path: str, down_dir: str,
                   title: str):
    """
    Generic download table saving to the item relative path

    :param files:
    :param descriptions:
    :param dwn_rel_path: the path of the download relative to the page
    :param item_path: the absolute path to the item csv table file (to save to)
    :param down_dir: the path to the download directory
    :param title: the title for the download table
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
    # storage of the last modified times
    last_modified = dict()
    # loop around files and get the paths
    for it, filename in enumerate(files):
        # get the basename
        basename = os.path.basename(filename)
        # get the in path
        in_paths[basename] = filename
        # get the reference path (for the link)
        ref_paths[basename] = dwn_rel_path + basename
        # get the outpath
        out_paths[basename] = os.path.join(down_dir, basename)
        # get the descriptions
        descs[basename] = descriptions[it]
        # get the last modified time of the filename
        if os.path.exists(filename):
            last_mod = Time(os.path.getmtime(filename), format='unix').iso
            last_modified[basename] = last_mod
        else:
            last_modified[basename] = 'N/A'
    # --------------------------------------------------------------------------
    # construct the stats table
    # --------------------------------------------------------------------------
    # start with a download dictionary
    down_dict = dict(Description=[], Value=[], Uploaded=[])
    # loop around files
    for basename in in_paths:
        # add the rdb file
        down_dict['Description'].append(descs[basename])
        down_dict['Value'].append(_make_download(basename, ref_paths[basename]))
        down_dict['Uploaded'].append(last_modified[basename])
        # copy the file from in path to out path
        #   if file is already here then don't copy
        if in_paths[basename] != out_paths[basename]:
            shutil.copy(in_paths[basename], out_paths[basename])
    # --------------------------------------------------------------------------
    # change the columns names
    down_dict2 = dict()
    down_dict2[title] = down_dict['Description']
    down_dict2['File URL'] = down_dict['Value']
    down_dict2['Uploaded'] = down_dict['Uploaded']
    # --------------------------------------------------------------------------
    # convert to table
    down_table = Table(down_dict2)
    # write to file as csv file
    down_table.write(item_path, format='ascii.csv', overwrite=True)


def create_file_list(files: List[str], path: str):
    """
    Writes a list of files to disk
    """
    # if file exists remove it
    if os.path.exists(path):
        os.remove(path)
    # sort files alphabetically
    files = np.sort(files)
    # open file
    with open(path, 'w') as filelist:
        # loop around files
        for filename in files:
            # write to file
            filelist.write(filename + '\n')


def spec_plot(spec_props: Dict[str, Any], plot_path: str, plot_title: str):
    # get parameters from props
    mjd = spec_props['mjd']
    ext_y = spec_props['EXT_Y']
    ext_h = spec_props['EXT_H']
    ext_y_label = spec_props['EXT_Y_LABEL']
    ext_h_label = spec_props['EXT_H_LABEL']
    # --------------------------------------------------------------------------
    # setup the figure
    fig, frame = plt.subplots(1, 1, figsize=(12, 6))
    # set background color
    frame.set_facecolor(PLOT_BACKGROUND_COLOR)
    # --------------------------------------------------------------------------
    # Top plot SNR Y
    # --------------------------------------------------------------------------
    # # plot the CCF RV points
    frame.plot_date(mjd.plot_date, ext_y, fmt='.', alpha=0.5,
                    label=ext_y_label)
    frame.plot_date(mjd.plot_date, ext_h, fmt='.', alpha=0.5,
                    label=ext_h_label)
    frame.legend(loc=0)
    frame.grid(which='both', color='lightgray', ls='--')
    frame.set(xlabel='Date', ylabel='EXT SNR')
    # --------------------------------------------------------------------------
    # add title
    plt.suptitle(plot_title)
    plt.subplots_adjust(hspace=0.1, left=0.1, right=0.99)
    # save figure and close the plot
    plt.savefig(plot_path)
    plt.close()


def spec_stats_table(spec_props: Dict[str, Any], stat_path: str, title: str):
    from apero.core.math import estimate_sigma
    # get parameters from props
    num_raw = spec_props['NUM_RAW_FILES']
    num_pp = spec_props['NUM_PP_FILES']
    num_ext = spec_props['NUM_EXT_FILES']
    num_tcorr = spec_props['NUM_TCORR_FILES']
    ext_y = spec_props['EXT_Y']
    ext_h = spec_props['EXT_H']
    first_raw = spec_props['FIRST_RAW'].iso
    last_raw = spec_props['LAST_RAW'].iso
    first_pp = spec_props['FIRST_PP'].iso
    last_pp = spec_props['LAST_PP'].iso
    first_ext = spec_props['FIRST_EXT'].iso
    last_ext = spec_props['LAST_EXT'].iso
    first_tcorr = spec_props['FIRST_TCORR'].iso
    last_tcorr = spec_props['LAST_TCORR'].iso
    # --------------------------------------------------------------------------
    # Calculate stats
    # --------------------------------------------------------------------------
    # average SNR
    med_snr_y = np.nanmedian(ext_y)
    med_snr_h = np.nanmedian(ext_h)
    # RMS of SNR
    rms_snr_y = estimate_sigma(ext_y)
    rms_snr_h = estimate_sigma(ext_h)
    # --------------------------------------------------------------------------
    # construct the stats table
    # --------------------------------------------------------------------------
    # start with a stats dictionary
    stat_dict = dict(Description=[], Value=[])
    # add number of raw files
    stat_dict['Description'].append('Number raw files [first, last]')
    stat_dict['Value'].append(f'{num_raw} [{first_raw}, {last_raw}]')
    # add number of pp files
    stat_dict['Description'].append('Number pp files [first, last]')
    stat_dict['Value'].append(f'{num_pp} [{first_pp}, {last_pp}]')
    # add number of ext files
    stat_dict['Description'].append('Number ext files [first, last]')
    stat_dict['Value'].append(f'{num_ext} [{first_ext}, {last_ext}]')
    # add number of tcorr files
    stat_dict['Description'].append('Number tcorr files [first, last]')
    stat_dict['Value'].append(f'{num_tcorr} [{first_tcorr}, {last_tcorr}]')
    # add the SNR in Y
    stat_dict['Description'].append('Median SNR Y')
    value = '{:.2f} :math:`\pm` {:.2f} m/s'.format(med_snr_y, rms_snr_y)
    stat_dict['Value'].append(value)
    # add the SNR in H
    stat_dict['Description'].append('Median SNR H')
    value = '{:.2f} :math:`\pm` {:.2f} m/s'.format(med_snr_h, rms_snr_h)
    stat_dict['Value'].append(value)
    # --------------------------------------------------------------------------
    # change the columns names
    stat_dict2 = dict()
    stat_dict2[title] = stat_dict['Description']
    stat_dict2[' '] = stat_dict['Value']
    # --------------------------------------------------------------------------
    # convert to table
    stat_table = Table(stat_dict2)
    # write to file as csv file
    stat_table.write(stat_path, format='ascii.csv', overwrite=True)


def lbl_plot(lbl_props: Dict[str, Any], plot_path: str,
             plot_title: str) -> Dict[str, Any]:
    # setup the figure
    fig, frame = plt.subplots(2, 1, figsize=(12, 6), sharex='all')
    # get parameters from props
    plot_date = lbl_props['plot_date']
    vrad = lbl_props['vrad']
    svrad = lbl_props['svrad']
    snr_h = lbl_props['snr_h']
    snr_h_label = lbl_props['SNR_H_LABEL']
    # set background color
    frame[0].set_facecolor(PLOT_BACKGROUND_COLOR)
    frame[1].set_facecolor(PLOT_BACKGROUND_COLOR)
    # --------------------------------------------------------------------------
    # Top plot LBL RV
    # --------------------------------------------------------------------------
    # plot the points
    frame[0].plot_date(plot_date, vrad, fmt='.', alpha=0.5,
                       color='green', ls='None')
    # plot the error bars
    frame[0].errorbar(plot_date, vrad, yerr=svrad,
                      marker='o', alpha=0.5, color='green', ls='None')
    # find percentile cuts that will be expanded by 150% for the ylim
    pp = np.nanpercentile(vrad, [10, 90])
    diff = pp[1] - pp[0]
    central_val = np.nanmean(pp)
    # used for plotting but also for the flagging of outliers
    ylim = [central_val - 1.5 * diff, central_val + 1.5 * diff]
    # length of the arrow flagging outliers
    l_arrow = (ylim[1] - ylim[0]) / 10.0
    # flag the low outliers
    low = vrad < ylim[0]
    # get the x and y values of the outliers to be looped over within
    # the arrow plotting
    xpoints = np.array(plot_date[low], dtype=float)
    x_range = np.nanmax(plot_date) - np.nanmin(plot_date)
    for ix in range(len(xpoints)):
        frame[0].arrow(xpoints[ix], ylim[0] + l_arrow * 2, 0, -l_arrow,
                       color='red', head_width=0.01 * x_range,
                       head_length=0.25 * l_arrow, alpha=0.5)
    # same as above for the high outliers
    high = vrad > ylim[1]
    xpoints = np.array(plot_date[high], dtype=float)
    for ix in range(len(xpoints)):
        frame[0].arrow(xpoints[ix], ylim[1] - l_arrow * 2, 0, l_arrow,
                       color='red', head_width=0.01 * x_range,
                       head_length=0.25 * l_arrow, alpha=0.5)
    # setting the plot
    frame[0].set(ylim=ylim)
    frame[0].set(title=plot_title)
    frame[0].grid(which='both', color='lightgray', linestyle='--')
    frame[0].set(ylabel='Velocity [m/s]')
    # --------------------------------------------------------------------------
    # Bottom plot SNR
    # --------------------------------------------------------------------------
    # simple plot of the SNR in a sample order. You need to
    # update the relevant ketword for SPIRou
    frame[1].plot_date(plot_date, snr_h, fmt='.',
                       alpha=0.5, color='green', ls='None')
    frame[1].grid(which='both', color='lightgray', linestyle='--')
    frame[1].set(xlabel='Date')
    frame[1].set(ylabel=snr_h_label)
    plt.tight_layout()
    # --------------------------------------------------------------------------
    # save figure and close the plot
    plt.savefig(plot_path)
    plt.close()
    # some parameters are required later save them in a dictionary
    lbl_props['low'] = low
    lbl_props['high'] = high
    lbl_props['ylim'] = ylim
    # return the props
    return lbl_props


def lbl_stats_table(lbl_props: Dict[str, Any], stat_path: str, title: str):
    # get parameters from props
    rjd = lbl_props['rjd']
    vrad = lbl_props['vrad']
    svrad = lbl_props['svrad']
    low = lbl_props['low']
    high = lbl_props['high']
    vel_domain = lbl_props['ylim']
    # --------------------------------------------------------------------------
    # compute the stats
    # --------------------------------------------------------------------------
    # get the 25, 50 and 75 percentile of the velocity uncertainty
    p_sigma = np.nanpercentile(svrad, [25, 50, 75])
    # get the 25, 50 and 75 percentile of the velocity
    v_sigma = np.nanpercentile(abs(vrad - np.nanmedian(vrad)),
                               [25, 50, 75])
    # calculate the number of nights
    n_nights = len(np.unique(np.floor(rjd)))
    # calculate the systemetic velocity
    sys_vel = np.nanmedian(vrad)
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
    stat_dict['Value'].append(len(vrad))
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
    stat_dict['Value'].append('{:.2f} to {:.2f} m/s'.format(*vel_domain))
    # --------------------------------------------------------------------------
    # change the columns names
    stat_dict2 = dict()
    stat_dict2[title] = stat_dict['Description']
    stat_dict2[' '] = stat_dict['Value']
    # --------------------------------------------------------------------------
    # convert to table
    stat_table = Table(stat_dict2)
    # write to file as csv file
    stat_table.write(stat_path, format='ascii.csv', overwrite=True)


def ccf_plot(ccf_props: Dict[str, Any], plot_path: str, plot_title: str):
    # get parameters from props
    mjd = ccf_props['mjd']
    vrad = ccf_props['dv']
    svrad = ccf_props['sdv']
    rv_vec = ccf_props['rv_vec']
    y1_1sig = ccf_props['y1_1sig']
    y2_1sig = ccf_props['y2_1sig']
    y1_2sig = ccf_props['y1_2sig']
    y2_2sig = ccf_props['y2_2sig']
    med_ccf = ccf_props['med_ccf']
    has_fit = ccf_props['has_fit']
    fit = ccf_props['fit']
    xlim = ccf_props['xlim']
    ylim = ccf_props['ylim']
    # --------------------------------------------------------------------------
    # setup the figure
    fig, frame = plt.subplots(3, 1, figsize=(12, 9))
    # set background color
    frame[0].set_facecolor(PLOT_BACKGROUND_COLOR)
    frame[1].set_facecolor(PLOT_BACKGROUND_COLOR)
    frame[2].set_facecolor(PLOT_BACKGROUND_COLOR)
    # --------------------------------------------------------------------------
    # Top plot CCF RV
    # --------------------------------------------------------------------------
    # # plot the CCF RV points
    frame[0].plot_date(mjd.plot_date, vrad, fmt='.', alpha=0.5,
                       color='green')
    # plot the CCF RV errors
    frame[0].errorbar(mjd.plot_date, vrad, yerr=svrad, fmt='o',
                      alpha=0.5, color='green')
    # find percentile cuts that will be expanded by 150% for the ylim
    pp = np.nanpercentile(vrad, [10, 90])
    diff = pp[1] - pp[0]
    central_val = np.nanmean(pp)
    # used for plotting but also for the flagging of outliers
    if diff == 0:
        ylim = [0, 1]
    else:
        ylim = [central_val - 1.5 * diff, central_val + 1.5 * diff]
    # length of the arrow flagging outliers
    l_arrow = (ylim[1] - ylim[0]) / 10.0
    # flag the low outliers
    low = vrad < ylim[0]
    # get the x and y values of the outliers to be looped over within
    # the arrow plotting
    xpoints = np.array(mjd.plot_date[low], dtype=float)
    x_range = np.nanmax(mjd.plot_date) - np.nanmin(mjd.plot_date)
    for ix in range(len(xpoints)):
        frame[0].arrow(xpoints[ix], ylim[0] + l_arrow * 2, 0, -l_arrow,
                       color='red', head_width=0.01 * x_range,
                       head_length=0.25 * l_arrow, alpha=0.5)
    # same as above for the high outliers
    high = vrad > ylim[1]
    xpoints = np.array(mjd.plot_date[high], dtype=float)
    for ix in range(len(xpoints)):
        frame[0].arrow(xpoints[ix], ylim[1] - l_arrow * 2, 0, l_arrow,
                       color='red', head_width=0.01 * x_range,
                       head_length=0.25 * l_arrow, alpha=0.5)
    # setting the plot
    frame[0].set(ylim=ylim)
    frame[0].grid(which='both', color='lightgray', ls='--')
    frame[0].set(xlabel='Date', ylabel='Velocity [m/s]')
    # --------------------------------------------------------------------------
    # Middle plot median CCF
    # --------------------------------------------------------------------------
    # mask by xlim
    limmask = (rv_vec > xlim[0]) & (rv_vec < xlim[1])

    frame[1].fill_between(rv_vec[limmask], y1_2sig[limmask], y2_2sig[limmask],
                          color='orange', alpha=0.5)
    frame[1].fill_between(rv_vec[limmask], y1_1sig[limmask], y2_1sig[limmask],
                          color='red', alpha=0.5)
    frame[1].plot(rv_vec[limmask], med_ccf[limmask], alpha=1.0, color='black')
    if has_fit:
        frame[1].plot(rv_vec[limmask], fit[limmask], alpha=0.8,
                      label='Gaussian fit', ls='--')
    frame[1].legend(loc=0)
    frame[1].set(xlabel='RV [km/s]',
                 ylabel='Normalized CCF')
    frame[1].grid(which='both', color='lightgray', ls='--')

    # --------------------------------------------------------------------------
    # Bottom plot median CCF residuals
    # --------------------------------------------------------------------------
    if has_fit:
        frame[2].fill_between(rv_vec[limmask], y1_2sig[limmask] - fit[limmask],
                              y2_2sig[limmask] - fit[limmask], color='orange',
                              alpha=0.5, label='2-$\sigma$')
        frame[2].fill_between(rv_vec[limmask], y1_1sig[limmask] - fit[limmask],
                              y2_1sig[limmask] - fit[limmask], color='red',
                              alpha=0.5, label='1-$\sigma$')
        frame[2].plot(rv_vec[limmask], med_ccf[limmask] - fit[limmask],
                      alpha=0.8, label='Median residual')
        frame[2].legend(loc=0)
        frame[2].set(xlabel='RV [km/s]', ylabel='Residuals')
    else:
        frame[2].text(0.5, 0.5, 'No fit to CCF possible',
                      horizontalalignment='center')
        frame[2].legend(loc=0)
        frame[2].set(xlim=[0, 1], ylim=[0, 1], xlabel='RV [km/s]',
                     ylabel='Residuals')
    frame[2].grid(which='both', color='lightgray', ls='--')
    # --------------------------------------------------------------------------
    # add title
    plt.suptitle(plot_title)
    plt.subplots_adjust(hspace=0.1, left=0.1, right=0.99)
    # save figure and close the plot
    plt.savefig(plot_path)
    plt.close()


def ccf_stats_table(ccf_props: Dict[str, Any], stat_path: str, title: str):
    from apero.core.math import estimate_sigma
    # get parameters from props
    vrad = ccf_props['dv']
    fwhm = ccf_props['fwhm']
    # --------------------------------------------------------------------------
    # compute the stats
    # --------------------------------------------------------------------------
    # get the systemic velocity
    sys_vel = np.nanmedian(vrad)
    # get the error in systemic velocity
    err_sys_vel = estimate_sigma(vrad)
    # get the fwhm
    ccf_fwhm = np.nanmedian(fwhm)
    # get the error on fwhm
    err_ccf_fwhm = estimate_sigma(fwhm)
    # --------------------------------------------------------------------------
    # construct the stats table
    # --------------------------------------------------------------------------
    # start with a stats dictionary
    stat_dict = dict(Description=[], Value=[])
    # add systemic velocity
    stat_dict['Description'].append('CCF systemic velocity')
    value = '{:.2f} :math:`\pm` {:.2f} m/s'.format(sys_vel, err_sys_vel)
    stat_dict['Value'].append(value)
    # add fwhm
    stat_dict['Description'].append('CCF FWHM')
    value = '{:.2f} :math:`\pm` {:.2f} m/s'.format(ccf_fwhm, err_ccf_fwhm)
    stat_dict['Value'].append(value)
    # add number of files
    stat_dict['Description'].append('Number of CCF files')
    stat_dict['Value'].append(len(ccf_props['files']))
    # --------------------------------------------------------------------------
    # change the columns names
    stat_dict2 = dict()
    stat_dict2[title] = stat_dict['Description']
    stat_dict2[' '] = stat_dict['Value']
    # --------------------------------------------------------------------------
    # convert to table
    stat_table = Table(stat_dict2)
    # write to file as csv file
    stat_table.write(stat_path, format='ascii.csv', overwrite=True)


def time_series_stats_table(time_series_props: Dict[str, Any], stat_path: str):
    # get parameters from props
    columns = time_series_props['columns']
    # --------------------------------------------------------------------------
    # push columns into table
    stat_table = Table()
    for column in columns:
        stat_table[column] = time_series_props[column]
    # write to file as csv file
    stat_table.write(stat_path, format='ascii.csv', overwrite=True)


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


def _make_download(value: str, url: str) -> str:
    """
    Make a download link from a value and a url
    """
    return f':download:`{value} <{url}>`'


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
    ari_gsettings = dict(apero_profiles['settings'])
    header_settings = dict(apero_profiles['headers'])
    del apero_profiles['settings']
    del apero_profiles['headers']
    # deal with a reset
    if ari_gsettings['reset']:
        # remove working directory
        shutil.rmtree(ari_gsettings['working directory'], ignore_errors=True)
    # ----------------------------------------------------------------------
    # step 2: for each profile compile all stats
    all_apero_stats = dict()
    # loop around profiles from yaml file
    for apero_profile_name in apero_profiles:
        # sort out settings
        ari_settings = get_settings(ari_gsettings, apero_profile_name)
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
            apero_stats = compile_stats(ari_gsettings, ari_settings,
                                        apero_profile, header_settings)
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
    ari_settings = get_settings(ari_gsettings)
    # ----------------------------------------------------------------------
    # step 3: write markdown files
    write_markdown(ari_gsettings, ari_settings, all_apero_stats)
    # ----------------------------------------------------------------------
    # step 4: compile sphinx files
    compile_docs(ari_settings)
    # ----------------------------------------------------------------------
    # step 5: upload to hosting
    upload_docs(ari_gsettings, ari_settings)

# =============================================================================
