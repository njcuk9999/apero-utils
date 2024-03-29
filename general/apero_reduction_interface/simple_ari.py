#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-02-20 at 9:35

@author: cook
"""
import argparse
import copy
import glob
import os
import re
import shutil
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd
import yaml
from astropy.table import Table
from astropy.time import Time, TimeDelta
from astropy.io import fits
from astropy import units as uu
from scipy.optimize import curve_fit

# =============================================================================
# Define variables
# =============================================================================
# define which parellisation mode to run in
MULTI = 'Process'
# define the astrometric database column names to get
ASTROMETRIC_COLUMNS = ['OBJNAME', 'RA_DEG', 'DEC_DEG', 'TEFF', 'SP_TYPE']
ASTROMETRIC_DTYPES = [str, float, float, float, str]
# define the log database column names to get
LOG_COLUMNS = ['RECIPE', 'SHORTNAME', 'RUNSTRING', 'START_TIME', 'END_TIME',
               'PASSED_ALL_QC', 'ENDED', 'ERRORMSGS', 'LOGFILE']
LOG_TYPES = ['str', 'str', 'str', 'str', 'str', 'bool', 'bool', 'str', 'str']
# define which instruments have polar
HAS_POLAR = dict(SPIROU=True, NIRPS_HE=False, NIRPS_HA=False)
# define log levels to report
LOG_LEVELS = ['error', 'warning']
# list tables to load
TABLE_NAMES = ['OBJECT_TABLE', 'RECIPE_TABLE']
# Define output path
OUTPUT_PATH = '.'
# define sphinx directories
_DATA_DIR = 'data'
_RST_DIR = 'rst'
_HTML_DIR = '_build/html'
_OUT_DIR = 'output'
_ITEM_DIR = 'items'
_DOWN_DIR = 'downloads'
_OBJ_INDEX_DIR = 'objindex'
# Define output object data dir
_OBJ_OUT_DIR = 'objects'
# define path to local htpasswd file
_PASS_DIR = 'pass'
# define the column which is the object name
OBJECT_COLUMN = 'OBJNAME'

COLUMNS = dict()
COLUMNS['OBJNAME'] = 'OBJNAME'
COLUMNS['RA'] = 'RA [Deg]'
COLUMNS['DEC'] = 'Dec [Deg]'
COLUMNS['TEFF'] = 'Teff [K]'
COLUMNS['SPTYPE'] = 'SpT'
COLUMNS['DPRTYPE'] = 'DPRTYPE'
COLUMNS['RAW'] = 'raw files'
COLUMNS['PP'] = 'pp files'
COLUMNS['EXT'] = 'ext files'
COLUMNS['TCORR'] = 'tcorr files'
COLUMNS['CCF'] = 'ccf files'
COLUMNS['POL'] = 'pol files'
COLUMNS['efits'] = 'e.fits'
COLUMNS['tfits'] = 't.fits'
COLUMNS['vfits'] = 't.fits'
COLUMNS['pfits'] = 'p.fits'
COLUMNS['lbl'] = 'LBL'
COLUMNS['last_obs'] = 'Last observed'
COLUMNS['last_proc'] = 'Last Processed'
# polar columns (only used for instruments that use polarimetry)
POLAR_COLS = ['POL', 'pfits']
# time series column names
TIME_SERIES_COLS = ['Obs Dir', 'First obs mid',
                    'Last obs mid', 'Number of ext', 'Number of tcorr',
                    'Seeing', 'Airmass',
                    'Mean Exptime', 'Total Exptime', 'DPRTYPEs', None, None]
# define the lbl files
LBL_FILETYPES = ['lbl_rdb', 'lbl2_rdb', 'lbl_drift', 'lbl2_drift', 'lbl_fits']
LBL_FILENAMES = ['lbl_{0}_{1}.rdb', 'lbl2_{0}_{1}.rdb',
                 'lbl_{0}_{1}_drift.rdb', 'lbl2_{0}_{1}_drift.rdb',
                 'lbl_{0}_{1}.fits']
LBL_FILE_DESC = ['RDB file', 'RDB2 file', 'Drift file', 'Drift2 file',
                 'LBL RDB fits file']
LBL_DOWNLOAD = [True, True, True, True, False]
# define which one of these files to use as the main plot file
LBL_PLOT_FILE = 0
# define the LBL stat dir (inside the lbl directory)
LBL_STAT_DIR = 'lblstats'
# define the LBL stat files {0} is for the object name + template name
#  i.e. {objname}_{template}
LBL_STAT_FILES = dict()
LBL_STAT_FILES['LBL Diagnostic Plots'] = 'lbl_{0}_plots.pdf'
LBL_STAT_FILES['LBL BERV zp RDB file'] = 'lbl_{0}_bervzp.rdb'
LBL_STAT_FILES['LBL BERV zp RDB2 file'] = 'lbl2_{0}_bervzp.rdb'
LBL_STAT_FILES['LBL BERV Zp Diaganostic Plots'] = 'lbl_{0}_bervzp_plots.pdf'
LBL_STAT_FILES['LBL-PCA RDB file'] = 'lbl_{0}_PCAx.rdb'
LBL_STAT_FILES['LBL-PCA RDB2 file'] = 'lbl2_{0}_PCAx.rdb'
# define how many ccf files to use
MAX_NUM_CCF = 100
# object page styling
DIVIDER_COLOR = '#FFA500'
DIVIDER_HEIGHT = 6
PLOT_BACKGROUND_COLOR = '#FEFDE1'
# define the default column width
DEFAULT_COL_WIDTH = 10
# define table width info
DEFAULT_TABLE_WIDTH = 100
DEFAULT_TABLE_LENGTH = 6
# define rsync commands
RSYNC_CMD_IN = 'rsync -avuz -e "{SSH}" {USER}@{HOST}:{INPATH} {OUTPATH}'
RSYNC_CMD_OUT = 'rsync -avuz -e "{SSH}" {INPATH} {USER}@{HOST}:{OUTPATH}'
# define the html col names for each table
HTML_INCOL_NAMES = dict()
HTML_INCOL_NAMES['OBJECT_TABLE'] = list(COLUMNS.keys())
HTML_INCOL_NAMES['RECIPE_TABLE'] = LOG_COLUMNS
# define the html col names for each table
HTML_OUTCOL_NAMES = dict()
HTML_OUTCOL_NAMES['OBJECT_TABLE'] = list(COLUMNS.values())
HTML_OUTCOL_NAMES['RECIPE_TABLE'] = LOG_COLUMNS
# define the html col types for each table ('str' or 'list')
HTML_COL_TYPES = dict()
HTML_COL_TYPES['OBJECT_TABLE'] = ['str'] * len(HTML_OUTCOL_NAMES['OBJECT_TABLE'])
HTML_COL_TYPES['RECIPE_TABLE'] = ['str'] * len(HTML_OUTCOL_NAMES['RECIPE_TABLE'])
# define the prefilled request link
PREFILLED_REQUEST_LINK = ('https://docs.google.com/forms/d/e/1FAIpQLSev3v7EnHq'
                          '2KyQyVWlDikA8-tDzMINYZA0SkYz93vAanpIdhA/'
                          'viewform?usp=pp_url')
# each entry has an id - we define these in a dictionary
PREFILLED_RDICT = dict()
PREFILLED_RDICT['OBJNAMES'] = 'entry.193319269'
PREFILLED_RDICT['STARTDATE'] = 'entry.293060210'
PREFILLED_RDICT['ENDDATE'] = 'entry.299959352'
PREFILLED_RDICT['APERO_MODE'] = 'entry.923271952'
PREFILLED_RDICT['FIBER'] = 'entry.1651377065'
PREFILLED_RDICT['DPRTYPES'] = 'entry.298205572'
PREFILLED_RDICT['DRSOUTID'] = 'entry.459458944'
# define the request form fiber types
RDICT_FIBERS = ['Science fiber', 'Reference fiber', 'All fibers']
# define science dprtypes
SCIENCE_DPRTYPES = ['OBJ_FP', 'OBJ_SKY', 'OBJ_DARK', 'POLAR_FP', 'POLAR_DARK']


# =============================================================================
# Classes
# =============================================================================
class FileType:
    def __init__(self, name: str, block_kind: Optional[str] = None,
                 kw_output: Optional[str] = None, fiber: Optional[str] = None,
                 count: bool = True, chain: Optional[str] = None):
        self.name = str(name)
        self.block_kind = copy.deepcopy(block_kind)
        self.kw_output = copy.deepcopy(kw_output)
        self.fiber = copy.deepcopy(fiber)
        self.count = bool(count)
        self.chain = copy.deepcopy(chain)

        self.cond = None
        self.files = []
        self.qc_mask = None
        self.obsdirs = []
        self.processed_times = []
        self.first = None
        self.last = None
        self.num = 0
        self.num_passed = 0
        self.num_failed = 0

    def copy_new(self) -> 'FileType':
        # inherit values from self
        new = FileType(self.name, self.block_kind, self.kw_output, self.fiber,
                       self.count)
        return new

    def count_files(self, objname, indexdbm, logdbm,
                    filetypes: Dict[str, 'FileType']):
        # deal with not wanting to count in every case
        if not self.count:
            return
            # ------------------------------------------------------------------
        # Construct the condition for the query
        # ------------------------------------------------------------------
        self.cond = f'KW_OBJNAME="{objname}"'
        if self.block_kind is not None:
            self.cond += f' AND BLOCK_KIND="{self.block_kind}"'
        if self.kw_output is not None:
            self.cond += f' AND KW_OUTPUT="{self.kw_output}"'
        if self.fiber is not None:
            self.cond += f' AND KW_FIBER="{self.fiber}"'
        # ------------------------------------------------------------------
        # run counting conditions using indexdbm
        # ------------------------------------------------------------------
        # deal with chain value being zero (don't continue)
        if self.chain is not None:
            if self.chain in filetypes:
                # get the number of files in chain filetype
                if filetypes[self.chain].num:
                    return
        # get the files
        dbcols = 'ABSPATH,OBS_DIR,KW_PID,KW_MID_OBS_TIME,KW_DRS_DATE_NOW'
        findex_table = indexdbm.get_entries(dbcols, condition=self.cond)
        # get a mask of rows that passed QC (based on PID)
        if self.name != 'raw':
            mask = _filter_pids(findex_table, logdbm)
            self.files = np.array(findex_table['ABSPATH'])
            self.qc_mask = mask
            self.obsdirs = np.array(findex_table['OBS_DIR'])
            # add last processed time of all files for this object
            pdates = np.array(findex_table['KW_DRS_DATE_NOW']).astype(str)
            if len(pdates) > 0:
                times_it = np.max(Time(pdates, format='iso'))
                self.processed_times.append(times_it)
        else:
            self.files = np.array(findex_table['ABSPATH'])
            self.qc_mask = np.ones(len(self.files)).astype(bool)
            self.obsdirs = np.array(findex_table['OBS_DIR'])
        # get the first and last files in time
        mjdmids = np.array(findex_table['KW_MID_OBS_TIME']).astype(float)
        mjdmids = Time(mjdmids, format='mjd')
        if len(mjdmids) > 0:
            self.first = np.min(mjdmids)
            self.last = np.max(mjdmids)
        # count files
        self.num = len(self.files)
        self.num_passed = np.sum(self.qc_mask)
        self.num_failed = np.sum(~self.qc_mask)

    def get_files(self, qc: Optional[bool] = None, attr='files'):
        # get the value of the attribute
        vector = getattr(self, attr)
        # deal with no qc --> all files
        if qc is None:
            return vector
        # deal with no mask --> all files
        if self.qc_mask is None:
            return vector
        # deal with qc = True --> return only qc files
        if qc:
            return vector[self.qc_mask]
        # deal with qc = False --> return only non qc files
        else:
            return vector[~self.qc_mask]


class ObjectData:
    def __init__(self, objname, filetypes: Dict[str, FileType]):
        # get the object name
        self.objname = objname
        self.ra: Optional[float] = None
        self.dec: Optional[float] = None
        self.teff: Optional[float] = None
        self.sptype: Optional[str] = None
        self.dprtypes: Optional[str] = None
        # ---------------------------------------------------------------------
        # Add files as copy of filetype class
        # ---------------------------------------------------------------------
        self.filetypes: Dict[str, FileType] = dict()
        for key in filetypes:
            self.filetypes[key] = filetypes[key].copy_new()
        # flag for whether we have polar files
        self.has_polar = False
        # ---------------------------------------------------------------------
        # lbl parameters
        self.lbl_templates: Optional[List[str]] = None
        self.lbl_select: Optional[str] = None
        # add the lbl stat files into a dictionary
        self.lbl_stat_files = dict()
        # loop around lbl stats files and load lblfilekey dict in
        #   each dict should contain obj+temp combinations
        for lblfilekey in LBL_STAT_FILES:
            self.lbl_stat_files[lblfilekey] = dict()
        # ---------------------------------------------------------------------
        # last processed of all files for this object
        self.last_processed: Optional[Time] = None
        # ---------------------------------------------------------------------
        # get spectrum output parameters (for page integration)
        self.spec_plot_path: Optional[str] = None
        self.spec_stats_table: Optional[str] = None
        self.spec_rlink_table: Optional[str] = None
        self.spec_dwn_table: Optional[str] = None
        # get lbl output parameters (for page integration)
        self.lbl_combinations = []
        self.lbl_plot_path = dict()
        self.lbl_stats_table = dict()
        self.lbl_rlink_table = dict()
        self.lbl_dwn_table = dict()
        # get ccf output parameters (for page integration)
        self.ccf_plot_path: Optional[str] = None
        self.ccf_stats_table: Optional[str] = None
        self.ccf_rlink_table: Optional[str] = None
        self.ccf_dwn_table: Optional[str] = None
        # get time series output parameters (for page integration)
        self.time_series_plot_path: Optional[str] = None
        self.time_series_stats_table: Optional[str] = None
        self.time_series_rlink_table: Optional[str] = None
        self.time_series_dwn_table: Optional[str] = None
        # ---------------------------------------------------------------------
        # these are only added when we need them
        self.profile: Optional[Dict[str, Any]] = None
        self.settings: Optional[Dict[str, Any]] = None
        self.gsettings: Optional[Dict[str, Any]] = None
        self.headers: Optional[Dict[str, Any]] = None
        # ---------------------------------------------------------------------
        # parameters for the object index
        self.objurl: Optional[str] = None
        self.objpageref: Optional[str] = None
        # ---------------------------------------------------------------------
        # store the required header info
        self.header_dict = dict()

    def add_astrometrics(self, table):
        self.ra = table['RA_DEG']
        self.dec = table['DEC_DEG']
        self.teff = table['TEFF']
        self.sptype = table['SP_TYPE']

    def add_files_stats(self, indexdbm, logdbm):
        # loop around raw files
        for key in self.filetypes:
            # get iterations filetype
            filetype = self.filetypes[key]
            # count files and get lists of files
            filetype.count_files(self.objname, indexdbm, logdbm, self.filetypes)
            # if there are no entries we have no raw files for this object
            if key == 'raw':
                if filetype.num == 0:
                    return
        # ------------------------------------------------------------------
        # Add a dpr type column
        dprtypes = indexdbm.get_entries('KW_DPRTYPE',
                                        condition=self.filetypes['pp'].cond)
        self.dprtypes = ','.join(list(np.unique(dprtypes)))
        # ------------------------------------------------------------------
        # get all filetype last processing times
        all_last_processed = []
        for key in self.filetypes:
            if len(self.filetypes[key].processed_times) > 0:
                for _time in self.filetypes[key].processed_times:
                    all_last_processed.append(_time)
        # convert to Time
        if len(all_last_processed) > 0:
            all_last_processed = Time(np.array(all_last_processed))
            # get the last processed time of all files
            self.last_processed = np.max(all_last_processed)
        else:
            self.last_processed = None

    def add_settings(self, profile: dict, settings: dict,
                     gsettings: dict, headers: dict):
        self.profile = profile
        self.settings = settings
        self.gsettings = gsettings
        self.headers = headers

    def populate_header_dict(self):
        """
        Populate the header dictionary with the required header keys
        :return:
        """
        # loop around COUNT COLS and populate header dict
        for key in self.filetypes:
            # check file kind in headers
            if key not in self.headers:
                continue
            # get the header keys
            self.get_header_keys(self.headers[key], self.filetypes[key].files)

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

    def rlink(self, filetype: Optional[str] = None,
              startdate: Optional[str] = None,
              enddate: Optional[str] = None,
              apero_mode: Optional[str] = None,
              fiber: Optional[str] = None,
              dprtypes: Optional[str] = None,
              drsoutid: Optional[str] = None):
        """
        Create a link to pre-fill the request form

        :param filetype: str, the filetype (FileType.name) to use
        :param startdate: optional str, the start date DD/MM/YYYY
        :param enddate: optional str, the end date DD/MM/YYYY
        :param apero_mode: optional str, the apero mode to use (e.g.
                           nirps_he_online)
        :param fiber: optional str, the fibers to get should be one of
                      "Science fiber", "Reference fiber", "All fibers"
        :param dprtypes: optional str, the dprtypes
        :param drsoutid: optional str, the drsoutid (not used if filetype is
                         defined)
        :return:
        """
        # lets create a url
        url = PREFILLED_REQUEST_LINK
        # get the objname
        url += _url_addp(PREFILLED_RDICT['OBJNAMES'], self.objname)
        # add the start date
        if startdate is not None:
            url += _url_addp(PREFILLED_RDICT['STARTDATE'], startdate)
        # add the end date
        if enddate is not None:
            url += _url_addp(PREFILLED_RDICT['ENDDATE'], enddate)
        # add the apero mode
        if apero_mode is not None:
            url += _url_addp(PREFILLED_RDICT['APERO_MODE'], apero_mode)
        # add the fiber
        if fiber is not None and fiber in RDICT_FIBERS:
            url += _url_addp(PREFILLED_RDICT['FIBER'], fiber)
        # add the dprtypes
        if dprtypes is not None:
            url += _url_addp(PREFILLED_RDICT['DPRTYPES'], dprtypes)
        else:
            dprtypes = np.char.array(SCIENCE_DPRTYPES)
            dprtypes = list(dprtypes.strip())
            url += _url_addp(PREFILLED_RDICT['DPRTYPES'], dprtypes)
        # add the drsoutids using filetype
        if filetype in self.filetypes:
            # get the file instance class
            fileinst = self.filetypes[filetype]
            # add the drsoutids
            url += _url_addp(PREFILLED_RDICT['DRSOUTID'], fileinst.kw_output)
        elif drsoutid is not None:
            # add the drsoutids
            url += _url_addp(PREFILLED_RDICT['DRSOUTID'], drsoutid)
        # return the url
        return url

    def get_spec_parameters(self):
        #
        ext_files = self.filetypes['ext'].get_files()
        # don't go here if ext files are not present
        if len(ext_files) == 0:
            return
        # ---------------------------------------------------------------------
        # generate place to save figures
        item_save_path = self.settings['ITEMS']
        item_rel_path = f'../{_ITEM_DIR}/'
        down_save_path = self.settings['DOWNS']
        down_rel_path = f'../{_DOWN_DIR}/'
        # ---------------------------------------------------------------------
        # storage for spectrum values
        spec_props = dict()
        # get files to use
        spec_props['RAW'] = self.filetypes['raw']
        spec_props['EXT'] = self.filetypes['ext']
        spec_props['TCORR'] = self.filetypes['tcorr']
        spec_props['S1D'] = self.filetypes['s1d']
        spec_props['SC1D'] = self.filetypes['sc1d']
        # ---------------------------------------------------------------------
        # object properties
        spec_props['COORD_URL'] = f'indexing_{self.objname}'
        spec_props['RA'] = self.ra
        spec_props['Dec'] = self.dec
        spec_props['Teff'] = self.teff
        spec_props['Spectral Type'] = self.sptype
        spec_props['DPRTYPES'] = self.dprtypes
        # ---------------------------------------------------------------------
        # header dict alias
        hdict = self.header_dict
        ftypes = self.filetypes
        # get values for use in plot
        spec_props['mjd'] = Time(np.array(hdict['EXT_MJDMID']))
        spec_props['EXT_Y'] = np.array(hdict['EXT_Y'])
        ext_h = np.array(hdict['EXT_H'])
        spec_props['EXT_H'] = ext_h
        spec_props['EXT_Y_LABEL'] = self.headers['ext']['EXT_Y']['label']
        spec_props['EXT_H_LABEL'] = self.headers['ext']['EXT_H']['label']
        spec_props['NUM_RAW_FILES'] = ftypes['raw'].num_passed
        spec_props['NUM_PP_FILES'] = ftypes['pp'].num_passed
        spec_props['NUM_EXT_FILES'] = ftypes['ext'].num_passed
        spec_props['NUM_TCORR_FILES'] = ftypes['tcorr'].num_passed
        spec_props['NUM_PP_FILES_FAIL'] = ftypes['pp'].num_failed
        spec_props['NUM_EXT_FILES_FAIL'] = ftypes['ext'].num_failed
        spec_props['NUM_TCORR_FILES_FAIL'] = ftypes['tcorr'].num_failed
        # -----------------------------------------------------------------
        # add the first and last raw file type
        first_time = self.filetypes['raw'].first
        last_time = self.filetypes['raw'].last
        if first_time is not None:
            spec_props['FIRST_RAW'] = first_time.iso
        else:
            spec_props['FIRST_RAW'] = None
        if last_time is not None:
            spec_props['LAST_RAW'] = last_time.iso
        else:
            spec_props['LAST_RAW'] = None
        # Add first / last pp files
        if ftypes['pp'].num_passed > 0:
            spec_props['FIRST_PP'] = Time(np.min(hdict['PP_MJDMID'])).iso
            spec_props['LAST_PP'] = Time(np.max(hdict['PP_MJDMID'])).iso
            spec_props['LAST_PP_PROC'] = Time(np.max(hdict['PP_PROC'])).iso
            spec_props['PP_VERSION'] = ','.join(list(np.unique(hdict['PP_VERSION'])))
        else:
            spec_props['FIRST_PP'] = None
            spec_props['LAST_PP'] = None
            spec_props['LAST_PP_PROC'] = None
            spec_props['PP_VERSION'] = None
            # Add first / last ext files
        if ftypes['ext'].num_passed > 0:
            spec_props['FIRST_EXT'] = Time(np.min(hdict['EXT_MJDMID'])).iso
            spec_props['LAST_EXT'] = Time(np.max(hdict['EXT_MJDMID'])).iso
            spec_props['LAST_EXT_PROC'] = Time(np.max(hdict['EXT_PROC'])).iso
            spec_props['EXT_VERSION'] = ','.join(list(np.unique(hdict['EXT_VERSION'])))
        else:
            spec_props['FIRST_EXT'] = None
            spec_props['LAST_EXT'] = None
            spec_props['LAST_EXT_PROC'] = None
            spec_props['EXT_VERSION'] = None
        # Add first / last tcorr files
        if ftypes['tcorr'].num_passed > 0:
            spec_props['FIRST_TCORR'] = Time(np.min(hdict['TCORR_MJDMID'])).iso
            spec_props['LAST_TCORR'] = Time(np.max(hdict['TCORR_MJDMID'])).iso
            spec_props['LAST_TCORR_PROC'] = Time(np.max(hdict['TCORR_PROC'])).iso
            spec_props['TCORR_VERSION'] = ','.join(list(np.unique(hdict['TCORR_VERSION'])))
        else:
            spec_props['FIRST_TCORR'] = None
            spec_props['LAST_TCORR'] = None
            spec_props['LAST_TCORR_PROC'] = None
            spec_props['TCORR_VERSION'] = None
        # -----------------------------------------------------------------
        # standard request keyword args
        rkwargs = dict(fiber='Science fiber',
                       dprtypes=SCIENCE_DPRTYPES,
                       apero_mode=self.settings['CPN'])
        # add the links to request data
        spec_props['RLINK_EXT_E2DSFF'] = self.rlink(filetype='ext', **rkwargs)
        spec_props['RLINK_EXT_S1D_V'] = self.rlink(filetype='s1d', **rkwargs)
        spec_props['RLINK_TELLU_OBJ'] = self.rlink(filetype='tcorr', **rkwargs)
        spec_props['RLINK_TELLU_S1DV'] = self.rlink(filetype='sc1d', **rkwargs)
        spec_props['RLINK_TELLU_TEMP'] = self.rlink(drsoutid='TELLU_TEMP',
                                                    **rkwargs)
        spec_props['RLINK_TELLU_TEMP_S1D'] = self.rlink(drsoutid='TELLU_TEMP_S1DV',
                                                        **rkwargs)
        spec_props['RLINK_DRS_POST_E'] = self.rlink(filetype='efiles', **rkwargs)
        spec_props['RLINK_DRS_POST_T'] = self.rlink(filetype='tfiles', **rkwargs)
        spec_props['RLINK_DRS_POST_S'] = self.rlink(drsoutid='DRS_POST_S',
                                                    **rkwargs)
        spec_props['RLINK_DRS_POST_V'] = self.rlink(filetype='vfiles', **rkwargs)
        spec_props['RLINK_DRS_POST_P'] = self.rlink(filetype='pfiles', **rkwargs)
        # -----------------------------------------------------------------
        # we have to match files (as ext_files, tcorr_files and raw_files may
        #   be different lengths)
        matched = False
        # get the median snr
        med_snr = np.nanmedian(ext_h)
        n_ext_h = abs(ext_h - med_snr)
        # sort all snr by closest to the median
        all_snr_pos = list(np.argsort(n_ext_h))
        # set these up
        pos_ext, pos_raw, pos_s1d, pos_sc1d = None, None, None, None
        file_ext = 'NoFile'
        # loop until we match
        while not matched and len(all_snr_pos) > 0:
            # Find the closest to the median
            pos_ext = all_snr_pos[0]
            # set the filename
            file_ext = spec_props['EXT'].get_files()[pos_ext]
            # Find the matching raw file
            pos_raw = _match_file(reffile=file_ext,
                                  files=spec_props['RAW'].get_files(qc=True))
            pos_s1d = _match_file(reffile=file_ext,
                                  files=spec_props['S1D'].get_files(qc=True))
            pos_sc1d = _match_file(reffile=file_ext,
                                   files=spec_props['SC1D'].get_files(qc=True))
            # three conditions that have to be met for a match
            cond1 = pos_raw is not None
            cond2 = pos_s1d is not None
            cond3 = pos_sc1d is not None
            # we only stop is a match is found
            if cond1 and cond2 and cond3:
                matched = True
            else:
                all_snr_pos = all_snr_pos[1:]
        # ---------------------------------------------------------------------
        cond1 = pos_raw is None
        cond2 = pos_s1d is None
        cond3 = pos_sc1d is None
        # deal with case we have no matching raw file we have a big problem
        #   extracted file cannot exist without a raw file
        if cond1 or cond2 or cond3:
            if cond1:
                print(f'No raw file matching {file_ext}. '
                      f'This should not be possible')
            if cond2 or cond3:
                print(f'No s1d file matching {file_ext}')
            return
        # ---------------------------------------------------------------------
        # get the extracted spectrum for the spectrum with the highest SNR
        ext_file = spec_props['S1D'].get_files(qc=True)[pos_s1d]
        ext_table = load_table(ext_file, hdu=1)
        # get wavelength masks for plotting
        wavemap = ext_table['wavelength']
        limits = self.gsettings['SpecWave']
        wavemask0 = (wavemap > limits['limit0'][0])
        wavemask0 &= (wavemap < limits['limit0'][1])
        wavemask1 = (wavemap > limits['limit1'][0])
        wavemask1 &= (wavemap < limits['limit1'][1])
        wavemask2 = (wavemap > limits['limit2'][0])
        wavemask2 &= (wavemap < limits['limit2'][1])
        wavemask3 = (wavemap > limits['limit3'][0])
        wavemask3 &= (wavemap < limits['limit3'][1])
        # ---------------------------------------------------------------------
        # push into spec_props
        spec_props['WAVE'] = np.array(ext_table['wavelength'])
        spec_props['EXT_SPEC'] = np.array(ext_table['flux'])
        spec_props['EXT_SPEC_ERR'] = np.array(ext_table['eflux'])
        spec_props['WAVEMASK0'] = wavemask0
        spec_props['WAVEMASK1'] = wavemask1
        spec_props['WAVEMASK2'] = wavemask2
        spec_props['WAVEMASK3'] = wavemask3
        spec_props['WAVELIM0'] = limits['limit0']
        spec_props['WAVELIM1'] = limits['limit1']
        spec_props['WAVELIM2'] = limits['limit2']
        spec_props['WAVELIM3'] = limits['limit3']
        spec_props['MAX_SNR'] = np.round(spec_props['EXT_H'][pos_ext], 2)
        raw_file = spec_props['RAW'].get_files(qc=True)[pos_raw]
        spec_props['MAX_FILE'] = os.path.basename(raw_file)
        # ---------------------------------------------------------------------
        # deal with having telluric file
        if pos_sc1d is not None:
            # get the telluric corrected spectrum for the spectrum with the
            # highest SNR
            tcorr_file = spec_props['SC1D'].get_files(qc=True)[pos_sc1d]
            tcorr_table = load_table(tcorr_file, hdu=1)
            # push into spec_props
            spec_props['TCORR_SPEC'] = np.array(tcorr_table['flux'])
            spec_props['TCORR_SPEC_ERR'] = np.array(tcorr_table['eflux'])
        else:
            spec_props['TCORR_SPEC'] = None
            spec_props['TCORR_SPEC_ERR'] = None
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
        spec_stats_table(spec_props, stat_path, title='Spectrum Information')
        # -----------------------------------------------------------------
        # construct the header file
        ext_header_file = os.path.join(down_save_path,
                                       f'ext2d_header_{self.objname}_file.csv')
        create_header_file(spec_props['EXT'].get_files(qc=True),
                           self.headers, 'ext', self.header_dict,
                           ext_header_file)
        # -----------------------------------------------------------------
        # Create the request link table for this object
        # -----------------------------------------------------------------
        # get the rlink base name
        rlink_base_name = 'spec_rlink_' + self.objname + '.txt'
        # get the rlink table path
        rlink_item_path = os.path.join(item_save_path, rlink_base_name)
        # define the keys in spec_props that contain rlinks to add
        rlinks = ['RLINK_EXT_E2DSFF', 'RLINK_EXT_S1D_V',
                  'RLINK_TELLU_OBJ', 'RLINK_TELLU_S1DV',
                  'RLINK_TELLU_TEMP', 'RLINK_TELLU_TEMP_S1D',
                  'RLINK_DRS_POST_E', 'RLINK_DRS_POST_T',
                  'RLINK_DRS_POST_S', 'RLINK_DRS_POST_V',
                  'RLINK_DRS_POST_P']
        # define the rlink descriptions
        rlink_text = ['Extracted 2D spectra', 'Extracted 1D spectra',
                      'Telluric corrected 2D spectra',
                      'Telluric corrected 1D spectra',
                      '2D Template (after telluric correction)',
                      '1D Template (after telluric correction)',
                      'Packaged extraction files (e.fits)',
                      'Packaged telluric corrected files (t.fits)',
                      'Packaged s1d extract+tcorr (s.fits)',
                      'Packaged velocity files (v.fits)',
                      'Packaged polarimetry files (p.fits)']
        # compute the rlink table
        create_request_link_table(spec_props, rlink_item_path, rlinks,
                                  rlink_text)
        # -----------------------------------------------------------------
        # Create the file lists for this object
        # -----------------------------------------------------------------
        # construct the save path for ext files (2D)
        ext2d_file = os.path.join(down_save_path,
                                  f'ext2d_{self.objname}_file_list.txt')
        create_file_list(spec_props['EXT'].get_files(qc=True), ext2d_file)
        # construct the save path for ext files (1D)
        ext1d_file = os.path.join(down_save_path,
                                  f'ext1d_{self.objname}_file_list.txt')
        create_file_list(spec_props['S1D'].get_files(qc=True), ext1d_file)
        # construct the save path for the tcorr files (2D)
        tcorr2d_file = os.path.join(down_save_path,
                                    f'tcorr2d_{self.objname}_file_list.txt')
        create_file_list(spec_props['TCORR'].get_files(qc=True), tcorr2d_file)
        # construct the save path for the tcorr files (1D)
        tcorr1d_file = os.path.join(down_save_path,
                                    f'tcorr1d_{self.objname}_file_list.txt')
        create_file_list(spec_props['SC1D'].get_files(qc=True), tcorr1d_file)
        # -----------------------------------------------------------------
        # construct the download table
        # -----------------------------------------------------------------
        # get the download base name
        dwn_base_name = 'spec_download_' + self.objname + '.txt'
        # get the download table path
        dwn_item_path = os.path.join(item_save_path, dwn_base_name)
        # define the download files
        down_files = [ext2d_file, ext1d_file, tcorr2d_file, tcorr1d_file,
                      ext_header_file]
        # define the download descriptions
        down_descs = ['Extracted 2D spectra', 'Extracted 1D spectra',
                      'Telluric corrected 2D spectra',
                      'Telluric corrected 1D spectra',
                      'Extracted 2D header file']
        # compute the download table
        download_table(down_files, down_descs, dwn_item_path, down_rel_path,
                       down_save_path, title='Spectrum Downloads')
        # -----------------------------------------------------------------
        # update the paths
        self.spec_plot_path = item_rel_path + plot_base_name
        self.spec_stats_table = item_rel_path + stat_base_name
        self.spec_rlink_table = item_rel_path + rlink_base_name
        self.spec_dwn_table = item_rel_path + dwn_base_name

    def get_lbl_parameters(self):
        """
        Get the LBL parameters for this object

        :return:
        """
        # get lbl rdb files
        lbl_files = dict()
        for filetype in LBL_FILETYPES:
            if self.filetypes[filetype].num > 0:
                lbl_files[filetype] = self.filetypes[filetype].get_files()

        # don't go here is lbl rdb files are not present
        if len(lbl_files) == 0:
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
        # loop around templates (found previously)
        for template in self.lbl_templates:
            # make objtmp string
            lbl_objtmp = f'{self.objname}_{template}'
            # set each template to have a dictionary of files
            lbl_objtmps[lbl_objtmp] = dict()
            # loop around each lbl file type
            for lbl_filetype in lbl_files:
                # file to be found
                matched_file = None
                # find the file that matches the template
                for lbl_file in lbl_files[lbl_filetype]:
                    if lbl_objtmp in lbl_file:
                        matched_file = lbl_file
                        break
                # append to list if we found a matching objname+template file
                if matched_file is not None:
                    lbl_objtmps[lbl_objtmp][lbl_filetype] = matched_file
        # ---------------------------------------------------------------------
        # def the plot file
        plot_file = LBL_FILETYPES[LBL_PLOT_FILE]
        # loop around the objname+template combinations
        for lbl_objtmp in lbl_objtmps:
            # deal with the plot file not existing
            if plot_file not in lbl_objtmps[lbl_objtmp]:
                # must set these to None if no LBL files
                self.lbl_plot_path[lbl_objtmp] = None
                self.lbl_stats_table[lbl_objtmp] = None
                self.lbl_rlink_table[lbl_objtmp] = None
                self.lbl_dwn_table[lbl_objtmp] = None
                continue
            # get the plot file for this objname+template
            lbl_pfile = lbl_objtmps[lbl_objtmp][plot_file]
            # load rdb file
            rdb_table = load_table(lbl_pfile, format='ascii.rdb')
            # get the values required
            lbl_props['rjd'] = np.array(rdb_table['rjd'])
            lbl_props['vrad'] = np.array(rdb_table['vrad'])
            lbl_props['svrad'] = np.array(rdb_table['svrad'])
            lbl_props['plot_date'] = np.array(rdb_table['plot_date'])
            lbl_props['snr_h'] = np.array(rdb_table[ext_h_key])
            lbl_props['SNR_H_LABEL'] = self.headers['LBL']['EXT_H']['label']
            lbl_props['RESET_RV'] = np.array(rdb_table['RESET_RV']).astype(bool)
            lbl_props['NUM_RESET_RV'] = np.sum(rdb_table['RESET_RV'])
            # get the lbl header key
            lbl_version_hdrkey = self.headers['LBL']['LBL_VERSION']['key']
            # find the lbl fits file
            for it, lblfilekey in enumerate(LBL_FILETYPES):
                if lblfilekey == 'lbl_fits':
                    # get the lbl file
                    lbl_file = lbl_objtmps[lbl_objtmp][lblfilekey]
                    # get the header
                    lbl_hdr = fits.getheader(lbl_file)
                    # if we have the lbl version header key add it
                    if lbl_version_hdrkey in lbl_hdr:
                        lbl_props['version'] = lbl_hdr[lbl_version_hdrkey]
            # -----------------------------------------------------------------
            # standard request keyword args
            rkwargs = dict(fiber='Science fiber',
                           dprtypes=SCIENCE_DPRTYPES,
                           apero_mode=self.settings['CPN'])
            # add the links to request data
            lbl_props['RLINK_LBL_FITS'] = self.rlink(filetype='lbl.fits',
                                                     **rkwargs)
            for filetype in LBL_FILETYPES:
                lbl_props[f'RLINK_{filetype}'] = self.rlink(filetype=filetype,
                                                            **rkwargs)

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
            # Create the request link table for this object
            # -----------------------------------------------------------------
            # get the rlink base name
            rlink_base_name = 'lbl_rlink_' + self.objname + '.txt'
            # get the rlink table path
            rlink_item_path = os.path.join(item_save_path, rlink_base_name)
            # define the keys in spec_props that contain rlinks to add
            rlinks = ['RLINK_LBL_FITS']
            for filetype in LBL_FILETYPES:
                rlinks += [f'RLINK_{filetype}']

            # define the rlink descriptions
            rlink_text = ['LBL fits files']
            for filetype in LBL_FILETYPES:
                rlink_text += [f'{filetype} files']
            # compute the rlink table
            create_request_link_table(lbl_props, rlink_item_path, rlinks,
                                      rlink_text)
            # -----------------------------------------------------------------
            # construct the download table
            # -----------------------------------------------------------------
            # get the download base name
            dwn_base_name = 'lbl_download_' + lbl_objtmp + '.txt'
            # get the download table path
            item_path = os.path.join(item_save_path, dwn_base_name)
            # define the download files
            down_files = []
            # define the download descriptions
            down_descs = []
            # add the lbl files to the down files
            for it, lblfilekey in enumerate(LBL_FILETYPES):
                # check that we want to add to download files
                if not LBL_DOWNLOAD[it]:
                    continue
                # check for file in lbl_objtmps
                if lblfilekey not in lbl_objtmps[lbl_objtmp]:
                    continue
                # get lbl file
                lbl_file = lbl_objtmps[lbl_objtmp][lblfilekey]
                # add to down_files
                down_files.append(lbl_file)
                down_descs.append(LBL_FILE_DESC[it])
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
            self.lbl_rlink_table[lbl_objtmp] = item_rel_path + rlink_base_name
            self.lbl_dwn_table[lbl_objtmp] = item_rel_path + dwn_base_name
        # ---------------------------------------------------------------------
        # set the lbl combinations
        self.lbl_combinations = lbl_objtmps.keys()

    def get_ccf_parameters(self):

        from apero.core.math import normal_fraction
        # get ccf files
        ccf_files = self.filetypes['ccf'].get_files(qc=True)
        # don't go here is lbl rdb files are not present
        if len(ccf_files) == 0:
            return
        # ---------------------------------------------------------------------
        # generate place to save figures
        item_save_path = self.settings['ITEMS']
        item_rel_path = f'../{_ITEM_DIR}/'
        down_save_path = self.settings['DOWNS']
        down_rel_path = f'../{_DOWN_DIR}/'
        # -----------------------------------------------------------------
        # alias to header dict
        hdict = self.header_dict
        # storage for ccf values
        ccf_props = dict()
        # get values for use in plot
        ccf_props['mjd'] = Time(np.array(hdict['CCF_MJDMID']))
        dv_vec = hdict['CCF_DV'].to(uu.m / uu.s).value
        ccf_props['dv'] = np.array(dv_vec)
        sdv_vec = hdict['CCF_SDV'].to(uu.m / uu.s).value
        ccf_props['sdv'] = np.array(sdv_vec)
        fwhm_vec = hdict['CCF_FWHM'].to(uu.m / uu.s).value
        ccf_props['fwhm'] = np.array(fwhm_vec)
        ccf_props['masks'] = np.array(hdict['CCF_MASK'])
        ccf_props['files'] = self.filetypes['ccf'].get_files(qc=True)
        ccf_props['files_failed'] = self.filetypes['ccf'].get_files(qc=False)
        # -----------------------------------------------------------------
        ccf_props['FIRST_CCF'] = Time(np.min(hdict['CCF_MJDMID'])).iso
        ccf_props['LAST_CCF'] = Time(np.max(hdict['CCF_MJDMID'])).iso
        ccf_props['LAST_CCF_PROC'] = Time(np.max(hdict['CCF_PROC'])).iso
        # -----------------------------------------------------------------
        # ccf version
        ccf_props['CCF_VERSION'] = ','.join(list(np.unique(hdict['CCF_VERSION'])))
        # -----------------------------------------------------------------
        # standard request keyword args
        rkwargs = dict(fiber='Science fiber',
                       dprtypes=SCIENCE_DPRTYPES,
                       apero_mode=self.settings['CPN'])
        # add the links to request data
        ccf_props['RLINK_CCF'] = self.rlink(filetype='ccf', **rkwargs)
        # -----------------------------------------------------------------
        # select ccf files to use
        ccf_props = choose_ccf_files(ccf_props)
        # load the first file to get the rv vector
        ccf_table0 = load_table(ccf_props['select_files'][0], format='fits',
                                hdu=1)
        # get the rv vector
        ccf_props['rv_vec'] = ccf_table0['RV']
        # storage for the CCF vectors
        all_ccf = np.zeros((len(ccf_props['select_files']), len(ccf_table0)))
        # loop around all other files, load them and load into all_ccf
        for row, select_file in enumerate(ccf_props['select_files']):
            table_row = load_table(select_file, format='fits', hdu=1)
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
        # set the plot title
        plot_title = f'CCF {self.objname} [mask={ccf_props["chosen_mask"]}]'
        # plot the lbl figure
        ccf_plot(ccf_props, plot_path, plot_title=plot_title)
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
        # Create the request link table for this object
        # -----------------------------------------------------------------
        # get the rlink base name
        rlink_base_name = 'ccf_rlink_' + self.objname + '.txt'
        # get the rlink table path
        rlink_item_path = os.path.join(item_save_path, rlink_base_name)
        # define the keys in spec_props that contain rlinks to add
        rlinks = ['RLINK_CCF']
        # define the rlink descriptions
        rlink_text = ['CCF files']
        # compute the rlink table
        create_request_link_table(ccf_props, rlink_item_path, rlinks,
                                  rlink_text)
        # -----------------------------------------------------------------
        # Create the file lists for this object
        # -----------------------------------------------------------------
        # construct the save path for ccf files
        ccf_file = os.path.join(down_save_path,
                                f'ccf_{self.objname}_file_list.txt')
        create_file_list(ccf_files, ccf_file)
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
        self.ccf_rlink_table = item_rel_path + rlink_base_name
        self.ccf_dwn_table = item_rel_path + dwn_base_name

    def get_time_series_parameters(self):
        # get ext files
        ftypes = self.filetypes
        ext_files_all = ftypes['ext'].get_files(qc=True)
        tcorr_files_all = ftypes['tcorr'].get_files(qc=True)
        # don't go here is lbl rdb files are not present
        if len(ext_files_all) == 0:
            return
        # ---------------------------------------------------------------------
        # generate place to save figures
        item_save_path = self.settings['ITEMS']
        item_rel_path = f'../{_ITEM_DIR}/'
        down_save_path = self.settings['DOWNS']
        down_rel_path = f'../{_DOWN_DIR}/'
        # -----------------------------------------------------------------
        # storage for ccf values
        ts_props = dict()
        # get labels
        snr_y_label = self.headers['ext']['EXT_Y']['label']
        snr_y_label = snr_y_label.replace(r'$\mu$', 'u')
        snr_h_label = self.headers['ext']['EXT_H']['label']
        snr_h_label = snr_h_label.replace(r'$\mu$', 'u')
        ext_col = 'ext_files'
        tcorr_col = 'tcorr_files'
        rlink_ext_col = 'Request ext files'
        rlink_tcorr_col = 'Request tcorr files'
        # ---------------------------------------------------------------------
        # construct the stats table
        # ---------------------------------------------------------------------
        # columns
        ts_props['columns'] = TIME_SERIES_COLS[0:9]
        ts_props['columns'] += [snr_y_label, snr_h_label]
        ts_props['columns'] += [TIME_SERIES_COLS[9]]
        ts_props['columns'] += [ext_col, tcorr_col]
        ts_props['columns'] += [rlink_ext_col, rlink_tcorr_col]
        # get values for use in time series table
        for time_series_col in TIME_SERIES_COLS:
            ts_props[time_series_col] = []
        ts_props[snr_y_label] = []
        ts_props[snr_h_label] = []
        ts_props[ext_col] = []
        ts_props[tcorr_col] = []
        ts_props[rlink_ext_col] = []
        ts_props[rlink_tcorr_col] = []
        # get values from self.header_dict
        mjd_vec = np.array(self.header_dict['EXT_MJDMID'])
        seeing_vec = np.array(self.header_dict['EXT_SEEING'])
        airmass_vec = np.array(self.header_dict['EXT_AIRMASS'])
        exptime_vec = np.array(self.header_dict['EXT_EXPTIME'])
        snry_vec = np.array(self.header_dict['EXT_Y'])
        snyh_vec = np.array(self.header_dict['EXT_H'])
        dprtype_vec = np.array(self.header_dict['EXT_DPRTYPE'])
        # get the obs dirs for files that passed qc
        obs_dirs_ext = ftypes['ext'].get_files(qc=True, attr='obsdirs')
        obs_dirs_tcorr = ftypes['tcorr'].get_files(qc=True, attr='obsdirs')
        # get unique object directories (for this object)
        u_obs_dirs = np.unique(obs_dirs_ext)
        # loop around observation directories
        for obs_dir in u_obs_dirs:
            # create a mask for this observation directory
            obs_mask_ext = obs_dirs_ext == obs_dir
            obs_mask_tcorr = obs_dirs_tcorr == obs_dir
            # get the first and last mjd for this observation directory
            first_time = Time(np.min(mjd_vec[obs_mask_ext]))
            last_time = Time(np.max(mjd_vec[obs_mask_ext]))
            first_iso = first_time.iso
            last_iso = last_time.iso
            # get the number of observations for this observation
            num_obs_ext = str(np.sum(obs_mask_ext))
            num_obs_tcorr = str(np.sum(obs_mask_tcorr))
            # get the seeing for this observation directory
            # TODO: nanmean for all below?
            seeing = np.mean(seeing_vec[obs_mask_ext])
            seeing = '{:.3f}'.format(seeing)
            # get the airmass for this observation directory
            airmass = np.mean(airmass_vec[obs_mask_ext])
            airmass = '{:.3f}'.format(airmass)
            # get the mean exposure time
            exptime = np.mean(exptime_vec[obs_mask_ext])
            exptime = '{:.3f}'.format(exptime)
            # get the total exposure time
            texptime = np.sum(exptime_vec[obs_mask_ext])
            texptime = '{:.3f}'.format(texptime)
            # get the mean snr_y
            snry = np.mean(snry_vec[obs_mask_ext])
            snry = '{:.3f}'.format(snry)
            # get the mean snr_h
            snyh = np.mean(snyh_vec[obs_mask_ext])
            snyh = '{:.3f}'.format(snyh)
            # get the dprtypes
            dprtype = ','.join(list(np.unique(dprtype_vec[obs_mask_ext])))
            # -----------------------------------------------------------------
            # Create the ext and tellu for this object
            # -----------------------------------------------------------------
            ext_files = ext_files_all[obs_mask_ext]
            tcorr_files = tcorr_files_all[obs_mask_tcorr]
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
            ts_props[TIME_SERIES_COLS[0]].append(obs_dir)
            ts_props[TIME_SERIES_COLS[1]].append(first_iso)
            ts_props[TIME_SERIES_COLS[2]].append(last_iso)
            ts_props[TIME_SERIES_COLS[3]].append(num_obs_ext)
            ts_props[TIME_SERIES_COLS[4]].append(num_obs_tcorr)
            ts_props[TIME_SERIES_COLS[5]].append(seeing)
            ts_props[TIME_SERIES_COLS[6]].append(airmass)
            ts_props[TIME_SERIES_COLS[7]].append(exptime)
            ts_props[TIME_SERIES_COLS[8]].append(texptime)
            ts_props[snr_y_label].append(snry)
            ts_props[snr_h_label].append(snyh)
            ts_props[TIME_SERIES_COLS[9]].append(dprtype)
            ts_props[ext_col].append(ext_value)
            ts_props[tcorr_col].append(tcorr_value)
            # -----------------------------------------------------------------
            # standard request keyword args
            rkwargs = dict(fiber='Science fiber',
                           dprtypes=SCIENCE_DPRTYPES,
                           apero_mode=self.settings['CPN'])
            # get the date YYYY-MM-DD format
            rlink_start = first_time.strftime('%Y-%m-%d')
            rlink_end = last_time.strftime('%Y-%m-%d')
            # deal with rlink_start and end being the same
            if rlink_end == rlink_start:
                tdelta = TimeDelta(1 * uu.day)
                rlink_end = (last_time + tdelta).strftime('%Y-%m-%d')
            # add the links to request data
            time_series_ext_rlink = self.rlink(filetype='ext',
                                               startdate=rlink_start,
                                               enddate=rlink_end, **rkwargs)
            time_series_tcorr_rlink = self.rlink(filetype='tcorr',
                                                 startdate=rlink_start,
                                                 enddate=rlink_end, **rkwargs)
            # -----------------------------------------------------------------
            # Create the request link table for this object
            # -----------------------------------------------------------------
            # add the ext rlink
            if time_series_ext_rlink is not None:
                rargs = ['Extracted 2D files', time_series_ext_rlink]
                ts_props[rlink_ext_col].append('`{0} <{1}>`_'.format(*rargs))
            # add the corr rlink
            if time_series_tcorr_rlink is not None:
                rargs = ['Telluric corrected 2D files', time_series_tcorr_rlink]
                ts_props[rlink_tcorr_col].append('`{0} <{1}>`_'.format(*rargs))
        # -----------------------------------------------------------------
        # construct the stats
        # -----------------------------------------------------------------
        # get the stats base name
        time_series_base_name = 'time_series_stat_' + self.objname + '.txt'
        # get the stat path
        stat_path = os.path.join(item_save_path, time_series_base_name)
        # compute the stats
        time_series_stats_table(ts_props, stat_path)
        # -----------------------------------------------------------------
        # update the paths
        self.time_series_plot_path = None
        self.time_series_stats_table = item_rel_path + time_series_base_name
        self.time_series_dwn_table = None


# =============================================================================
# Define functions
# =============================================================================
def get_args():
    """
    Define the command line arguments

    :return: argparse namespace object containing arguments
    """
    parser = argparse.ArgumentParser(description='Simple ARI interface')
    # add obs dir
    parser.add_argument('profile', type=str, default='None',
                        help='The profile yaml to use')
    parser.add_argument('--filter', type=str, default='None',
                        help='The specific profile to use in yaml')
    # test mode - do not run apero recipes just test
    parser.add_argument('--debug', type=bool, default=False,
                        help='Run in debug mode (no parallel filters objs)')
    # load arguments with parser
    args = parser.parse_args()
    # return arguments
    return args


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
    param_settings['OBJ_INDEX'] = os.path.join(working_dir, _OBJ_INDEX_DIR)
    param_settings['OBJ_INDEX_ITEM'] = os.path.join(param_settings['OBJ_INDEX'],
                                                    _ITEM_DIR)
    param_settings['OBJ_INDEX_DOWN'] = os.path.join(param_settings['OBJ_INDEX'],
                                                    _DOWN_DIR)
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
    if not os.path.exists(param_settings['OBJ_INDEX']):
        os.makedirs(param_settings['OBJ_INDEX'])
    if not os.path.exists(param_settings['OBJ_INDEX_ITEM']):
        os.makedirs(param_settings['OBJ_INDEX_ITEM'])
    if not os.path.exists(param_settings['OBJ_INDEX_DOWN']):
        os.makedirs(param_settings['OBJ_INDEX_DOWN'])
    # return all parameter settings
    return param_settings


def compile_stats(gsettings: dict, settings: dict, profile: dict,
                  headers: dict) -> dict:
    """
    Compile the stats for a given profile

    :param gsettings: dict, global settings
    :param settings: dict, profile settings
    :param profile: dict, the profile to compile stats for
    :param headers: dict, header keys
    :return: dict, the stats for the profile
    """
    # set up storage for the output dictionary
    profile_stats: Dict[str, Any] = dict()
    profile_stats['TABLES']: Dict[str, Any] = dict()
    profile_stats['OPROPS']: Dict[str, Any] = dict()
    # deal with updating the path (as DRS_UCONFIG) from apero profile
    update_apero_profile(profile)
    # get paths to tables
    object_table_file = os.path.join(settings['DATA'], 'OBJECT_TABLE.fits')
    recipe_table_file = os.path.join(settings['DATA'], 'RECIPE_TABLE.fits')
    # message_table_file = os.path.join(settings['DATA'], 'MESSAGE_TABLE.fits')
    # get skip criteria
    skip_obj_table = profile['skip obj table']
    skip_recipe_table = profile['skip recipe table']
    # skip_msg_table = profile['skip msg table']
    # ------------------------------------------------------------------
    # deal with skipping object table
    if skip_obj_table and os.path.exists(object_table_file):
        profile_stats['TABLES'][TABLE_NAMES[0]] = load_table(object_table_file)
    elif skip_obj_table:
        profile_stats['TABLES'][TABLE_NAMES[0]] = None
    else:
        # get the object table (astropy table)
        object_classes = compile_apero_object_table(gsettings)
        # add the lbl count
        object_classes = add_lbl_count(profile, object_classes)
        # add the object pages
        object_classes = add_obj_pages(gsettings, settings, profile, headers,
                                       object_classes)
        # make object table
        object_table = make_obj_table(object_classes)
        # add final object table to profile stats
        profile_stats['TABLES'][TABLE_NAMES[0]] = object_table
        profile_stats['OPROPS'] = object_classes
    # ------------------------------------------------------------------
    # deal with skipping recipe table
    if skip_recipe_table and os.path.exists(recipe_table_file):
        profile_stats['TABLES'][TABLE_NAMES[1]] = load_table(recipe_table_file)
    elif skip_recipe_table:
        profile_stats['TABLES'][TABLE_NAMES[1]] = None
    else:
        # get the recipe log table
        profile_stats['TABLES'][TABLE_NAMES[1]] = compile_apero_recipe_table()
    # ------------------------------------------------------------------
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
    # ------------------------------------------------------------------
    # cannot import until after DRS_UCONFIG is set
    from apero.base import base
    from apero.core import constants
    from apero.core.constants import param_functions
    # reload DPARAMS and IPARAMS
    base.DPARAMS = base.load_database_yaml()
    base.IPARAMS = base.load_install_yaml()
    # ------------------------------------------------------------------
    # invalidate cache
    param_functions.CONFIG_CACHE = dict()
    # make sure parameters is reloaded (and not cached)
    _ = constants.load(cache=False)


def compile_obj_index_page(gsettings: dict, settings: dict,
                           oprops: Dict[str, Dict[str, ObjectData]]):
    from apero.tools.module.documentation import drs_markdown
    from apero.core import constants
    from apero.core.core import drs_database
    from apero.core.utils import drs_startup
    from apero.base.base import TQDM as tqdm
    # get the parameter dictionary of constants from apero
    params = constants.load()
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # print progress
    # generate place to save figures
    item_save_path = settings['OBJ_INDEX_ITEM']
    item_rel_path = f'{_OBJ_INDEX_DIR}/{_ITEM_DIR}/'
    down_save_path = settings['OBJ_INDEX_DOWN']
    down_rel_path = f'{_OBJ_INDEX_DIR}/{_DOWN_DIR}/'
    # -------------------------------------------------------------------------
    # load object database
    objdbm = drs_database.AstrometricDatabase(params)
    objdbm.load_db()
    # get all objects
    objnames = objdbm.get_entries('OBJNAME')
    # storage for outputs
    objdict = dict()
    # loop around objects and create a section for each
    for objname in tqdm(objnames):
        entry = dict()
        entry['profile_items'] = []
        entry['profile_names'] = []
        entry['find_files'] = []
        entry['find_descs'] = []
        # ---------------------------------------------------------------------
        # add profile reference links for this object
        # ---------------------------------------------------------------------
        for apero_profile_name in oprops:
            # look in this profile
            oprops_profile = oprops[apero_profile_name]
            # deal with no reduction for this profile
            if len(oprops_profile) == 0:
                continue
            # skip missing objects
            if objname not in oprops_profile:
                continue
            # get object class for this objname
            object_class = oprops_profile[objname]
            # add to entry
            entry['profile_items'].append(object_class.objpageref)
            entry['profile_names'].append(apero_profile_name)
        # ---------------------------------------------------------------------
        # find finder charts
        # ---------------------------------------------------------------------
        if gsettings['find directory'] not in [None, 'None', 'Null', '']:
            # get find directory
            find_path = gsettings['find directory']
            # look for objname in this directory
            find_files, find_descs = find_finder_charts(find_path, objname)
            # push to entry
            entry['find_files'] = find_files
            entry['find_descs'] = find_descs

        # ---------------------------------------------------------------------
        # generate finder chart download table
        # ---------------------------------------------------------------------
        if len(entry['find_files']) > 0:
            # make download table and store table to add
            fd_args = [entry, objname, item_save_path, item_rel_path,
                       down_save_path, down_rel_path]
            entry['find_table'] = make_finder_download_table(*fd_args)
        # deal with no finder charts
        else:
            entry['find_table'] = None
        # ---------------------------------------------------------------------
        # add to stroage
        objdict[objname] = entry
    # -------------------------------------------------------------------------
    # create ARI index page
    obj_index_page = drs_markdown.MarkDownPage('object_index')
    # add title
    obj_index_page.add_title('APERO Reduction Interface (ARI) '
                             'Object Index Page')
    # -------------------------------------------------------------------------
    # Add basic text
    # construct text to add
    obj_index_page.add_text('Object Index')
    obj_index_page.add_newline()
    obj_index_page.add_text('Object by object index. '
                            'Links to all profiles and finding charts')
    obj_index_page.add_newline()
    obj_index_page.add_text('Please note: Your object may be under another '
                            'name. Please check `here <https://docs.google.com/'
                            'spreadsheets/d/'
                            '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/'
                            'edit?usp=sharing>`_, the name displayed in '
                            'ARI will be the first column [OBJNAME]')
    obj_index_page.add_newline()
    obj_index_page.add_text('If you have any issues please report using '
                            '`this sheet <https://docs.google.com/spreadsheets/d/1Ea_WEFTlTCbth'
                            'R24aaQm4KaleIteLuXLgn4RiNBnEqs/edit?usp=sharing>`_.')
    obj_index_page.add_newline()
    # -------------------------------------------------------------------------
    # loop around objects and create a section for each
    for objname in objdict:
        # get this iterations entry
        entry = objdict[objname]
        # add reference to this section
        obj_index_page.add_reference(f'indexing_{objname}')
        # add section
        obj_index_page.add_section(objname)
        # add table of contents
        if len(entry['profile_items']) > 0:
            obj_index_page.add_newline()
            obj_index_page.add_table_of_contents(items=entry['profile_items'],
                                                 sectionname=None)
        else:
            obj_index_page.add_newline()
            obj_index_page.add_text('Currently not reduced under any profile')
        # add the finder chart table
        if entry['find_table'] is not None:
            obj_index_page.add_newline()
            # add the finder chart table
            obj_index_page.add_csv_table('', entry['find_table'],
                                         cssclass='csvtable2')
        else:
            obj_index_page.add_newline()
            obj_index_page.add_text('Currently no finder chart')
        obj_index_page.add_newline()
    # -------------------------------------------------------------------------
    # save index page
    obj_index_file = os.path.join(settings['WORKING'], 'obj_index.rst')
    obj_index_page.write_page(obj_index_file)


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
    index_page.add_text('If you have any issues please report using '
                        '`this sheet <https://docs.google.com/spreadsheets/d/1Ea_WEFTlTCbth'
                        'R24aaQm4KaleIteLuXLgn4RiNBnEqs/edit?usp=sharing>`_.')
    index_page.add_newline()
    index_page.add_text('Last updated: {0} [UTC]'.format(Time.now()))
    index_page.add_newline()
    # -------------------------------------------------------------------------
    # add table of contents
    index_page.add_table_of_contents(profile_files)
    # -------------------------------------------------------------------------
    index_page.add_section('Objects index')
    index_page.add_newline()
    index_page.add_text('Object by object index. '
                        'Links to all profiles and finding charts')
    index_page.add_newline()
    index_page.add_text('Please note this includes objects not currently '
                        'observed.')
    index_page.add_newline()
    index_page.add_table_of_contents(items=['obj_index.rst'],
                                     sectionname=None)
    # -------------------------------------------------------------------------
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
        # -----------------------------------------------------------------
        # Add basic text
        # construct text to add
        profile_page.add_text(f'This is the APERO Reduction Interface (ARI) '
                              f'for the reduction: {cprofile_name}')
        profile_page.add_newline()
        profile_page.add_text('Please note: Your object may be under another '
                              'name. Please check `here <https://docs.google.com/'
                              'spreadsheets/d/'
                              '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/'
                              'edit?usp=sharing>`_, the name displayed in ARI will '
                              'be the first column [OBJNAME]')
        profile_page.add_newline()
        profile_page.add_text('If you have any issues please report using '
                              '`this sheet <https://docs.google.com/spreadsheets/d/1Ea_WEFTlTCbth'
                              'R24aaQm4KaleIteLuXLgn4RiNBnEqs/edit?usp=sharing>`_')
        profile_page.add_newline()
        profile_page.add_text('Last updated: {0} [UTC]'.format(Time.now()))
        profile_page.add_newline()
        # -----------------------------------------------------------------
        # store the reference name for profile page table of contents
        table_files = []
        # urls
        table_urls = dict()
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
            title = f'{table_title} ({cprofile_name})'
            table_page.add_title(title)
            # -----------------------------------------------------------------
            # Add basic text
            # construct text to add
            table_page.add_text(f'This is the APERO Reduction Interface (ARI) '
                                f'for the reduction: {cprofile_name}')
            table_page.add_newline()
            table_page.add_text('Please note: Your object may be under another '
                                'name. Please check `here <https://docs.google.com/'
                                'spreadsheets/d/'
                                '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/'
                                'edit?usp=sharing>`_, the name displayed in ARI will '
                                'be the first column [OBJNAME]')
            table_page.add_newline()
            table_page.add_text('If you have any issues please report using '
                                '`this sheet <https://docs.google.com/spreadsheets/d/1Ea_WEFTlTCbth'
                                'R24aaQm4KaleIteLuXLgn4RiNBnEqs/edit?usp=sharing>`_')
            table_page.add_newline()
            table_page.add_text('Last updated: {0} [UTC]'.format(Time.now()))
            table_page.add_newline()

            # save table page
            table_markdown_file = os.path.basename(table_filename).lower()
            table_markdown_file = table_markdown_file.replace('.csv', '.rst')

            # deal with column widths for this file type
            if table is not None and len(table) > 0:
                # if table_name == 'OBJECT_TABLE':
                #     # add the csv version of this table
                #     table_page.add_csv_table('', f'../{_DATA_DIR}/' +
                #                              table_filename,
                #                              cssclass='csvtable2')

                if table_name == 'OBJECT_TABLE':
                    # add the csv version of this table
                    table_page.add_csv_table('', f'../{_DATA_DIR}/' +
                                             table_filename, cssclass='csvtable2')
                    # store the reference name for profile page table of contents
                    #  these are in the same dir as profile_page so do not add the
                    #  rst sub-dir
                    table_files.append(os.path.basename(table_markdown_file))
                elif table_name == 'RECIPE_TABLE':
                    # add the recipe tables
                    add_recipe_tables(settings, table, table_name)
                    # create a URL to link the pages
                    table_url_file = table_filename.replace('.csv', '.html')
                    table_url_file = table_url_file.lower()
                    # add link to a set of links
                    table_urls[title] = f'../tables/{table_url_file}'
            else:
                # if we have no table then add a message
                table_page.add_text('No table created.')
                # store the reference name for profile page table of contents
                #  these are in the same dir as profile_page so do not add the
                #  rst sub-dir
                table_files.append(os.path.basename(table_markdown_file))
            # write table page
            print('Writing table page: {0}'.format(table_markdown_file))
            table_page.write_page(os.path.join(settings['RST'],
                                               table_markdown_file))

        # add table of contents to profile page
        profile_page.add_table_of_contents(table_files)

        # add list of urls
        if len(table_urls) > 0:
            # add a section
            profile_page.add_newline()
            profile_page.add_text('Other: ')
            profile_page.add_newline()
            # add the urls as a list
            for table_url in table_urls:
                # get url from table urls
                url = table_urls[table_url]
                # add the url
                profile_page.lines += [f'* `{table_url} <{url}>`_']
                # add a new line
                profile_page.add_newline()
            profile_page.add_newline()
        # save profile page
        profile_page.write_page(os.path.join(settings['RST'], 'profile.rst'))


def recipe_date_table(table: Table, table_name: str,
                      ) -> Tuple[Table, List[str], List[str], np.ndarray]:
    """
    Create the date table which links to each date page

    :param table:
    :param table_name:
    :return:
    """
    # define the columns (passed back to main code)
    date_colnames = ['DATE', 'NUM_TOTAL', 'NUM_FAIL', 'LINK']
    date_coltypes = ['str', 'str', 'str', 'url']
    table_dates = []
    # dictionary for storage
    date_dict = dict()
    # add columns
    for col in date_colnames:
        date_dict[col] = []
    # get table name
    table_filename = table_name.lower()
    # loop around rows in table
    for row in range(len(table)):
        # convert start time into a YYYY-MM-DD
        date = table['START_TIME'][row].split(' ')[0]
        # get the date given to this row
        table_dates.append(date)
    # make the table dates a numpy array
    table_dates = np.array(table_dates)
    # get a list of unique dates
    unique_dates = list(set(table_dates))

    for unique_date in unique_dates:
        # get a mask for this date in the original table
        mask = table_dates == unique_date
        # count the number of entries
        num_entries = np.sum(mask)
        # count the number of errors (False or 0)
        num_errors = np.sum(table['ENDED'][mask] == 'False')
        num_errors += np.sum(table['ENDED'][mask] == 0)
        # create html url link
        html_file = f'{table_filename}_{unique_date.replace("-", "_")}.html'
        # deal with populating date table
        date_dict['DATE'].append(unique_date)
        date_dict['NUM_TOTAL'].append(num_entries)
        date_dict['NUM_FAIL'].append(num_errors)
        date_dict['LINK'].append(html_file)
    # convert date_dict into table
    date_table = Table(date_dict)
    # sort by date (newest first)
    sortmask = np.argsort(date_dict['DATE'])[::-1]
    date_table = date_table[sortmask]
    # return the date_table, colnames, coltypes and the date value for each col
    return date_table, date_colnames, date_coltypes, table_dates


def add_recipe_tables(settings: Dict[str, Any], table: Table, table_name: str):
    # import from apero
    from apero.tools.module.error import error_html
    # set html body
    # Take directly from one of the sphinx pages (this is a massive hack)
    html_body1 = """
    <div class="pageheader">

    <ul>
    <li><a title="Home" href="http://apero.exoplanets.ca">
        <i class="fa fa-home fa-3x" aria-hidden="true"></i></a></li>
    <li><a title="install" href="http://apero.exoplanets.ca/user/general/installation">
        <i class="fa fa-cog fa-3x" aria-hidden="true"></i></a></li>
    <li><a title="github" href="https://github.com/njcuk9999/apero-drs">
        <i class="fa fa-git-square fa-3x" aria-hidden="true"></i></a></li>
    <li><a title="download paper" href="https://ui.adsabs.harvard.edu/abs/2022PASP..134k4509C">
        <i class="fa fa-file-pdf-o fa-3x" aria-hidden="true"></i></a></li>
    <li><a title="UdeM" href="http://apero.exoplanets.ca/main/misc/udem.html">
        <i class="fa fa-university fa-3x" aria-hidden="true"></i></a></li>
  </ul>

    <div>
    <a href="http://apero.exoplanets.ca">
      <img src="../../_static/images/apero_logo.png" alt="APERO" />
    </a>
    <br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A PipelinE to Reduce Observations
    </div>

    </div>
    
    <div class="document">
      <div class="documentwrapper">
        <div class="bodywrapper">
          <div class="body" role="main">
    
    <h1>{TITLE}</h1>
    <p> Note the date is the date processed NOT the observation directory 
        (or night directory)</p>
    <br>
    <p><a href="../rst/profile.html">Back to profile page ({PROFILE})</a></p>
    <br>
    <p> 
    A list of known errors can be found 
    <a href="https://docs.google.com/spreadsheets/d/15Gu_aY6h9Esw1uTF8Y5JCHl6m7191AviJNTPkbeTiQE/edit?usp=sharing">here</a>
    <br>
    Please report any errors missing.
    </p>
    <br>
    
    
    """

    html_body2 = """
            <div class="clearer"></div>
          </div>
        </div>
      </div>
      <div class="clearer"></div>
    </div>
    """
    # set html table class
    table_class = 'class="csvtable2 docutils align-default"'
    # css to include
    ccs_files = ['../../_static/pygments.css',
                 '../../_static/bizstyle.css',
                 '../../_static/apero.css']
    # profile name
    cprofile_name = settings['CPN']
    # make path
    table_path = os.path.join(settings['WORKING'], 'extras',
                              cprofile_name, 'tables')
    # get table name
    table_filename = table_name.lower()
    # split the table into sub-tables based on start date
    dout = recipe_date_table(table, table_name)
    date_table, date_colnames, date_coltypes, table_dates = dout
    # loop around sub tables
    for date in date_table['DATE']:
        # find all rows in table that conform to this date
        tablemask = table_dates == date
        # make the table
        subtable = table[tablemask]
        # get the filename to create
        subtable_filename = f'{table_filename}_{date.replace("-", "_")}'
        # get html col names
        html_out_col_names = HTML_OUTCOL_NAMES[table_name]
        html_col_types = HTML_COL_TYPES[table_name]
        # convert table to outlist
        tout = error_html.table_to_outlist(subtable, html_out_col_names,
                                           out_types=html_col_types)
        outlist, t_colnames, t_coltype = tout

        # convert outlist to a html/javascript table
        html_table = error_html.filtered_html_table(outlist, t_colnames,
                                                    t_coltype, clean=False,
                                                    log=False,
                                                    table_class=table_class)
        # deal with path not existing
        if not os.path.exists(table_path):
            os.makedirs(table_path)
        # construct local path to save html to
        subtable_html = os.path.join(table_path, f'{subtable_filename}.html')
        # build html page

        html_title = 'Recipe log for {0} ({1})'.format(date, cprofile_name)
        html_body1_filled = html_body1.format(TITLE=html_title,
                                              PROFILE=cprofile_name)

        html_body1_filled += f"""
        <p>
        <a href="recipe_table.html">Back to recipe log ({cprofile_name})</a>
        </p>
        <br>
        """

        html_content = error_html.full_page_html(html_body1=html_body1_filled,
                                                 html_table=html_table,
                                                 html_body2=html_body2,
                                                 css=ccs_files)
        # write html page
        with open(subtable_html, 'w') as wfile:
            wfile.write(html_content)

    # -------------------------------------------------------------------------
    # make recipe table
    # -------------------------------------------------------------------------
    # construct local path to save html to
    table_html = os.path.join(table_path, f'{table_filename.lower()}.html')
    # convert table to outlist
    tout = error_html.table_to_outlist(date_table, date_colnames,
                                       out_types=date_coltypes)
    outlist, t_colnames, t_coltype = tout
    # convert outlist to a html/javascript table
    html_table = error_html.filtered_html_table(outlist,
                                                t_colnames,
                                                t_coltype,
                                                clean=False, log=False,
                                                table_class=table_class)
    # build html page
    html_title = 'Recipe log ({0})'.format(cprofile_name)
    html_body1_filled = html_body1.format(TITLE=html_title,
                                          PROFILE=cprofile_name)
    html_content = error_html.full_page_html(html_body1=html_body1_filled,
                                             html_table=html_table,
                                             html_body2=html_body2,
                                             css=ccs_files)
    # write html page
    with open(table_html, 'w') as wfile:
        wfile.write(html_content)


def compile_docs(gsettings: dict, settings: dict):
    """
    Compile the documentation

    :param gsettings: dict, the global settings
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
    # get reset
    reset = gsettings['reset']
    # ------------------------------------------------------------------
    # get list of files to copy
    copy_files = ['conf.py', 'make.bat', 'Makefile']
    # loop around files and copy
    for copy_file in copy_files:
        outpath = os.path.join(settings['WORKING'], copy_file)
        # remove file if it exists
        if os.path.exists(outpath):
            os.remove(outpath)
        # copy conf.py make.bat and Makefile to the working directory
        shutil.copy(__file__.replace('simple_ari.py', copy_file), outpath)
    # get _static directory
    static_outdir = os.path.join(settings['WORKING'], '_static')
    # deal with static_outdir existing
    if os.path.exists(static_outdir) and reset:
        shutil.rmtree(static_outdir)
    # copy the static directory as well
    shutil.copytree(__file__.replace('simple_ari.py', '_static'),
                    os.path.join(settings['WORKING'], '_static'),
                    dirs_exist_ok=True)
    shutil.copytree(__file__.replace('simple_ari.py', '_templates'),
                    os.path.join(settings['WORKING'], '_templates'),
                    dirs_exist_ok=True)
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
    # copy extras
    # ------------------------------------------------------------------
    wlog(params, '', 'Copying extra files')
    extras_dir = os.path.join(settings['WORKING'], 'extras')
    # copy everything from extras into output directory
    shutil.copytree(extras_dir, os.path.join(settings['OUT']),
                    dirs_exist_ok=True)
    # ------------------------------------------------------------------
    # change back to current directory
    os.chdir(cwd)


def sync_docs(gsettings: dict, settings: dict):
    # must import here (so that os.environ is set)
    # noinspection PyPep8Naming
    from apero.core import constants
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
    # make sure we copy contents not directory
    if not out_dir.endswith(os.sep):
        out_dir += os.sep
    # get rsync dict
    rdict = dict()
    rdict['SSH'] = gsettings['ssh']['options']
    rdict['USER'] = gsettings['ssh']['user']
    rdict['HOST'] = gsettings['ssh']['host']
    rdict['INPATH'] = os.path.join(gsettings['ssh']['directory'],
                                   f'ari/{instrument}/')
    rdict['OUTPATH'] = out_dir
    # print command to rsync
    wlog(params, '', RSYNC_CMD_IN.format(**rdict))
    # run command (will require password)
    os.system(RSYNC_CMD_IN.format(**rdict))
    # ------------------------------------------------------------------
    # rsync compiled files
    rdict = dict()
    rdict['SSH'] = gsettings['ssh']['options']
    rdict['USER'] = gsettings['ssh']['user']
    rdict['HOST'] = gsettings['ssh']['host']
    rdict['INPATH'] = os.path.join(gsettings['ssh']['directory'],
                                   f'ari/profile/')
    rdict['OUTPATH'] = settings['WORKING'] + os.sep
    # print command to rsync
    wlog(params, '', RSYNC_CMD_IN.format(**rdict))
    # run command (will require password)
    os.system(RSYNC_CMD_IN.format(**rdict))


def upload_docs(gsettings: dict, settings: dict):
    """
    Upload the documentation to the web server

    :param gsettings: dict, global settings
    :param settings: dict, the settings for the documentation

    :return:
    """
    # must import here (so that os.environ is set)
    # noinspection PyPep8Naming
    from apero.core import constants
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
    wlog(params, '', RSYNC_CMD_OUT.format(**rdict))
    # run command (will require password)
    os.system(RSYNC_CMD_OUT.format(**rdict))
    # ------------------------------------------------------------------
    # rsync compiled files
    rdict = dict()
    rdict['SSH'] = gsettings['ssh']['options']
    rdict['USER'] = gsettings['ssh']['user']
    rdict['HOST'] = gsettings['ssh']['host']
    rdict['INPATH'] = settings['WORKING'] + os.sep
    rdict['OUTPATH'] = os.path.join(gsettings['ssh']['directory'],
                                    f'ari/profile/')
    # print command to rsync
    wlog(params, '', RSYNC_CMD_OUT.format(**rdict))
    # run command (will require password)
    os.system(RSYNC_CMD_OUT.format(**rdict))


# =============================================================================
# Functions for the apero stats
# =============================================================================
def compile_apero_object_table(gsettings) -> Dict[str, ObjectData]:
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
    # -------------------------------------------------------------------------
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
    # -------------------------------------------------------------------------
    # log progress
    wlog(params, '', 'Loading objects from astrometric database')
    # -------------------------------------------------------------------------
    # get the astrometric database from apero
    astrodbm = drs_database.AstrometricDatabase(params)
    astrodbm.load_db()
    # -------------------------------------------------------------------------
    # deal with filtering by object
    if gsettings['filter objects']:
        subconditions = []
        for objname in gsettings['objects']:
            subconditions.append(f'OBJNAME="{objname}"')
        condition = ' OR '.join(subconditions)
        condition = f'({condition})'
    else:
        condition = None
    # -------------------------------------------------------------------------
    # define types of files we want to count
    filetypes = dict()
    filetypes['raw'] = FileType('raw', block_kind='raw')
    filetypes['pp'] = FileType('pp', block_kind='tmp', chain='raw')
    filetypes['ext'] = FileType('ext', block_kind='red', chain='pp',
                                kw_output='EXT_E2DS_FF', fiber=science_fiber)
    filetypes['tcorr'] = FileType('tcorr', block_kind='red', chain='ext',
                                  kw_output='TELLU_OBJ',
                                  fiber=science_fiber)
    filetypes['ccf'] = FileType('ccf', block_kind='red', chain='tcorr',
                                kw_output='CCF_RV', fiber=science_fiber)
    filetypes['polar'] = FileType('polar', block_kind='red', chain='tcorr',
                                  kw_output='POL_DEG',
                                  fiber=science_fiber, count=has_polar)
    filetypes['efiles'] = FileType('efiles', block_kind='out', chain='pp',
                                   kw_output='DRS_POST_E')
    filetypes['tfiles'] = FileType('tfiles', block_kind='out', chain='ext',
                                   kw_output='DRS_POST_T')
    filetypes['vfiles'] = FileType('vfiles', block_kind='out', chain='tcorr',
                                   kw_output='DRS_POST_V')
    filetypes['pfiles'] = FileType('pfiles', block_kind='out', chain='tcorr',
                                   kw_output='DRS_POST_P', count=has_polar)
    filetypes['s1d'] = FileType('s1d', block_kind='red', chain='ext',
                                kw_output='EXT_S1D_V', fiber=science_fiber)
    filetypes['sc1d'] = FileType('sc1d', block_kind='red', chain='tcorr',
                                 kw_output='SC1D_V_FILE', fiber=science_fiber)
    # lbl files added as filetype but don't count in same was as other files
    filetypes['lbl.fits'] = FileType('lbl.fits', count=False)
    for filetype in LBL_FILETYPES:
        filetypes[filetype] = FileType(filetype, count=False)
    # -------------------------------------------------------------------------
    # log that we are loading
    # get the object table from the astrometric database
    object_table = astrodbm.get_entries(columns=','.join(ASTROMETRIC_COLUMNS),
                                        condition=condition)
    # -------------------------------------------------------------------------
    # create objects
    obj_classes = dict()
    for row, objname in enumerate(object_table[OBJECT_COLUMN]):
        obj_class = ObjectData(objname, filetypes)
        # add astrometric data
        obj_class.add_astrometrics(object_table.iloc[row])
        # set the has_polar key
        obj_class.has_polar = has_polar
        # append to storage
        obj_classes[objname] = obj_class
    # -------------------------------------------------------------------------
    # get the index database from file index database
    indexdbm = drs_database.FileIndexDatabase(params)
    indexdbm.load_db()
    # get the log database
    logdbm = drs_database.LogDatabase(params)
    logdbm.load_db()
    # -------------------------------------------------------------------------
    # log progress
    wlog(params, '', 'Compiling object files stats (this may take a while)')
    # -------------------------------------------------------------------------
    # for each object we run several counts
    # loop around objects in the object table
    for objname in tqdm(list(obj_classes.keys())):
        # get object class
        obj_class = obj_classes[objname]
        # add files stats
        obj_class.add_files_stats(indexdbm, logdbm)
    # -------------------------------------------------------------------------
    # sort objects by name
    obj_classes_sorted = dict()
    # loop through all objects sorted alphabetically
    for objname in np.sort(list(obj_classes.keys())):
        # get the object class for this object name
        obj_class = obj_classes[objname]
        # reject objects that have no raw files
        if obj_class.filetypes['raw'].num == 0:
            continue
        # add to sorted list of objects
        obj_classes_sorted[objname] = obj_class
    # return object table
    return obj_classes_sorted


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
    # sort columns based on LOG_COLUMNS (not sql order) and as a table
    out_log_table = Table()
    for c_it, col in enumerate(LOG_COLUMNS):
        # deal with bools
        if LOG_TYPES[c_it] == 'bool':
            # find null values
            mask = (log_table[col] == 1) | (log_table[col] == 0)
            # convert log_table column to bools replacing
            #   null values with 0 and then convert to string and fill
            #   null values with blanks
            vector = log_table[col]
            vector[~mask] = 0
            vector = np.array(vector).astype(bool)
            vector = np.array(vector).astype(str)
            vector[~mask] = ''
            out_log_table[col] = vector
        else:
            out_log_table[col] = np.array(log_table[col]).astype(str)

    # ------------------------------------------------------------------
    return out_log_table


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
                # noinspection PyBroadException
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
def add_lbl_count(profile: dict, object_classes: Dict[str, ObjectData]
                  ) -> Dict[str, ObjectData]:
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
    if lbl_path in params:
        lbl_path = params[lbl_path]
    # -------------------------------------------------------------------------
    # deal with no valid lbl path
    if lbl_path is None:
        return object_classes
    # deal with lbl path not existing
    if not os.path.exists(lbl_path):
        return object_classes
    # print that we are analysing lbl outputs
    wlog(params, '', 'Analysing LBL files')
    # -------------------------------------------------------------------------
    # loop around objects
    for objname in tqdm(object_classes):
        # get the object class for this objname
        object_class = object_classes[objname]
        # ---------------------------------------------------------------------
        # LBL RV files
        # ---------------------------------------------------------------------
        # for each object find all directories in lbl path that match this
        #   object name
        lblrv_dir = glob.glob(os.path.join(lbl_path, 'lblrv', f'{objname}_*'))
        # ---------------------------------------------------------------------
        # deal with no directories --> skip
        if len(lblrv_dir) == 0:
            continue
        # ---------------------------------------------------------------------
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
        _select = templates[select]
        _count = int(counts[select])
        # add the counts to the object class
        object_class.lbl_templates = templates
        object_class.lbl_select = _select
        # set the number of files
        object_class.filetypes['lbl.fits'].num = _count
        # ---------------------------------------------------------------------
        # LBL files
        # ---------------------------------------------------------------------
        # loop around all lbl files
        for it, filetype in enumerate(LBL_FILETYPES):
            lbl_files = []
            # loop around templates
            for template in templates:
                # push the object and template name into the glob
                lbl_file = LBL_FILENAMES[it].format(objname, template)
                # get all the lbl files for this object name
                lbl_file_path = os.path.join(lbl_path, 'lblrdb', lbl_file)
                # remove drift files from the lbl rdb files
                if os.path.exists(lbl_file_path):
                    lbl_files.append(lbl_file_path)
            # add list to the LBLRDB file dict for this object
            object_class.filetypes[filetype].files = lbl_files
            # add to object table
            object_class.filetypes[filetype].num = len(lbl_files)
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
            object_class.lbl_stat_files[lblfilekey] = dict()
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
                    object_class.lbl_stat_files[lblfilekey][objtmp] = lblspath
    # -------------------------------------------------------------------------
    # return the object table
    return object_classes


def choose_ccf_files(ccf_props: Dict[str, Any]) -> Dict[str, Any]:
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

    ccf_props['select_files'] = sfiles
    ccf_props['chosen_mask'] = chosen_mask
    # return ccf props
    return ccf_props


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
    # noinspection PyBroadException
    try:
        # noinspection PyTupleAssignmentBalance
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
def add_obj_page(it: int, key: str, profile: dict, gsettings: dict,
                 settings: dict, headers,
                 object_classes: Dict[str, ObjectData],
                 return_dict: Any = None) -> Dict[str, Any]:
    # get object
    object_class = object_classes[key]
    # deal with no return_dict
    if return_dict is None:
        return_dict = dict()
    # add settings to object class
    object_class.add_settings(profile, settings, gsettings, headers)
    # get the object name for this row
    objname = object_class.objname
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
    msg = '\tCreating page for {0} [{1} of {2}]'
    margs = [objname, it + 1, len(object_classes)]
    wlog(params, '', msg.format(*margs))
    # ---------------------------------------------------------------------
    # populate the header dictionary for this object instance
    # wlog(params, '', f'\t\tPopulating header dictionary')
    object_class.populate_header_dict()
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
    object_page.add_newline(nlines=3)
    # ---------------------------------------------------------------------
    # Spectrum section
    # ---------------------------------------------------------------------
    # print progress
    # wlog(params, '', f'\t\tCreating spectrum section')
    # add spectrum section
    objpage_spectrum(object_page, names[0], items[0], object_class)
    # ---------------------------------------------------------------------
    # LBL section
    # ---------------------------------------------------------------------
    # print progress
    # wlog(params, '', f'\t\tCreating LBL section')
    # add LBL section
    objpage_lbl(object_page, names[1], items[1], object_class)
    # ---------------------------------------------------------------------
    # CCF section
    # ---------------------------------------------------------------------
    # print progress
    # wlog(params, '', f'\t\tCreating CCF section')
    # add CCF section
    objpage_ccf(object_page, names[2], items[2], object_class)
    # ---------------------------------------------------------------------
    # Time series section
    # ---------------------------------------------------------------------
    # print progress
    # wlog(params, '', f'\t\tCreating time series section')
    # add time series section
    objpage_timeseries(object_page, names[3], items[3], object_class)
    # ---------------------------------------------------------------------
    # print progress
    # wlog(params, '', f'\t\tWriting to disk')
    # construct a path for the object name
    object_page_path = settings['OBJ_OUT']
    # construct the rst filename
    rst_filename = f'{objname}.rst'
    # save object page
    object_page.write_page(os.path.join(object_page_path, rst_filename))
    # get ref path
    obj_ref_page = os.path.join('.', clean_name, _OBJ_OUT_DIR, rst_filename)
    # ---------------------------------------------------------------------
    # return dictioanry of properties
    rprops = dict()
    rprops['OBJNAME'] = str(objname)
    rprops['OBJURL'] = str(obj_url)
    rprops['OBJPAGEREF'] = str(obj_ref_page)
    return_dict[it] = rprops
    # ------------------------------------------------------------------
    # print progress
    msg = '\tFinished creating page for {0} [{1} of {2}]'
    margs = [objname, it + 1, len(object_classes)]
    wlog(params, '', msg.format(*margs), colour='magenta')
    # things to return
    return return_dict


def add_obj_pages(gsettings: dict, settings: dict, profile: dict,
                  headers: dict, object_classes: Dict[str, ObjectData]
                  ) -> Dict[str, ObjectData]:
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
    # -------------------------------------------------------------------------
    # deal with no entries in object table
    if len(object_classes) == 0:
        # print progress
        wlog(params, '', 'No objects found in object table')
        # return empty table
        return object_classes
    # -------------------------------------------------------------------------
    # print progress
    wlog(params, 'info', 'Creating object pages')
    # set up the arguments for the multiprocessing
    args = [0, '', profile, gsettings, settings, headers, object_classes]
    # get the number of cores
    n_cores = gsettings.get('N_CORES', profile.get('N_CORES', None))
    if n_cores is None:
        raise ValueError('Must define N_CORES in settings or profile')
    # storage for results
    results_dict = dict()
    # -------------------------------------------------------------------------
    # deal with running on a single core
    if n_cores == 1:
        # change the object column to a url
        for it, key in enumerate(object_classes):
            # combine arguments
            itargs = [it, key] + args[2:]
            # run the pool
            results = add_obj_page(*itargs)
            # push result to result storage
            results_dict[key] = results[it]
    # -------------------------------------------------------------------------
    elif n_cores > 1:
        if MULTI == 'POOL':
            from multiprocessing import get_context
            # list of params for each entry
            params_per_process = []
            for it, key in enumerate(object_classes):
                itargs = [it, key] + args[2:]
                params_per_process.append(itargs)
            # start parellel jobs
            with get_context('spawn').Pool(n_cores, maxtasksperchild=1) as pool:
                results = pool.starmap(add_obj_page, params_per_process)
            # fudge back into return dictionary
            for row in range(len(results)):
                objname = results[row]['OBJNAME']
                # push into results dict
                results_dict[objname] = results[row]
        else:
            from multiprocessing import Process, Manager
            # split up groups
            group_iterations, group_keys = [], []
            all_iterations = list(range(len(object_classes)))
            all_keys = list(object_classes.keys())
            ngroups = int(np.ceil(len(object_classes)/n_cores))
            for group_it in range(ngroups):
                start = group_it * n_cores
                end = (group_it * n_cores) + n_cores
                iterations = all_iterations[start:end]
                keys = all_keys[start:end]
                # push into storage
                group_iterations.append(iterations)
                group_keys.append(keys)
            # start the process manager
            manager = Manager()
            return_dict = manager.dict()
            # do the multiprocessing
            for group_it in range(ngroups):
                jobs = []
                # loop around jobs in group
                for group_jt in range(len(group_iterations[group_it])):
                    group_args = [group_iterations[group_it][group_jt],
                                  group_keys[group_it][group_jt]] + args[2:]
                    group_args += [return_dict]
                    # get parallel process
                    process = Process(target=add_obj_page, args=group_args)
                    process.start()
                    jobs.append(process)
                # do not continue until finished
                for pit, proc in enumerate(jobs):
                    proc.join()
                # fudge back into return dictionary
                for row in return_dict.keys():
                    objname = str(return_dict[row]['OBJNAME'])
                    results_dict[objname] = dict(return_dict[row])
    # -------------------------------------------------------------------------
    # update object classes with results
    # -------------------------------------------------------------------------
    # replace object name with the object name + object url
    for key in object_classes:
        if key in results_dict:
            # get the object class for this key
            object_class = object_classes[key]
            # -----------------------------------------------------------------
            # update results
            # -----------------------------------------------------------------
            # This is where we add any results coming back from add_obj_page
            object_class.objurl = results_dict[key]['OBJURL']
            object_class.objpageref = results_dict[key]['OBJPAGEREF']
    # -------------------------------------------------------------------------
    # return the object table
    return object_classes


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
    if len(object_instance.filetypes['ext'].files) == 0:
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
    # add request table
    if object_instance.spec_rlink_table is not None:
        # add the stats table
        page.add_csv_table('', object_instance.spec_rlink_table,
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
    # get the first lbl files
    lbl_files = dict()
    for filetype in LBL_FILETYPES:
        if object_instance.filetypes[filetype].num > 0:
            lbl_files[filetype] = object_instance.filetypes[filetype].files
    # ------------------------------------------------------------------
    # deal with no spectrum found
    if len(lbl_files) == 0:
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
        # add request table
        if object_instance.lbl_rlink_table[objcomb] is not None:
            # add the stats table
            page.add_csv_table('', object_instance.lbl_rlink_table[objcomb],
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
    # get lbl rdb files
    ccf_files = object_instance.filetypes['ccf'].files
    # ------------------------------------------------------------------
    # deal with no spectrum found
    if len(ccf_files) == 0:
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
    # add request table
    if object_instance.ccf_rlink_table is not None:
        # add the stats table
        page.add_csv_table('', object_instance.ccf_rlink_table,
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
    if object_instance.time_series_dwn_table is not None:
        # add the stats table
        page.add_csv_table('', object_instance.time_series_dwn_table,
                           cssclass='csvtable2')


def make_obj_table(object_instances: Dict[str, ObjectData]) -> Optional[Table]:
    # storage dictionary for conversion to table
    table_dict = dict()
    # deal with no entries
    if len(object_instances) == 0:
        return None
    # get the first instance (for has polar)
    key = list(object_instances.keys())[0]
    object_class0 = object_instances[key]
    # start columns as empty lists
    for col in COLUMNS:
        if col in POLAR_COLS:
            if object_class0.has_polar:
                table_dict[col] = []
        else:
            table_dict[col] = []
    # -------------------------------------------------------------------------
    # loop around rows in table
    for key in object_instances:
        # get the class for this objname
        object_class = object_instances[key]
        # set the object name
        table_dict['OBJNAME'].append(object_class.objurl)
        # set the ra and dec
        table_dict['RA'].append(object_class.ra)
        table_dict['DEC'].append(object_class.dec)
        # set the Teff
        if object_class.teff is None:
            table_dict['TEFF'].append(np.nan)
        else:
            table_dict['TEFF'].append(object_class.teff)
        # set the SpT
        if object_class.sptype is None:
            table_dict['SPTYPE'].append('')
        else:
            table_dict['SPTYPE'].append(object_class.sptype)
        # set the dprtypes
        table_dict['DPRTYPE'].append(object_class.dprtypes)
        # set the raw number of files
        table_dict['RAW'].append(object_class.filetypes['raw'].num)
        # set the number of pp files (passed qc)
        table_dict['PP'].append(object_class.filetypes['pp'].num_passed)
        # set the number of ext files (passed qc)
        table_dict['EXT'].append(object_class.filetypes['ext'].num_passed)
        # set the number of tcorr files (passed qc)
        table_dict['TCORR'].append(object_class.filetypes['tcorr'].num_passed)
        # set the number of ccf files (passed qc)
        table_dict['CCF'].append(object_class.filetypes['ccf'].num_passed)
        # set the number of polar files (passed qc)
        if object_class.has_polar:
            table_dict['POL'].append(object_class.filetypes['polar'].num_passed)
        # set the number of e.fits files
        table_dict['efits'].append(object_class.filetypes['efiles'].num)
        # set the number of t.fits files
        table_dict['tfits'].append(object_class.filetypes['tfiles'].num)
        # set the number of v.fits files
        table_dict['vfits'].append(object_class.filetypes['vfiles'].num)
        # set the number of p.fits files
        if object_class.has_polar:
            table_dict['pfits'].append(object_class.filetypes['pfiles'].num)
        # set the number of lbl files
        table_dict['lbl'].append(object_class.filetypes['lbl.fits'].num)
        # set the last observed value raw file
        table_dict['last_obs'].append(object_class.filetypes['raw'].last.iso)
        # set the last processed value
        if object_class.last_processed is not None:
            table_dict['last_proc'].append(object_class.last_processed.iso)
        else:
            table_dict['last_proc'].append('')
    # -------------------------------------------------------------------------
    # finally convert this to a table but use the output column names
    out_table = Table()
    # loop around all columns
    for col in list(table_dict.keys()):
        # get the out table column name
        new_col = COLUMNS[col]
        # push values into out table
        out_table[new_col] = table_dict[col]
    # -------------------------------------------------------------------------
    # sort the column by the object name
    out_table.sort('OBJNAME')
    # return the table
    return out_table


# =============================================================================
# Plots, stat tables and download tables
# =============================================================================
def download_table(files: List[str], descriptions: List[str],
                   item_path: str, dwn_rel_path: str, down_dir: str,
                   title: str, versions: Optional[List[str]] = None):
    """
    Generic download table saving to the item relative path

    :param files: the list of files to add to download table
    :param descriptions: the list of descriptions for each file
    :param dwn_rel_path: the path of the download relative to the page
    :param item_path: the absolute path to the item csv table file (to save to)
    :param down_dir: the path to the download directory
    :param title: the title for the download table
    :param versions: the versions of the files (optional)
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
    # storage for the version
    out_versions = dict()
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
        # ---------------------------------------------------------------------
        # get the last modified time of the filename
        if os.path.exists(filename):
            last_mod = Time(os.path.getmtime(filename), format='unix').iso
            last_modified[basename] = last_mod
        else:
            last_modified[basename] = 'N/A'
        # ---------------------------------------------------------------------
        # deal with version
        if versions is not None:
            out_versions[basename] = versions[it]
        # if file is a fits file get the version from the header
        elif filename.endswith('.fits'):
            # get the header
            hdr = fits.getheader(filename)
            if 'VERSION' in hdr:
                out_versions[basename] = hdr['VERSION']
            else:
                out_versions[basename] = ''
        # otherwise we have no version info
        else:
            out_versions[basename] = ''

    # --------------------------------------------------------------------------
    # construct the stats table
    # --------------------------------------------------------------------------
    # start with a download dictionary
    down_dict = dict(Description=[], Value=[], Uploaded=[], Version=[])
    # flag for having versions
    has_version = False
    # loop around files
    for basename in in_paths:
        # add the rdb file
        down_dict['Description'].append(descs[basename])
        down_dict['Value'].append(_make_download(basename, ref_paths[basename]))
        down_dict['Uploaded'].append(last_modified[basename])
        down_dict['Version'].append(out_versions[basename])
        # check for version
        if out_versions[basename] != '':
            has_version = True
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
    # add version only if one file has version (with the has_version flag)
    if has_version:
        down_dict2['Version'] = down_dict['Version']
    # --------------------------------------------------------------------------
    # convert to table
    down_table = Table(down_dict2)
    # write to file as csv file
    down_table.write(item_path, format='ascii.csv', overwrite=True)


def create_request_link_table(props: Dict[str, Any], save_path: str,
                              rlinks: List[str], rlinkstxt: List[str]):
    # rlinks
    rlink_dict = dict()
    rlink_dict['Request'] = []
    # loop around rlinks
    for r_it, rlink in enumerate(rlinks):
        if rlink not in props:
            continue
        if props[rlink] is None:
            continue
        # get the rlink text
        rlinktxt = rlinkstxt[r_it]
        # construct the request link in markdown format
        rlinkstr = '`{0} <{1}>`_'.format(rlinktxt, props[rlink])
        # push into dictionary
        rlink_dict['Request'].append(rlinkstr)
    # --------------------------------------------------------------------------
    # convert to table
    rdict_table = Table(rlink_dict)
    # write to file as csv file
    rdict_table.write(save_path, format='ascii.csv', overwrite=True)


def create_file_list(files: List[str], path: str):
    """
    Writes a list of files to disk
    """
    # if file exists remove it
    if os.path.exists(path):
        # noinspection PyBroadException
        try:
            os.remove(path)
        except Exception as _:
            pass
    # sort files alphabetically
    files = np.sort(files)
    # open file
    with open(path, 'w') as filelist:
        # loop around files
        for filename in files:
            # write to file
            filelist.write(filename + '\n')


def create_header_file(files: List[str], headers: Dict[str, Any], filetype: str,
                       hdict: Dict[str, Any], filename: str):
    """
    Creates a header file (csv) from a dictionary of header keys

    :param files: list of str, the files to loop around
    :param headers: dict, a dictionary containing those headers to add for
                    each filetype
    :param filetype: str, a filetype (i.e. ext, pp, tcorr, ccf)
    :param hdict: dict, the header values
    :param filename: str, the filename to save the csv file to

    :return: None, writes file to disk
    """
    # check file kind in headers (if it isn't we don't create files)
    if filetype not in headers:
        return
    # storage for dict-->table
    tabledict = dict()
    # first column is the filenames
    tabledict['filename'] = [os.path.basename(filename) for filename in files]
    # loop around keys in header for this filetype and add to tabledict
    for keydict in headers[filetype]:
        tabledict[keydict] = hdict[keydict]
    # convert to table
    table = Table(tabledict)
    # write to file
    table.write(filename, format='ascii.csv', overwrite=True)


def spec_plot(spec_props: Dict[str, Any], plot_path: str, plot_title: str):
    # get parameters from props
    mjd = spec_props['mjd']
    ext_y = spec_props['EXT_Y']
    ext_h = spec_props['EXT_H']
    ext_y_label = spec_props['EXT_Y_LABEL']
    ext_h_label = spec_props['EXT_H_LABEL']
    wavemap = spec_props['WAVE']
    ext_spec = spec_props['EXT_SPEC']
    tcorr_spec = spec_props['TCORR_SPEC']
    wavemask0 = spec_props['WAVEMASK0']
    wavemask1 = spec_props['WAVEMASK1']
    wavemask2 = spec_props['WAVEMASK2']
    wavemask3 = spec_props['WAVEMASK3']
    max_file = spec_props['MAX_FILE']
    max_snr = spec_props['MAX_SNR']
    wavelim0 = spec_props['WAVELIM0']
    wavelim1 = spec_props['WAVELIM1']
    wavelim2 = spec_props['WAVELIM2']
    wavelim3 = spec_props['WAVELIM3']
    # --------------------------------------------------------------------------
    # setup the figure
    plt.figure(figsize=(12, 12))
    frame0 = plt.subplot2grid((3, 3), (0, 0), colspan=3, rowspan=1)
    frame1 = plt.subplot2grid((3, 3), (1, 0), colspan=3, rowspan=1)
    frame2a = plt.subplot2grid((3, 3), (2, 0), colspan=1, rowspan=1)
    frame2b = plt.subplot2grid((3, 3), (2, 1), colspan=1, rowspan=1)
    frame2c = plt.subplot2grid((3, 3), (2, 2), colspan=1, rowspan=1)

    # set background color
    frame0.set_facecolor(PLOT_BACKGROUND_COLOR)
    frame1.set_facecolor(PLOT_BACKGROUND_COLOR)
    frame2a.set_facecolor(PLOT_BACKGROUND_COLOR)
    frame2b.set_facecolor(PLOT_BACKGROUND_COLOR)
    frame2c.set_facecolor(PLOT_BACKGROUND_COLOR)
    # --------------------------------------------------------------------------
    # Top plot SNR Y
    # --------------------------------------------------------------------------
    # # plot the CCF RV points
    frame0.plot_date(mjd.plot_date, ext_y, fmt='.', alpha=0.5,
                     label=ext_y_label)
    frame0.plot_date(mjd.plot_date, ext_h, fmt='.', alpha=0.5,
                     label=ext_h_label)
    frame0.legend(loc=0, ncol=2)
    frame0.grid(which='both', color='lightgray', ls='--')
    frame0.set(xlabel='Date', ylabel='EXT SNR')

    # --------------------------------------------------------------------------
    # Middle plot - full spectra + tcorr
    # --------------------------------------------------------------------------
    title = (f'Spectrum closest to Median {ext_h_label}'
             f'     SNR:{max_snr}     File: {max_file}')

    frame1.plot(wavemap[wavemask0], ext_spec[wavemask0],
                color='k', label='Extracted Spectrum', lw=0.5)
    if tcorr_spec is not None:
        frame1.plot(wavemap[wavemask0], tcorr_spec[wavemask0],
                    color='r', label='Telluric Corrected', lw=0.5)
        frame1.set_ylim((0, 1.5 * np.nanpercentile(tcorr_spec, 99)))
    frame1.set(xlabel='Wavelength [nm]', ylabel='Flux', xlim=wavelim0)
    frame1.set_title(title, fontsize=10)
    frame1.legend(loc=0, ncol=2)
    frame1.grid(which='both', color='lightgray', ls='--')
    # --------------------------------------------------------------------------
    # Bottom plots - Y, J, H spectra + tcorr
    # --------------------------------------------------------------------------
    masks = [wavemask1, wavemask2, wavemask3]
    frames = [frame2a, frame2b, frame2c]
    limits = [wavelim1, wavelim2, wavelim3]
    # loop around masks and frames and plot the middle plots
    for it in range(len(masks)):
        frame, mask, wavelim = frames[it], masks[it], limits[it]
        frame.plot(wavemap[mask], ext_spec[mask],
                   color='k', label='Extracted Spectrum', lw=0.5)
        if tcorr_spec is not None:
            frame.plot(wavemap[mask], tcorr_spec[mask],
                       color='r', label='Telluric Corrected', lw=0.5)
            ymin_mask = 0.5 * np.nanpercentile(tcorr_spec[mask], 1)
            ymax_mask = 1.5 * np.nanpercentile(tcorr_spec[mask], 99)
            frame.set_ylim((ymin_mask, ymax_mask))
        if it == 0:
            frame.set_ylabel('Flux')
        frame.set(xlabel='Wavelength [nm]', xlim=wavelim)
        frame.set_title(f'Zoom {it + 1}', fontsize=10)
        frame.grid(which='both', color='lightgray', ls='--')

    # --------------------------------------------------------------------------
    # add title
    plt.suptitle(plot_title)
    plt.subplots_adjust(bottom=0.05, left=0.06, right=0.99, hspace=0.3,
                        top=0.95)
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

    num_pp_qc = spec_props['NUM_PP_FILES_FAIL']
    num_ext_qc = spec_props['NUM_EXT_FILES_FAIL']
    num_tcorr_qc = spec_props['NUM_TCORR_FILES_FAIL']

    ext_y = spec_props['EXT_Y']
    ext_h = spec_props['EXT_H']
    first_raw = spec_props['FIRST_RAW']
    last_raw = spec_props['LAST_RAW']
    first_pp = spec_props['FIRST_PP']
    last_pp = spec_props['LAST_PP']
    last_pp_proc = spec_props['LAST_PP_PROC']
    first_ext = spec_props['FIRST_EXT']
    last_ext = spec_props['LAST_EXT']
    last_ext_proc = spec_props['LAST_EXT_PROC']
    first_tcorr = spec_props['FIRST_TCORR']
    last_tcorr = spec_props['LAST_TCORR']
    last_tcorr_proc = spec_props['LAST_TCORR_PROC']
    version_pp = spec_props['PP_VERSION']
    version_ext = spec_props['EXT_VERSION']
    version_tcorr = spec_props['TCORR_VERSION']
    coord_url = spec_props['COORD_URL']
    ra, dec = spec_props['RA'], spec_props['Dec']
    teff, spt = spec_props['Teff'], spec_props['Spectral Type']
    dprtypes = spec_props['DPRTYPES']
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
    # Add RA and Dec
    stat_dict['Description'].append('Coordinates')
    coordstr = _make_url('[FINDER CHART]', coord_url)
    stat_dict['Value'].append(f'({ra},{dec})..... {coordstr}')
    # Add Teff
    stat_dict['Description'].append('Teff')
    stat_dict['Value'].append(teff)
    # Add Spectral type
    stat_dict['Description'].append('Spectral Type')
    stat_dict['Value'].append(spt)
    # Add dprtypes
    stat_dict['Description'].append('DPRTYPES')
    stat_dict['Value'].append(dprtypes)
    # -------------------------------------------------------------------------
    # add number of raw files
    # -------------------------------------------------------------------------
    stat_dict['Description'].append('Total number raw files')
    stat_dict['Value'].append(f'{num_raw}')
    stat_dict['Description'].append('First raw files')
    stat_dict['Value'].append(f'{first_raw}')
    stat_dict['Description'].append('Last raw files')
    stat_dict['Value'].append(f'{last_raw}')

    # -------------------------------------------------------------------------
    # add number of pp files
    # -------------------------------------------------------------------------
    stat_dict['Description'].append('Total number PP files')
    stat_dict['Value'].append(f'{num_pp + num_pp_qc}')
    stat_dict['Description'].append('Number PP files passed QC')
    stat_dict['Value'].append(f'{num_pp}')
    stat_dict['Description'].append('Number PP files failed QC')
    stat_dict['Value'].append(f'{num_pp_qc}')
    stat_dict['Description'].append('First pp file [Mid exposure]')
    stat_dict['Value'].append(f'{first_pp}')
    stat_dict['Description'].append('Last pp file [Mid exposure]')
    stat_dict['Value'].append(f'{last_pp}')
    stat_dict['Description'].append('Last processed [pp]')
    stat_dict['Value'].append(f'{last_pp_proc}')
    stat_dict['Description'].append('Version [pp]')
    stat_dict['Value'].append(f'{version_pp}')
    # -------------------------------------------------------------------------
    # add number of ext files
    # -------------------------------------------------------------------------
    stat_dict['Description'].append('Total number ext files')
    stat_dict['Value'].append(f'{num_ext + num_ext_qc}')
    stat_dict['Description'].append('Number ext files passed QC')
    stat_dict['Value'].append(f'{num_ext}')
    stat_dict['Description'].append('Number ext files failed QC')
    stat_dict['Value'].append(f'{num_ext_qc}')
    stat_dict['Description'].append('First ext file [Mid exposure]')
    stat_dict['Value'].append(f'{first_ext}')
    stat_dict['Description'].append('Last ext file [Mid exposure]')
    stat_dict['Value'].append(f'{last_ext}')
    stat_dict['Description'].append('Last processed [ext]')
    stat_dict['Value'].append(f'{last_ext_proc}')
    stat_dict['Description'].append('Version [ext]')
    stat_dict['Value'].append(f'{version_ext}')
    # -------------------------------------------------------------------------
    # add number of tcorr files
    # -------------------------------------------------------------------------
    stat_dict['Description'].append('Total number tcorr files')
    stat_dict['Value'].append(f'{num_tcorr + num_tcorr_qc}')
    stat_dict['Description'].append('Number tcorr files passed QC')
    stat_dict['Value'].append(f'{num_tcorr}')
    stat_dict['Description'].append('Number tcorr files failed QC')
    stat_dict['Value'].append(f'{num_tcorr_qc}')
    stat_dict['Description'].append('First tcorr file [Mid exposure]')
    stat_dict['Value'].append(f'{first_tcorr}')
    stat_dict['Description'].append('Last tcorr file [Mid exposure]')
    stat_dict['Value'].append(f'{last_tcorr}')
    stat_dict['Description'].append('Last processed [tcorr]')
    stat_dict['Value'].append(f'{last_tcorr_proc}')
    stat_dict['Description'].append('Version [tcorr]')
    stat_dict['Value'].append(f'{version_tcorr}')
    # -------------------------------------------------------------------------
    # add the SNR in Y
    stat_dict['Description'].append('Median SNR Y')
    value = r'{:.2f} :math:`\pm` {:.2f}'.format(med_snr_y, rms_snr_y)
    stat_dict['Value'].append(value)
    # add the SNR in H
    stat_dict['Description'].append('Median SNR H')
    value = r'{:.2f} :math:`\pm` {:.2f}'.format(med_snr_h, rms_snr_h)
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
    reset_mask = lbl_props['RESET_RV']
    # sort data by date
    sort = np.argsort(plot_date)
    plot_date = plot_date[sort]
    vrad = vrad[sort]
    svrad = svrad[sort]
    snr_h = snr_h[sort]
    reset_mask = reset_mask[sort]
    # set background color
    frame[0].set_facecolor(PLOT_BACKGROUND_COLOR)
    frame[1].set_facecolor(PLOT_BACKGROUND_COLOR)
    # --------------------------------------------------------------------------
    # Top plot LBL RV
    # --------------------------------------------------------------------------
    # plot the points
    frame[0].plot_date(plot_date[~reset_mask], vrad[~reset_mask], fmt='.',
                       alpha=0.5, color='green', ls='None')
    frame[0].plot_date(plot_date[reset_mask], vrad[reset_mask], fmt='.',
                       alpha=0.5, color='purple', ls='None')
    # plot the error bars
    frame[0].errorbar(plot_date[~reset_mask], vrad[~reset_mask],
                      yerr=svrad[~reset_mask],
                      marker='o', alpha=0.5, color='green', ls='None',
                      label='Good')
    frame[0].errorbar(plot_date[reset_mask], vrad[reset_mask],
                      yerr=svrad[reset_mask],
                      marker='o', alpha=0.5, color='purple', ls='None',
                      label='Possibly bad (reset rv)')
    # find percentile cuts that will be expanded by 150% for the ylim
    pp = np.nanpercentile(vrad, [10, 90])
    diff = pp[1] - pp[0]
    central_val = np.nanmean(pp)
    # used for plotting but also for the flagging of outliers
    ylim = [central_val - 1.5 * diff, central_val + 1.5 * diff]
    # length of the arrow flagging outliers
    l_arrow = 0.05 * (ylim[1] - ylim[0])
    # store the bad points
    bad_points = []
    # set the arrow properties
    arrowprops = dict(arrowstyle='<-', linewidth=2, color='red')
    arrow = None
    # --------------------------------------------------------------------------
    # flag the low outliers
    low = vrad < ylim[0]
    # get the x and y values of the outliers to be looped over within
    # the arrow plotting
    xpoints = np.array(plot_date[low], dtype=float)
    # x_range = np.nanmax(plot_date) - np.nanmin(plot_date)
    for ix in range(len(xpoints)):
        bad_points.append(ix)
        arrow = frame[0].annotate('',
                                  xy=(xpoints[ix], ylim[0] + l_arrow),
                                  xytext=(xpoints[ix], ylim[0] - l_arrow * 2),
                                  xycoords='data', textcoords='data',
                                  arrowprops=arrowprops)

        # frame[0].arrow(xpoints[ix], ylim[0] + l_arrow * 2, 0, -l_arrow,
        #                color='red', head_width=0.01 * x_range,
        #                head_length=0.25 * l_arrow, alpha=0.5, label='Outliers')
    # same as above for the high outliers
    high = vrad > ylim[1]
    xpoints = np.array(plot_date[high], dtype=float)
    for ix in range(len(xpoints)):
        bad_points.append(ix)

        arrow = frame[0].annotate('',
                                  xy=(xpoints[ix], ylim[1] - l_arrow * 2),
                                  xytext=(xpoints[ix], ylim[1] + l_arrow),
                                  xycoords='data', textcoords='data',
                                  arrowprops=arrowprops)

        # frame[0].arrow(xpoints[ix], ylim[1] - l_arrow * 2, 0, l_arrow,
        #                color='red', head_width=0.01 * x_range,
        #                head_length=0.25 * l_arrow, alpha=0.5, label='Outliers')
    # --------------------------------------------------------------------------
    # setting the plot
    frame[0].set(ylim=ylim)
    frame[0].set(title=plot_title)
    frame[0].grid(which='both', color='lightgray', linestyle='--')
    frame[0].set(ylabel='Velocity [m/s]')
    # only keep one unique labels for legend
    handles, labels = [], []
    raw_handles, raw_labels = frame[0].get_legend_handles_labels()
    for it in range(len(raw_labels)):
        if raw_labels[it] not in labels:
            handles.append(raw_handles[it])
            labels.append(raw_labels[it])
    # --------------------------------------------------------------------------
    # Create a custom legend handle for the arrows
    if arrow is not None:
        arrow_handle = ArrowHandler()
        arrow_handle.arrowprops = arrowprops
        handler_map = {tuple: arrow_handle}
        handles.append((arrow,))
        labels.append('Outliers')
        # add legend
        frame[0].legend(handles, labels, loc=0, handler_map=handler_map)
    else:
        frame[0].legend(handles, labels, loc=0)
    # --------------------------------------------------------------------------
    # Bottom plot SNR
    # --------------------------------------------------------------------------
    # simple plot of the SNR in a sample order. You need to
    # update the relevant ketword for SPIRou
    frame[1].plot_date(plot_date[~reset_mask], snr_h[~reset_mask], fmt='.',
                       alpha=0.5, color='green', ls='None', label='Good')
    frame[1].plot_date(plot_date[reset_mask], snr_h[reset_mask], fmt='.',
                       alpha=0.5, color='purple', ls='None',
                       label='Possibily bad (reset rv)')
    # over plot the bad points from above
    if len(bad_points) > 0:
        bad_points = np.array(bad_points)
        frame[1].plot_date(plot_date[bad_points], snr_h[bad_points], fmt='.',
                           alpha=0.5, color='red', ls='None', label='Outliers')

    # add properties
    frame[1].grid(which='both', color='lightgray', linestyle='--')
    frame[1].set(xlabel='Date')
    frame[1].set(ylabel=snr_h_label)
    # add legend
    frame[1].legend(loc=0)
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
    version_lbl = lbl_props['version']
    num_reset_rv = lbl_props['NUM_RESET_RV']
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
    stat_dict['Description'].append('RV Uncertainty lbl.rdb (25, 50, 75 percentile)')
    stat_dict['Value'].append('{:.2f}, {:.2f}, {:.2f} m/s'.format(*p_sigma))
    # add the absolute deviation
    stat_dict['Description'].append('RV Absolute Deviation lbl.rdb (25, 50, 75 '
                                    'percentile)')
    stat_dict['Value'].append('{:.2f}, {:.2f}, {:.2f} m/s'.format(*v_sigma))
    # add the number of measurements
    stat_dict['Description'].append('Number of lbl.rdb Measurements')
    stat_dict['Value'].append(len(vrad))
    # add the spurious low points
    stat_dict['Description'].append('Number of lbl.rdb Spurious Low Points')
    stat_dict['Value'].append(np.sum(low))
    # add the spurious high points
    stat_dict['Description'].append('Number of lbl.rdb Spurious High Points')
    stat_dict['Value'].append(np.sum(high))
    # add the number of nights
    stat_dict['Description'].append('Number of Nights')
    stat_dict['Value'].append(n_nights)
    # add the number of reset RV points
    stat_dict['Description'].append('Number of Reset RV Points')
    stat_dict['Value'].append(num_reset_rv)
    # add the systemic velocity
    stat_dict['Description'].append('Systemic Velocity')
    stat_dict['Value'].append('{:.2f} m/s'.format(sys_vel))
    # add the Velocity domain considered
    stat_dict['Description'].append('Velocity Domain considered valid')
    stat_dict['Value'].append('{:.2f} to {:.2f} m/s'.format(*vel_domain))
    # add the LBL version
    stat_dict['Description'].append('LBL Version')
    stat_dict['Value'].append(version_lbl)
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

    # sort data by mjd.plot_date
    sort = np.argsort(mjd.plot_date)
    mjd = mjd[sort]
    vrad = vrad[sort]
    svrad = svrad[sort]
    # ylim = ccf_props['ylim']
    # --------------------------------------------------------------------------
    # setup the figure
    fig, frame = plt.subplots(4, 1, figsize=(12, 12))
    # set background color
    frame[0].set_facecolor(PLOT_BACKGROUND_COLOR)
    frame[1].set_facecolor(PLOT_BACKGROUND_COLOR)
    frame[2].set_facecolor(PLOT_BACKGROUND_COLOR)
    frame[3].set_facecolor(PLOT_BACKGROUND_COLOR)
    # --------------------------------------------------------------------------
    # Top plot CCF RV
    # --------------------------------------------------------------------------
    # # plot the CCF RV points
    frame[0].plot_date(mjd.plot_date, vrad, fmt='.', alpha=0.5,
                       color='green', label='Good')
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
    l_arrow = 0.05 * (ylim[1] - ylim[0])
    # set the arrow properties
    arrowprops = dict(arrowstyle='<-', linewidth=2, color='red')
    arrow = None
    # --------------------------------------------------------------------------
    # flag the low outliers
    low = vrad < ylim[0]
    # get the x and y values of the outliers to be looped over within
    # the arrow plotting
    xpoints = np.array(mjd.plot_date[low], dtype=float)
    # x_range = np.nanmax(mjd.plot_date) - np.nanmin(mjd.plot_date)
    for ix in range(len(xpoints)):
        arrow = frame[0].annotate('',
                                  xy=(xpoints[ix], ylim[0] - l_arrow),
                                  xytext=(xpoints[ix], ylim[0] + l_arrow * 2),
                                  xycoords='data', textcoords='data',
                                  arrowprops=arrowprops)

        # frame[0].arrow(xpoints[ix], ylim[0] + l_arrow * 2, 0, -l_arrow,
        #                color='red', head_width=0.01 * x_range,
        #                head_length=0.25 * l_arrow, alpha=0.5, label='Outliers')
    # same as above for the high outliers
    high = vrad > ylim[1]
    xpoints = np.array(mjd.plot_date[high], dtype=float)
    for ix in range(len(xpoints)):
        arrow = frame[0].annotate('',
                                  xy=(xpoints[ix], ylim[1] - l_arrow * 2),
                                  xytext=(xpoints[ix], ylim[1] + l_arrow),
                                  xycoords='data', textcoords='data',
                                  arrowprops=arrowprops)

        # frame[0].arrow(xpoints[ix], ylim[1] - l_arrow * 2, 0, l_arrow,
        #                color='red', head_width=0.01 * x_range,
        #                head_length=0.25 * l_arrow, alpha=0.5, label='Outliers')
    # --------------------------------------------------------------------------
    # setting the plot
    frame[0].set(ylim=ylim)
    frame[0].grid(which='both', color='lightgray', ls='--')
    frame[0].set(xlabel='Date', ylabel='Velocity [m/s]')
    # only keep one unique labels for legend
    handles, labels = [], []
    raw_handles, raw_labels = frame[0].get_legend_handles_labels()
    for it in range(len(raw_labels)):
        if raw_labels[it] not in labels:
            handles.append(raw_handles[it])
            labels.append(raw_labels[it])
    # --------------------------------------------------------------------------
    # Create a custom legend handle for the arrows
    if arrow is not None:
        arrow_handle = ArrowHandler()
        arrow_handle.arrowprops = arrowprops
        handler_map = {tuple: arrow_handle}
        handles.append((arrow,))
        labels.append('Outliers')
        # add legend
        frame[0].legend(handles, labels, loc=0, handler_map=handler_map)
    else:
        frame[0].legend(handles, labels, loc=0)
    # --------------------------------------------------------------------------
    # Middle plot median CCF
    # --------------------------------------------------------------------------
    # mask by xlim
    limmask = (rv_vec > xlim[0]) & (rv_vec < xlim[1])

    frame[1].fill_between(rv_vec[limmask], y1_2sig[limmask], y2_2sig[limmask],
                          color='orange', alpha=0.4)
    frame[1].fill_between(rv_vec[limmask], y1_1sig[limmask], y2_1sig[limmask],
                          color='red', alpha=0.4)
    frame[1].plot(rv_vec[limmask], med_ccf[limmask], alpha=1.0, color='black')
    if has_fit:
        frame[1].plot(rv_vec[limmask], fit[limmask], alpha=0.8,
                      label='Gaussian fit', ls='--')
    frame[1].legend(loc=0)
    frame[1].set(xlabel='RV [km/s]',
                 ylabel='Normalized CCF')
    frame[1].grid(which='both', color='lightgray', ls='--')

    # --------------------------------------------------------------------------
    # Middle plot median CCF residuals
    # --------------------------------------------------------------------------
    if has_fit:
        frame[2].fill_between(rv_vec[limmask], y1_2sig[limmask] - fit[limmask],
                              y2_2sig[limmask] - fit[limmask], color='orange',
                              alpha=0.4, label=r'2-$\sigma$')
        frame[2].fill_between(rv_vec[limmask], y1_1sig[limmask] - fit[limmask],
                              y2_1sig[limmask] - fit[limmask], color='red',
                              alpha=0.4, label=r'1-$\sigma$')
        frame[2].plot(rv_vec[limmask], med_ccf[limmask] - fit[limmask],
                      alpha=0.8, label='Median residual')
        frame[2].legend(loc=0, ncol=3)
        frame[2].set(xlabel='RV [km/s]', ylabel='Residuals [to fit]')
    else:
        frame[2].text(0.5, 0.5, 'No fit to CCF possible',
                      horizontalalignment='center')
        frame[2].legend(loc=0, ncol=3)
        frame[2].set(xlim=[0, 1], ylim=[0, 1], xlabel='RV [km/s]',
                     ylabel='Residuals')
    frame[2].grid(which='both', color='lightgray', ls='--')
    # --------------------------------------------------------------------------
    # Bottom plot median CCF residuals
    # --------------------------------------------------------------------------
    if has_fit:
        frame[3].fill_between(rv_vec[limmask],
                              y1_2sig[limmask] - med_ccf[limmask],
                              y2_2sig[limmask] - med_ccf[limmask], color='orange',
                              alpha=0.4, label=r'2-$\sigma$')
        frame[3].fill_between(rv_vec[limmask],
                              y1_1sig[limmask] - med_ccf[limmask],
                              y2_1sig[limmask] - med_ccf[limmask], color='red',
                              alpha=0.4, label=r'1-$\sigma$')
        frame[3].plot(rv_vec[limmask], med_ccf[limmask] - med_ccf[limmask],
                      alpha=0.8, label='Median residual')
        frame[3].legend(loc=0, ncol=3)
        frame[3].set(xlabel='RV [km/s]', ylabel='Residuals [To Median]')
    else:
        frame[3].text(0.5, 0.5, 'No fit to CCF possible',
                      horizontalalignment='center')
        frame[3].legend(loc=0, ncol=3)
        frame[3].set(xlim=[0, 1], ylim=[0, 1], xlabel='RV [km/s]',
                     ylabel='Residuals [To Median]')
    frame[3].grid(which='both', color='lightgray', ls='--')
    # --------------------------------------------------------------------------
    # add title
    plt.suptitle(plot_title)
    plt.subplots_adjust(hspace=0.2, left=0.1, right=0.99, bottom=0.05,
                        top=0.95)
    # save figure and close the plot
    plt.savefig(plot_path)
    plt.close()


def ccf_stats_table(ccf_props: Dict[str, Any], stat_path: str, title: str):
    from apero.core.math import estimate_sigma
    # get parameters from props
    vrad = ccf_props['dv']
    fwhm = ccf_props['fwhm']
    num_ccf = len(ccf_props['files'])
    num_ccf_qc = len(ccf_props['files_failed'])
    first_ccf = ccf_props['FIRST_CCF']
    last_ccf = ccf_props['LAST_CCF']
    last_ccf_proc = ccf_props['LAST_CCF_PROC']
    version_ccf = ccf_props['CCF_VERSION']
    chosen_mask = ccf_props['chosen_mask']
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
    # add mask used
    stat_dict['Description'].append('Mask used')
    stat_dict['Value'].append(chosen_mask)
    # add systemic velocity
    stat_dict['Description'].append('CCF systemic velocity')
    value = r'{:.2f} :math:`\pm` {:.2f} m/s'.format(sys_vel, err_sys_vel)
    stat_dict['Value'].append(value)
    # add fwhm
    stat_dict['Description'].append('CCF FWHM')
    value = r'{:.2f} :math:`\pm` {:.2f} m/s'.format(ccf_fwhm, err_ccf_fwhm)
    stat_dict['Value'].append(value)
    # add number of files
    stat_dict['Description'].append('Number of CCF files Total')
    stat_dict['Value'].append(num_ccf + num_ccf_qc)
    # add number of files
    stat_dict['Description'].append('Number of CCF passed QC')
    stat_dict['Value'].append(num_ccf)
    # add number of ccf files failed
    stat_dict['Description'].append('Number CCF files failed QC')
    stat_dict['Value'].append(num_ccf_qc)
    # add times
    stat_dict['Description'].append('First ccf file [Mid exposure]')
    stat_dict['Value'].append(f'{first_ccf}')
    stat_dict['Description'].append('Last ccf file [Mid exposure]')
    stat_dict['Value'].append(f'{last_ccf}')
    stat_dict['Description'].append('Last processed [ccf]')
    stat_dict['Value'].append(f'{last_ccf_proc}')
    stat_dict['Description'].append('Version [ccf]')
    stat_dict['Value'].append(f'{version_ccf}')
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
# Object index functions
# =============================================================================
def find_finder_charts(path: str, objname: str) -> Tuple[List[str], List[str]]:
    """
    Find finder charts for this object
    :param path: str, the path to the finder charts directory
    :param objname: str, the object name to locate
    :return:
    """
    # expected directory name
    expected_dir = os.path.join(path, objname)
    # deal with no directory --> no finder files
    if not os.path.exists(expected_dir):
        return [], []
    # storage for list of files
    list_of_files = []
    list_of_descs = []
    # loop around filenames
    for filename in os.listdir(expected_dir):
        # only include APERO finder charts
        if not filename.startswith('APERO_finder_chart_'):
            continue
        # only include pdf files
        if not filename.endswith('.pdf'):
            continue
        # get the finder desc
        description = filename.split(objname)[-1]
        description = description.strip('_')
        description = description.replace('_', '-').replace('.pdf', '')
        # append to list
        list_of_files.append(os.path.join(expected_dir, filename))
        list_of_descs.append(description)
    # return the list of files
    return list_of_files, list_of_descs


def make_finder_download_table(entry, objname, item_save_path, item_rel_path,
                               down_save_path, down_rel_path):
    # get the download base name
    dwn_base_name = f'finder_download_{objname}.txt'
    # get the download table path
    item_path = os.path.join(item_save_path, dwn_base_name)
    # compute the download table
    download_table(entry['find_files'], entry['find_descs'],
                   item_path, down_rel_path,
                   down_save_path, title='Finder charts')
    return item_rel_path + dwn_base_name


# =============================================================================
# Worker functions
# =============================================================================
def debug_mode(debugmode: bool, gsettings: dict):
    # switch off some options in debug mode
    if debugmode:
        gsettings['N_CORES'] = 1
        gsettings['filter objects'] = False
    return gsettings


def get_profiles_from_store(settings: dict, profiles: dict,
                            reprocess: dict) -> Tuple[dict, dict]:
    # construct profiles yaml file name
    path = settings['WORKING']
    profile_yaml_file = os.path.join(path, 'all_profiles.yaml')
    # deal with all profiles existing
    if os.path.exists(profile_yaml_file):
        # read all_profiles.yaml
        all_profiles = read_yaml_file(profile_yaml_file)
        # loop around all profiles and add missing ones to list
        for profile_name in all_profiles:
            if profile_name not in profiles:
                profiles[profile_name] = all_profiles[profile_name]
                reprocess[profile_name] = False
    # re-create the yaml file
    write_yaml(profiles, profile_yaml_file)
    # return updated profiles
    return profiles, reprocess


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


def load_table(filename, **kwargs):
    # attempt to load astropy table
    try:
        return Table.read(filename, **kwargs)
    except Exception as e:
        emsg = 'Cannot load table from: {0}\n\t{1}: {2}'
        eargs = [filename, type(e), str(e)]
        raise ValueError(emsg.format(*eargs))


def save_stats(settings: dict, stats: dict):
    """
    Save the stats for a given profile

    :param settings: dict, the settings for the profile:
    :param stats: dict, the stats to save

    :return:
    """
    # list tables to load
    tables = TABLE_NAMES
    # loop around tables and load them
    for table_name in tables:
        # deal with no table
        if stats['TABLES'][table_name] is None:
            continue
        # otherwise we assume it is an astropy table
        else:
            table = stats['TABLES'][table_name]
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
    profile_stats['TABLES'] = dict()
    profile_stats['OPROPS'] = dict()
    # list tables to load
    tables = ['OBJECT_TABLE', 'RECIPE_TABLE']
    # loop around tables and load them
    for table_name in tables:
        # construct filename for table
        table_filename = f'{table_name}.fits'
        # construct full path
        table_path = os.path.join(settings['DATA'], table_filename)
        # might not be on disk
        if not os.path.exists(table_path):
            profile_stats['TABLES'][table_name] = None
            continue
        # load the table
        profile_stats['TABLES'][table_name] = load_table(table_path)
    # return the profile stats
    return profile_stats


def read_yaml_file(pfilename: str) -> dict:
    """
    Read the yaml file  for the list of profiles using pyyaml

    :param pfilename: str, the filename of the yaml file

    :return:
    """
    # read the yaml file "profiles.yaml" for the list of profiles using pyyaml
    with open(pfilename, 'r') as stream:
        # try to load the yaml file
        try:
            profiles = yaml.safe_load(stream)
        # deal with a yaml error (for now just raise the error)
        except yaml.YAMLError as exc:
            raise exc
    # return the profiles
    return profiles


def write_yaml(dictionary: dict, filename: str):
    """
    Write a dictionary to a yaml file

    :param dictionary:
    :param filename:
    :return:
    """
    # remove file if it exists
    if os.path.exists(filename):
        os.remove(filename)
    # save yaml file
    with open(filename, 'w') as yfile:
        yaml.dump(dictionary, yfile)


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
        if "fallback" in keydict:
            rawvalue = keydict["fallback"]
            if rawvalue is None:
                if dtype == 'float':
                    rawvalue = np.nan
                elif dtype == 'str':
                    rawvalue = "N.A."
                else:
                    # Hard to define an unambiguous value for other types.
                    # Leaving as unsupported until the need arises.
                    raise TypeError('Null (None) fallback value unsupported'
                                    f' for dtype {dtype} (key {key})'
                                    f'\n\tFile: {filename}'
                    )
        else:
            raise ValueError(f'HeaderKey: {key} not found in header'
                             ' and no fallback specified in config'
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


def _match_file(reffile: str, files: List[str]):
    """
    Using a ref file split at the _pp level and try to locate the position
    of this file in the files list

    If ref file does not have _pp we assume it is a raw file and just remove
    the .fits and search based on this

    :param reffile: str,
    :param files:
    :return:
    """
    # force ref to be a basename
    refbasename = os.path.basename(reffile)
    # get the search string
    if '_pp' not in refbasename:
        searchstr = refbasename.split('.fits')[0]
    else:
        searchstr = refbasename.split('_pp')[0]

    # look for each string in list (correct entry should start with it)
    for it, filename in enumerate(files):
        # get basename of filename
        basename = os.path.basename(filename)
        # search for the string, if we find it we can stop here
        if basename.startswith(searchstr):
            return it
    # if we get to here we did not find the string - return None
    return None


def _filter_pids(findex_table: pd.DataFrame, logdbm: Any) -> np.ndarray:
    """
    Filter file index database by pid to find those that passed

    :param findex_table: Table, the file index database
    :param logdbm: the drs_database.LogDatabase instance

    :return: numpy 1D array, a True/False mask the same length as findex_table
    """
    # assume everything failed
    passed = np.zeros(len(findex_table)).astype(bool)
    if len(findex_table) == 0:
        return passed
    # get the pids for these files
    pids = np.array(findex_table['KW_PID'])
    # get the pid
    pid_conds = []
    for pid in pids:
        pid_conds.append(f'PID="{pid}"')
    pid_condition = ' OR '.join(pid_conds)
    # need to crossmatch again recipe log database for QC
    ltable = logdbm.get_entries('PID,PASSED_ALL_QC', condition=pid_condition)
    # get the columsn from the log table
    all_pids = np.array(ltable['PID'])
    # get the passed all qc value
    passed_all_qc = np.array(ltable['PASSED_ALL_QC'])
    # if value is None assume it passed
    null_mask = ~(passed_all_qc == 1)
    null_mask &= ~(passed_all_qc == 0)
    passed_all_qc[null_mask] = 1
    # push into True and False
    all_pass = passed_all_qc.astype(bool)
    # need to loop around all files
    for row in range(len(findex_table)):
        # get the pid for this row
        pid = findex_table['KW_PID'].iloc[row]
        # find all rows that have this pid
        mask = all_pids == pid
        # deal with no entries
        # noinspection PyTypeChecker
        if len(mask) == 0:
            continue
        # if all rows pass qc passed = 1
        if np.sum(all_pass[mask]):
            passed[row] = True
    # return the passed mask
    return passed


def _url_addp(key: str, value: Union[str, List[str], None],
              fmt: str = '&{0}={1}') -> str:
    # deal with value being None
    if value is None:
        return ''
    # deal with a str or a list
    if isinstance(value, str):
        values = [value]
    else:
        values = value
    # storage for output string
    outstr = ''
    # loop around all values
    for value in values:
        if ' ' in value:
            outstr += fmt.format(key, value.replace(' ', '+'))
        else:
            outstr += fmt.format(key, value)
    # return outstr
    return outstr


class ArrowHandler(HandlerPatch):
    def __init__(self):
        # get super call
        super().__init__()
        # set arrow props
        self.arrowprops = None

    def create_artists(self, legend, orig_handle, xdescent, ydescent,
                       width, height, fontsize, trans):
        arrow = FancyArrowPatch((0, 3), (25, 3), **self.arrowprops)
        arrow.set_mutation_scale(10)
        return [arrow]


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # get the args
    _args = get_args()
    profile_filename = _args.profile
    # deal with debug
    if _args.debug:
        debug = True
    else:
        debug = False
    # step 1: read the yaml file
    apero_profiles = read_yaml_file(profile_filename)
    # get global settings
    ari_gsettings = dict(apero_profiles['settings'])
    header_settings = dict(apero_profiles['headers'])
    del apero_profiles['settings']
    del apero_profiles['headers']
    # -------------------------------------------------------------------------
    # deal with reprocessing
    _reprocess = dict()
    # list current profiles
    current_profiles = list(apero_profiles.keys())
    # if filter in this list we remove all other profiles
    if _args.filter not in [None, 'None'] and _args.filter in current_profiles:
        for _profile in current_profiles:
            if _profile == _args.filter:
                _reprocess[_profile] = True
            else:
                _reprocess[_profile] = False
    else:
        for _profile in current_profiles:
            _reprocess[_profile] = True
    # -------------------------------------------------------------------------
    # deal with a reset
    if ari_gsettings['reset']:
        # remove working directory
        shutil.rmtree(ari_gsettings['working directory'], ignore_errors=True)
    # get settings for sync
    ari_settings = get_settings(ari_gsettings)
    # deal with a sync from server
    if ari_gsettings['sync']:
        # step 6: upload to hosting
        print('=' * 50)
        print('Syncing docs...')
        print('=' * 50)
        sync_docs(ari_gsettings, ari_settings)
    # deal with debug mode
    ari_gsettings = debug_mode(debug, ari_gsettings)
    # ----------------------------------------------------------------------
    # look for profiles not in apero_profiles that we have to add
    #  (i.e. they are on the server)
    apero_profiles, _reprocess = get_profiles_from_store(ari_settings,
                                                         apero_profiles,
                                                         _reprocess)
    # ----------------------------------------------------------------------
    # step 2: for each profile compile all stats
    all_apero_stats = dict()
    all_apero_oprops = dict()
    # loop around profiles from yaml file
    for _apero_profname in apero_profiles:
        # sort out settings
        ari_settings = get_settings(ari_gsettings, _apero_profname)
        # add profile name to settings
        apero_profiles[_apero_profname]['profile name'] = _apero_profname
        # we reprocess if the file does not exist or if REPROCESS is True
        if _reprocess[_apero_profname]:
            # print progress
            print('=' * 50)
            print('Compiling stats for profile: {0}'.format(_apero_profname))
            print('=' * 50)
            # get profile
            apero_profile = apero_profiles[_apero_profname]
            # compile stats
            apero_stats = compile_stats(ari_gsettings, ari_settings,
                                        apero_profile, header_settings)
            # add to all_apero_stats
            all_apero_stats[_apero_profname] = apero_stats['TABLES']
            all_apero_oprops[_apero_profname] = apero_stats['OPROPS']
            # -----------------------------------------------------------------
            # Save stats to disk
            save_stats(ari_settings, apero_stats)
        # otherwise we load the stats from disk
        else:
            # print progress
            print('=' * 50)
            print('Loading stats for profile: {0}'.format(_apero_profname))
            print('=' * 50)
            # load stats
            apero_stats = load_stats(ari_settings)
            # add to all_apero_stats
            all_apero_stats[_apero_profname] = apero_stats['TABLES']
            all_apero_oprops[_apero_profname] = apero_stats['OPROPS']
    # ----------------------------------------------------------------------
    # sort out settings
    ari_settings = get_settings(ari_gsettings)
    # ----------------------------------------------------------------------
    # write object index page
    print('=' * 50)
    print('Compling object index page...')
    print('=' * 50)
    compile_obj_index_page(ari_gsettings, ari_settings, all_apero_oprops)
    # ----------------------------------------------------------------------
    # step 4: write markdown files
    print('=' * 50)
    print('Writing markdown files...')
    print('=' * 50)
    write_markdown(ari_gsettings, ari_settings, all_apero_stats)
    # ----------------------------------------------------------------------
    # step 5: compile sphinx files
    print('=' * 50)
    print('Compiling docs...')
    print('=' * 50)
    compile_docs(ari_gsettings, ari_settings)
    # ----------------------------------------------------------------------
    # step 6: upload to hosting
    print('=' * 50)
    print('Uploading docs...')
    print('=' * 50)
    upload_docs(ari_gsettings, ari_settings)

# =============================================================================
