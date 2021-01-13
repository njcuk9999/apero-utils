"""
Test abstract class and factory for APERO.

@author: vandal
"""
import os
import glob
from datetime import datetime
from abc import ABC, abstractmethod
from typing import Optional, Union, List, Tuple

import numpy as np
import pandas as pd
from astropy.table import Table
from jinja2 import Environment, FileSystemLoader, select_autoescape

import apero_tests_func as atf
from apero.core import constants


class Test(ABC):
    """Test."""

    def __init__(self, inst: str = 'SPIROU', setup: Optional[str] = None):
        """__init__.

        :param inst: Instrument name used in APERO (Default is SPIROU).
        :type inst: str
        :param setup: Path to apero setup. Using DRS_UCONFIG if None.
        :type setup: Optional[str]
        """
        self._params = constants.load(inst)

        # setup path
        if setup is None:
            self._setup = os.environ['DRS_UCONFIG']
        else:
            self._setup = setup

        self._instrument = self.params['INSTRUMENT']  # instrument

        self._date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # test date

    # =========================================================================
    # Properties
    # =========================================================================
    @property
    def setup(self) -> str:
        """setup.

        :rtype: str
        """
        return self._setup

    @property
    def instrument(self) -> str:
        """instrument.

        :rtype: str
        """
        return self._instrument

    @property
    def params(self) -> constants.param_functions.ParamDict:
        """params.

        :rtype: constants.param_functions.ParamDict
        """
        return self._params

    @property
    def date(self) -> str:
        """date.

        :rtype: str
        """
        return self._date

    # =========================================================================
    # Public methods
    # =========================================================================
    def gen_html(self, html_dict: dict):
        """Generate HTML summary from jinja2 template.

        :param html_dict: dict with all values used in html template
        :type html_dict: dict
        """

        # Jinja2 env
        # TODO: Use package or define path in a main/init file?
        env = Environment(
            loader=FileSystemLoader('../templates'),
            autoescape=select_autoescape(['html', 'xml'])
        )

        # Create template for test
        template = env.get_template('.'.join([self.test_id, 'html']))

        html_text = template.render(html_dict)

        output_path = os.path.join('..',
                                   'out',
                                   self.test_id,
                                   '.'.join([self.test_id, 'html']),
                                   )
        with open(output_path, 'w') as f:
            f.write(html_text)

    # =========================================================================
    # Abstract methods common to all tests
    # =========================================================================
    @property
    @abstractmethod
    def name(self) -> str:
        """Test full unique name."""

    @property
    @abstractmethod
    def test_id(self) -> str:
        """Test short name (ID)."""

    @property
    @abstractmethod
    def recipe(self) -> str:
        """Recipe checked by the tests"""

    @property
    @abstractmethod
    def output_list(self) -> List[str]:
        """List of output string patterns

        :return: output_list
        :rtype: List[str]
        """

    @property
    @abstractmethod
    def fibers(self):
        """Fibers to search in outputs"""

    @abstractmethod
    def runtest(self):
        """Method called to run the tests on APERO products"""


class CalibTest(Test):
    """CalibTest."""

    def __init__(self,
                 inst: str = 'SPIROU',
                 setup: Optional[str] = None,
                 logdir: Union[str, list] = 'night'):
        """__init__.

        :param inst:
        :type inst: str
        :param setup:
        :type setup: Optional[str]
        :param logdir:
        :type logdir: str
        """

        super().__init__(inst=inst, setup=setup)

        # List of all reduced nights

        # Where logs can be found
        if logdir == 'night':
            logdirs = self.reduced_nights.copy()
        elif isinstance(logdir, str):
            logdirs = [logdir]
        elif isinstance(logdir, list):
            logdirs = logdir
        else:
            raise TypeError('logdir must be a list or a string')

        # Series of output files
        self._output_files = self._load_output()

        # Handle logs (this creates log_df and missing_logfits properties)
        self._log_df, self._missing_logs = self._load_log(logdirs)


        # Get aligned counts for logs
        self._output_num_align, self._log_num_align = self._align_counts()

        # Load calib db
        self._master_calib_df = self._load_calibdb_entries()
        self._output_calibdb = self._load_calibdb_output()

    # =========================================================================
    # General properties
    # =========================================================================
    @property
    def reduced_path(self) -> str:
        """Path to reduced directory

        :return: reduced_path
        :rtype: str
        """
        return self.params['DRS_DATA_REDUC'].rstrip(os.path.sep)

    @property
    def calibdb_path(self) -> str:
        """Path do calibDB directory.

        :return: calibdb_path
        :rtype: str
        """
        return self.params['DRS_CALIB_DB'].rstrip(os.path.sep)

    @property
    def reduced_nights(self) -> List[str]:
        """List of reduced night directories

        :return: reduced_nights
        :rtype: List[str]
        """
        return atf.list_nights(self.reduced_path)

    @property
    def output_files(self) -> pd.Series:
        """output_files.

        :return: output_files
        :rtype: pd.Series
        """
        return self._output_files

    @property
    def log_df(self) -> pd.DataFrame:
        """Dataframe with log information

        :return: log_df
        :rtype: pd.DataFrame
        """
        return self._log_df

    @property
    def master_calib_df(self) -> pd.DataFrame:
        """master_calib_df.

        :rtype: pd.DataFrame
        """
        return self._master_calib_df

    @property
    def output_calibdb(self) -> pd.Series:
        """output_calibdb.

        :rtype: pd.Series
        """
        return self._output_calibdb

    @property
    def calib_output_dict(self) -> dict:
        return dict(zip(self.calibdb_list, self.output_list))

    # =========================================================================
    # Properties derived from outputs
    # =========================================================================
    @property
    def output_night_num(self) -> pd.Series:
        """output_night_num.

        :returns: multi-index series with count of each output per night
        :rtype: pd.Series
        """
        ind_name = self.output_files.index.names
        return self.output_files.groupby(ind_name).size()

    @property
    def output_num_total(self) -> pd.Series:
        """Total umber of outputs for each pattern.

        :rtype: pd.Series
        """
        return self.output_files.groupby('PATTERN').size()

    @property
    def output_num_unique(self) -> pd.Series:
        """output_num_unique.

        :returns: series with total unique count of each output
        :rtype: pd.Series
        """
        return self.output_files.groupby('PATTERN').nunique()

    @property
    def log_num_night(self) -> pd.Series:
        """Number of log entries per night.

        :rtype: pd.Series
        """
        return self.log_df.groupby('DIRECTORY').size()

    @property
    def log_tot_num(self) -> int:
        """Total number of log entries.

        :rtype: int
        """
        return self.log_df.shape[0]

    @property
    def log_num_align(self) -> pd.Series:
        """log_num_align.

        emtpy nights fille with zeros

        :returns: series of log count with indexes aligned to MultiIndex output
        :rtype: pd.Series
        """
        return self._log_num_align

    @property
    def output_num_align(self) -> pd.Series:
        """log_num_align.

         series of log count with indexes aligned with logs
        (empty filled with 0)

        :rtype: pd.Series
        """
        return self._output_num_align

    @property
    def master_mask(self) -> pd.Series:
        return self.log_df['RUNSTRING'].str.contains('--master')

    @property
    def master_nights(self) -> np.ndarray:
        return self.master_mask[self.master_mask].index.values

    @property
    def master_recipe_num_logfits(self) -> int:
        return  self.master_nights.size

    @property
    def recipe_num_logfits(self) -> int:
        return self.log_tot_num - self.master_recipe_num_logfits

    @property
    def output_missing(self) -> pd.DataFrame:
        """get_missing_output."""

        # Get missing outputs that are in log
        log_mask = self.log_num_align > 0
        output_mask = self.output_num_align == 0
        missing_mask = log_mask & output_mask
        output_missing = self.output_num_align.index.to_frame()[missing_mask]
        output_missing = output_missing.reset_index(drop=True)

        return output_missing

    @property
    def tot_num_entry(self) -> pd.Series:
        """Total number of entries in calibdb.

        :rtype: pd.Series
        """
        return self.master_calib_df.groupby('KEY').size()

    @property
    def master_num_entry(self) -> pd.Series:
        """Number of 'master' entries in calibdb.

        :rtype: pd.Series
        """
        master_mask = self.master_calib_df.MASTER
        master_calib_group = self.master_calib_df[master_mask].groupby('KEY')

        return master_calib_group.size()

    @property
    def output_num_entry(self) -> pd.Series:
        """Total number of entries per output in calibdb (minus 'master' ones).

        :rtype: pd.Series
        """
        return self.tot_num_entry - self.master_num_entry

    @property
    def output_num_calibdb(self) -> pd.Series:
        """Number of outputs in calibdb.

        :rtype: pd.Series
        """
        return self.output_calibdb.groupby('KEY').size()

    @property
    def calib_missing_mask(self) -> pd.Series:
        """calib_missing_mask.

        :rtype: pd.Series
        """
        return ~self.master_calib_df.FILE.isin(self.output_calibdb)

    @property
    def log_qc_failed(self) -> pd.DataFrame:
        """log_qc_failed.

        :rtype: pd.DataFrame
        """
        return self.log_df[~self.log_df.PASSED_ALL_QC]

    @property
    def log_ended_false(self) -> pd.DataFrame:
        """log_ended_false.

        :rtype: pd.DataFrame
        """
        return self.log_df[~self.log_df.ENDED]

    # =========================================================================
    # Generators methods to derive things that are too heavy to calculate time
    # =========================================================================
    def _load_output(self) -> pd.Series:
        """Get output files
        Organize output files by output type and night using two indexing
        levels in pandas

        :param fibers: Fiber name(s) (str or list) to iterate files over
        :type fibers: Optional[Union[list,str]]

        :return: Output series of all output files. Each file has two indices
                 (output pattern, night date)
        :rtype: pd.Series
        """

        inds = pd.MultiIndex.from_product(  # multi-index output and night
                [self.output_list, self.reduced_nights],
                names=['PATTERN', 'DIRECTORY'],
                )
        sep = os.path.sep
        paths = pd.Series(self.reduced_path  # Series of path per night+output
                          + sep + inds.get_level_values('DIRECTORY')
                          + sep + inds.get_level_values('PATTERN'),
                          index=inds
                          )

        # Expand patterns
        if self.fibers is None:
            files = paths.apply(glob.glob)  # Get file list for each pattern
            files = files.explode()         # Pop lists in individual cells
        else:
            files = paths.apply(CalibTest._fiberglob, fibers=self.fibers)
            for _ in self.fibers:
                files = files.explode()  # Nested lists so explode Nfiber times

        # Keep filename without path
        files = files.dropna()
        files = files.apply(os.path.basename)

        return files

    @staticmethod
    def _fiberglob(pattern: str,
                   fibers: Union[List[str], Tuple[str]] = ('AB', 'C'),
                   ) -> list:
        """Equivalent of glob but iterating over multiple fibers:

        :param pattern:
        :type pattern: str
        :param fibers:
        :type fibers: Union[List[str], Tuple[str]]

        :returns: Nested list of glob outputs (one per fiber).
        :rtype: list
        """
        return [glob.glob(pattern.format(FIBER=f)) for f in fibers]

    def _load_log(self, logdirs: List[str]) -> Tuple[pd.DataFrame, list]:
        """_gen_log_df.

        Generate a dataframe with log for each directory, with DIRECTORY as
        index. Only logs of the recipe being tested are kept.
        """
        # Get all logs in a dataframe
        allpaths = [os.path.join(self.reduced_path, ld, 'log.fits')
                    for ld in logdirs]
        paths = [p for p in allpaths if os.path.isfile(p)]
        dfs = [Table.read(f).to_pandas() for f in paths]
        log_df = pd.concat(dfs)

        # Decode strings (fits from astropy are encoded)
        for col in log_df.columns:
            # try/except for non-string columns
            try:
                log_df[col] = log_df[col].str.decode('utf-8')
            except AttributeError:
                pass

        # Use DIRECTORY as index and keep onl relevant entries
        log_df = log_df.set_index(['DIRECTORY'])
        log_df = log_df.loc[log_df.RECIPE == self.recipe]

        # Store missing logs in another list
        missing_logs = [p for p in allpaths if not os.path.isfile(p)]

        return log_df, missing_logs

    def _load_calibdb_entries(self) -> pd.DataFrame:
        """_load_calibdb.

        :rtype: pd.DataFrame
        """
        # Read master_calib_{instument}.txt
        master_calib_path = os.path.join(  # Full path: dir and file
                                self.calibdb_path,
                                'master_calib_{}.txt'.format(self.instrument)
                                )
        f = open(master_calib_path, "r")
        master_calib_txt = f.read()

        # Find where DRS processed starts (keep first entry with -1)
        nprocessed = master_calib_txt.index('# DRS processed')
        index_start = master_calib_txt[:nprocessed].count('\n') - 1

        # Load dataframe and keep only badpix entries
        master_calib_df = pd.read_csv(
                master_calib_path,
                header=index_start,
                delimiter=' ',
                usecols=[0, 1, 2, 3],
                names=['KEY', 'MASTER', 'NIGHT', 'FILE'],
                index_col=0,  # Use KEY as index
                )
        if self.fibers is not None:
            calib_inds = [k.format(FIBER=f)
                          for k in self.calibdb_list
                          for f in self.fibers
                          ]
        else:
            calib_inds = self.calibdb_list

        master_calib_df.MASTER = master_calib_df.MASTER.astype(bool)
        master_calib_df = master_calib_df.loc[calib_inds]

        # Replace fibers by {FIBER} to match other tests
        # NOTE: Assumes fibers are at the end with \\b
        if self.fibers is not None:
            master_calib_df.index = master_calib_df.index.str.replace(
                    '_('+'|'.join(self.fibers)+')\\b',
                    '_{FIBER}'
                    )

        return master_calib_df

    def _load_calibdb_output(self) -> pd.Series:
        """_load_calibdb_output.

        :rtype: pd.Series
        """
        # Load output patterns that are in calibdb
        calib_ind = pd.Index(self.calibdb_list, name='KEY')
        out_patterns = pd.Series(self.output_list[:len(self.calibdb_list)],
                                 index=calib_ind)
        sep = os.path.sep
        output_calib_paths = pd.Series(self.calibdb_path + sep + out_patterns,
                                       index=calib_ind
                                       )

        # Expand patterns
        if self.fibers is None:
            files = output_calib_paths.apply(glob.glob)  # Get file list for each pattern
            files = files.explode()         # Pop lists in individual cells
        else:
            files = output_calib_paths.apply(CalibTest._fiberglob,
                                             fibers=self.fibers)
            for _ in self.fibers:
                files = files.explode()  # Nested lists so explode Nfiber times

        output_calibdb = files.apply(os.path.basename)

        return output_calibdb


    def _align_counts(self) -> Tuple[pd.DataFrame]:
        """Align output and log, fill with zeros where not intersecting.

        :param output_num: output counts per night and output
        :param log_num: log counts per night

        :returns: aligned output and log counts, in this order
        :rtype: Tuple[pd.DataFrame]
        """
        log_num_align, output_num_align = self.log_num_night.align(
                                                        self.output_night_num,
                                                        fill_value=0,
                                                        )
        # Filling might switch to float
        log_num_align = log_num_align.astype(int)
        output_num_align = output_num_align.astype(int)

        return (output_num_align, log_num_align)

    # =========================================================================
    # Abstract methods common to all calib tests
    # =========================================================================
    @property
    @abstractmethod
    def calibdb_list(self) -> List[str]:
        """List of calibDB entries

        :return: calibdb_list
        :rtype: List[str]
        """
