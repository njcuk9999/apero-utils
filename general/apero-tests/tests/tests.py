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
from astropy.io import fits
from jinja2 import Environment, FileSystemLoader, select_autoescape
from apero.core import constants

from . import utils as ut
from . import TEMPLATEDIR, OUTDIR


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
        env = Environment(
            loader=FileSystemLoader(TEMPLATEDIR),
            autoescape=select_autoescape(['html', 'xml'])
        )

        # Create template for test
        template = env.get_template('.'.join([self.test_id, 'html']))

        html_text = template.render(html_dict)

        output_path = os.path.join(OUTDIR,
                                   self.test_id,
                                   '.'.join([self.test_id, 'html']))
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

    # pylint: disable=too-many-public-methods
    # Can remove a few very simple properties with v0.7 after making sure
    # there are no major changes with the SQL databases vs fits file
    # (there should not be as we will probably still use a pd dataframe)

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
        # Missing logs to _ because not used
        self._log_df, _ = self._load_log(logdirs)

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
        return ut.list_nights(self.reduced_path)

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
        """Mapping from calibdb entries to output patterns.

        If there are less entries or outputs, only the first entries are
        used.

        :rtype: dict
        """
        return dict(zip(self.calibdb_list, self.output_list))

    # =========================================================================
    # Properties derived from outputs
    # =========================================================================
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
        """master_mask.

        :rtype: pd.Series
        """
        return self.log_df['RUNSTRING'].str.contains('--master')

    @property
    def master_nights(self) -> np.ndarray:
        """master_nights.

        :rtype: np.ndarray
        """
        return self.master_mask[self.master_mask].index.values

    @property
    def master_recipe_num_logfits(self) -> int:
        """master_recipe_num_logfits.

        :rtype: int
        """
        return  self.master_nights.size

    @property
    def recipe_num_logfits(self) -> int:
        """recipe_num_logfits.

        :rtype: int
        """
        return self.log_tot_num - self.master_recipe_num_logfits

    @property
    def output_missing(self) -> pd.DataFrame:
        """output_missing.

        :rtype: pd.DataFrame
        """

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
        f = open(master_calib_path, 'r')
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
        log_num_night = self.log_df.groupby('DIRECTORY').size()
        ind_name = self.output_files.index.names
        output_num_night = self.output_files.groupby(ind_name).size()
        log_num_align, output_num_align = log_num_night.align(
                                                        output_num_night,
                                                        fill_value=0,
                                                        )
        # Filling might switch to float
        log_num_align = log_num_align.astype(int)
        output_num_align = output_num_align.astype(int)

        return (output_num_align, log_num_align)

    # =========================================================================
    # Checks and stop common to all calibs
    # =========================================================================
    def check_qc(self, ncheck: int = 0) -> dict:
        """check_qc."""

        # Check passed QC and ENDED, can then access easily later
        num_logfits_qc_failed = self.log_qc_failed.shape[0]

        # check: QC
        if num_logfits_qc_failed == 0:
            comments_check_qc = ''
            inspect_check_qc = ''
        else:
            comments_check_qc = 'One or more recipe have failed QC.'
            log_reset = self.log_qc_failed.reset_index()
            data_dict_check_qc = {
                    'Night': log_reset.DIRECTORY.values,
                    'QC_STRING': self.log_qc_failed.QC_STRING.values,
                    }
            inspect_check_qc = ut.inspect_table(
                    self.test_id,
                    f'check{ncheck}',
                    data_dict_check_qc,
                    'Nights that Failed Quality Control'
                    )

        return comments_check_qc, inspect_check_qc

    def check_ended(self, ncheck: int = 0) -> dict:
        """check_ended."""
        num_logfits_ended_false = self.log_ended_false.shape[0]

        # check: ENDED
        if num_logfits_ended_false == 0:
            comments_check_ended = ''
            inspect_check_ended = ''
        else:
            comments_check_ended = 'One or more recipe have failed to finish.'
            log_reset = self.log_ended_false.reset_index()
            data_dict_check_ended = {
                    'Night': log_reset.DIRECTORY.values,
                    'ERRORS': self.log_ended_false.ERRORS.values,
                    'LOGFILE': self.log_ended_false.LOGFILE.values,
                    }
            inspect_check_ended = ut.inspect_table(
                    self.test_id,
                    f'check{ncheck}',
                    data_dict_check_ended,
                    'Nights that Failed to Finish'
                    )

        return comments_check_ended, inspect_check_ended

    def check_qc_plot(self, ncheck: int = 0) -> dict:
        """check_qc_plot.

        :rtype: dict
        """
        qc_names = self.log_df.QC_NAMES.str.split(r'\|\|', expand=True).iloc[0]
        qc_values = self.log_df.QC_VALUES.str.split(r'\|\|', expand=True)
        qc_values.columns = qc_names
        # NOTE: .convert_dtypes will do in pd versions >= 1.0.0
        float_mask = ~qc_values.isin(['True', 'False']).any()
        qc_values = qc_values.loc[:, float_mask].astype(float)

        data_dict_check_qc_plot = {'Night': qc_values.index.tolist()}
        for key, series in qc_values.iteritems():
            data_dict_check_qc_plot[key] = series.tolist()

        inspect_check_qc_plot = ut.inspect_plot(
                    self.test_id,
                    f'check{ncheck}',
                    data_dict_check_qc_plot,
                    f'{self.recipe}.py Quality Control'
                    )

        return inspect_check_qc_plot

    def check_duplicates(self) -> Tuple[str, pd.DataFrame]:
        """check_duplicates."""

        output_dup_mask = self.log_num_align > self.output_num_align
        if output_dup_mask.any():
            output_dup_ind = output_dup_mask[output_dup_mask]
            output_dup_ind = output_dup_ind.index.get_level_values(
                                                                'DIRECTORY'
                                                                ).unique()

            # Check if '--master' in runstring
            # Get subset of master in output_dup_ind
            if self.master_recipe_num_logfits > 0:
                comments_check_dup = ('An additional {0} recipe with '
                        '--master=True was called in the master '
                        'directory {1}.'.format(self.master_recipe_num_logfits,
                                                self.master_nights)
                        )

            # Check for non-master duplicates for each pattern
            dir_kwd = 'DIRECTORY'
            master_mask_group = self.master_mask.astype(int).groupby(dir_kwd)
            master_num_night = master_mask_group.sum()
            dup_num_night = master_mask_group.size()
            dup_not_master = dup_num_night - master_num_night
            log_dup_align, output_dup_align = dup_not_master.align(
                                        self.output_num_align[output_dup_mask]
                                            )
            true_dup = log_dup_align - output_dup_align  # number of dups
            true_dup = true_dup[true_dup > 0]  # Keep only non-zero
            true_dup.name = 'COUNT'
            true_dup = true_dup.reset_index()

        else:
            comments_check_dup = ''

        return comments_check_dup, true_dup

    def stop_output_log(self, dup: pd.DataFrame, nstop: int = 0) -> dict:
        """stop_output_log.

        Unique output == log?

        :param dup: true duplicate outputs (i.e. not master)
        :type dup: pd.DataFrame
        :param nstop: stop id number
        :type nstop: int
        :rtype: dict
        """
        if (self.output_num_unique == self.recipe_num_logfits).all():
            color = 'Lime'
            result = 'Yes'
            comment = ''
            inspect = ''
            data_dict = {}

        elif (self.output_num_unique < self.recipe_num_logfits).any():

            color = 'Yellow'
            result = 'No'
            # if missing output
            # NOTE: could there be both missing output and duplicate?
            if self.output_missing.shape[0] > 0:

                comment = 'One or more output were not produced.'
                # NOTE: might need arrays instead of list
                data_dict = {
                        'Night': self.output_missing.DIRECTORY.values,
                        'File name': self.output_missing.PATTERN.values,
                         }
                inspect = ut.inspect_table(
                        self.test_id,
                        f'stop{nstop}',
                        data_dict,
                        f'Missing Outputs in {self.reduced_path}'
                        )
            # if duplicates
            else:
                comment = ('Recipe called multiple times producing the '
                           'same outputs.')
                data_dict = {
                        'Night': self.dup.DIRECTORY.values,
                        'File name': dup.PATTERN.values,
                        'Occurrence': dup.COUNT.values,
                        }

                inspect_msg = ('Same Recipe Called Twice '
                               'or More Producing the Same Outputs in '
                               f'{self.reduced_path}'
                               )
                inspect = ut.inspect_table(
                        self.test_id,
                        f'stop{nstop}',
                        data_dict,
                        inspect_msg
                        )

        else:
            color = 'Red'
            result = 'No'
            comment = ('The number of unique output files should always be '
                       'smaller than or equal to the number of recipes '
                       'called.')
            inspect = ''
            data_dict = {}

        stop_dict = {
                'data': data_dict,
                'comment': comment,
                'inspect': inspect,
                'color': color,
                'result': result,
                }

        return stop_dict

    def stop_calibdb(self, calib_dup: pd.DataFrame, nstop: int = 0) -> dict:
        """stop_calibdb.

        Calibdb output == calibdb entries?

        :param calib_dup: duplicates in calibdb
        :type calib_dup: pd.DataFrame
        :param nstop: stop id number
        :type nstop: int
        :rtype: dict
        """
        if (self.output_num_calibdb == self.output_num_entry).all():
            color = 'Lime'
            result = 'Yes'
            comment = ''
            inspect = ''
            data_dict = {}

        elif (self.output_num_calibdb < self.output_num_entry).any():

            color = 'Yellow'
            result = 'No'

            # NOTE: Could we have both missing and dup?
            if self.calib_missing_mask.any():

                # missing outputs
                comment = 'One or more outputs are not in the calibDB.'
                missing_mask = self.calib_missing_mask
                missing_calibdb_output = self.master_calib_df[missing_mask]
                data_dict = {
                        'Night': missing_calibdb_output.NIGHT.values,
                        'File name': missing_calibdb_output.FILE.values,
                        }
                inspect = ut.inspect_table(
                                self.test_id,
                                f'stop{nstop}',
                                data_dict,
                                f'Missing Output in {self.calibdb_path}'
                                )
            else:

                # duplicates
                comment = ('Some entries in '
                           f'master_calib_{self.instrument}.txt '
                           'are identical.')
                data_dict = {
                        'Night': calib_dup.NIGHT.values,
                        'File name': calib_dup.FILE.values,
                        'Occurrence': calib_dup.COUNT.values,
                         }
                inspect = ut.inspect_table(
                        self.test_id,
                        f'stop{nstop}',
                        data_dict,
                        ('Duplicate Entries in '
                         f'master_calib_{self.instrument}.txt'
                         )
                        )

        else:
            color = 'Red'
            result = 'No'
            comment = ('The calibDB should not have more output files '
                       'than what was produced.')
            inspect = ''
            data_dict = {}

        stop_dict = {
                'data': data_dict,
                'comment': comment,
                'inspect': inspect,
                'color': color,
                'result': result,
                }

        return stop_dict

    def get_missing_previous_calib(self) -> pd.DataFrame:
        """get_missing_calib.

        :rtype: pd.DataFrame
        """

        # Load header of all output files (one output pattern only)
        full_paths = (self.reduced_path
                      + os.path.sep
                      + self.output_files.index.get_level_values('DIRECTORY')
                      + os.path.sep
                      + self.output_files)
        headers = full_paths.loc[self.output_list[0]].apply(fits.getheader)
        header_df = pd.DataFrame(headers.tolist(), index=headers.index)

        # Keep only matching PIDs
        # Not good: header_df = header_df[header_df.DRSPID.isin(self.log_df.PID)]
        log_pid_dir = self.log_df.reset_index().set_index('PID').DIRECTORY
        log_nights = log_pid_dir.loc[header_df.DRSPID]  # log nights for PIDs
        header_df = header_df[(log_nights == header_df.index).values]

        # Keep only calib columns
        used_calibs = [p
                       for p in self.previous_calibs
                       if p in header_df.columns]
        header_df = header_df[used_calibs]  # Keep calibs

        # Get masks (used and exists) and project condition on nights (axis=1)
        none_mask = (header_df == 'None')  # calib not used
        prefix = (self.reduced_path + os.path.sep
                  + header_df.index + os.path.sep)
        isfile_mask = header_df.radd(prefix, axis=0).applymap(os.path.isfile)
        missing_mask = ~(isfile_mask | none_mask)

        # Check nights where 1) used and 2) does not exists for each output
        missing_calib_all = header_df[missing_mask]
        missing_calib_list = [missing_calib_all[k]
                              for k in missing_calib_all.columns]
        missing_calib = pd.concat(missing_calib_list).sort_index()
        missing_calib = missing_calib.dropna()  # 2D mask yields NaNs if false
        missing_calib.name = 'FILE'
        missing_calib = missing_calib.reset_index()
        missing_calib = missing_calib.rename(columns={'DIRECTORY': 'LOC_DIR'})

        # Find calibration nights used
        pattern = (self.reduced_path
                   + os.path.sep + '*'
                   + os.path.sep + missing_calib.FILE)
        calib_path = pattern.apply(glob.glob).str[0]  # First glob for each
        calib_dir = calib_path.str.split(os.path.sep).str[-2]
        missing_calib['CALIB_DIR'] = calib_dir

        return missing_calib

    @staticmethod
    def check_previous_calib(missing_calib, ncheck: int = 0) -> dict:
        """check_previous_calib.

        :param missing_calib:
        :param ncheck:
        :type ncheck: int
        :rtype: dict
        """

        if missing_calib.shape[0] == 0:
            comments_missing_calib = ''
            inspect_missing_calib = ''
        else:
            comments_missing_calib = ('One or more recipe outputs'
                                      ' used the bad pixel calibrations from'
                                      ' another night')
            data_dict_missing_calib = {
                    'Recipe Night': missing_calib.LOC_DIR.values,
                    'Calibration file name': missing_calib.FILE.values,
                    'Calibration Night': missing_calib.CALIB_DIR.values,
                    }
            inspect_missing_calib = ut.inspect_table(
                        self.test_id,
                        f'check{ncheck}',
                        data_dict_missing_calib,
                        ('Night where the Recipe Output '
                         'used Calibrations from Another Night')
                        )

        return comments_missing_calib, inspect_missing_calib

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

    @property
    @abstractmethod
    def previous_calibs(self):
        """List of previous calib entries.

        :rtype: List[str]
        """
