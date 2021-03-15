"""
Test abstract class and factory for APERO.

@author: vandal
"""
import glob
import os
from abc import ABC, abstractmethod
from datetime import datetime
from typing import Optional, Union, List, Tuple

import pandas as pd
from apero.core import constants
from apero.core.instruments.spirou.recipe_definitions import recipes
from astropy.io import fits
from astropy.table import Table
from jinja2 import Environment, FileSystemLoader, select_autoescape

from .log import Log
from . import TEMPLATEDIR, OUTDIR
from . import utils as ut

RECIPE_DICT = dict(zip(list(map(lambda x: x.name, recipes)), recipes))


class Test(ABC):
    """Test."""

    @property
    @abstractmethod
    def output_list(self) -> List[str]:
        """List of output string patterns

        :return: output_list
        :rtype: List[str]
        """

    @property
    @abstractmethod
    def fibers(self) -> Optional[List[str]]:
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
                 logdir: Union[str, list] = 'night',
                 ):
        """__init__.

        :param inst:
        :type inst: str
        :param setup:
        :type setup: Optional[str]
        :param logdir:
        :type logdir: str
        """

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
        log_df, _ = self._load_log(logdirs)
        recipes = [self.recipe, f'cal_extract_{self.instrument.lower()}']
        self._log_all = Log(log_df, recipe=recipes)             # All recipes
        self._log_recipe = Log(log_df, recipe=recipes[0])   # recipe in test
        self._log_extract = Log(log_df, recipe=recipes[1])  # extracts called

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
    def log_all(self) -> Log:
        """log with all entries for test.

        :rtype: Log
        """
        return self._log_all

    @property
    def log_extract(self) -> Log:
        """log with entries of the called extract recipes.

        :rtype: Log
        """
        return self._log_extract

    @property
    def log_recipe(self) -> Log:
        """log with entries of the test recipe (self.recipe).

        :rtype: Log
        """
        return self._log_recipe

    @property
    def log(self) -> Log:
        """log entries used to compare with outputs.

        This is log_recipe if extract calls are used and log_recipe otherwise.
        Useful for checks that require comparison with output

        :rtype: Log
        """
        return self._log_extract if self.calls_extract else self._log_recipe

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
        """Total number of outputs for each pattern.
        :rtype: pd.Series
        """
        if not self.output_files.empty:
            return self.output_files.groupby('PATTERN').size()
        else:
            return pd.Series(0, index=self.output_list)

    @property
    def output_num_unique(self) -> pd.Series:
        """output_num_unique.

        :returns: series with total unique count of each output
        :rtype: pd.Series
        """
        if not self.output_files.empty:
            return self.output_files.groupby('PATTERN').nunique()
        else:
            return pd.Series(0, index=self.output_list)

    @property
    def output_missing(self) -> pd.DataFrame:
        """output_missing.

        :rtype: pd.DataFrame
        """

        # Align log and output
        log_num_align, output_num_align = self.log.get_align_count(
                                                        self.output_num_night
                                                        )

        # Get missing outputs that are in log
        log_mask = log_num_align > 0
        output_mask = output_num_align == 0
        missing_mask = log_mask & output_mask
        output_missing = output_num_align.index.to_frame()[missing_mask]
        output_missing = output_missing.reset_index(drop=True)

        return output_missing

    @property
    def tot_num_entry(self) -> pd.Series:
        """Total number of entries in calibdb.

        :rtype: pd.Series
        """
        if not self.master_calib_df.empty:
            return self.master_calib_df.groupby('KEY').size()
        else:
            return pd.Series(0, index=self.calibdb_list)

    @property
    def master_num_entry(self) -> pd.Series:
        """Number of 'master' entries in calibdb.

        :rtype: pd.Series
        """
        if not self.master_calib_df.empty:
            master_mask = self.master_calib_df.MASTER
            if master_mask.sum() == 0:
                return pd.Series(0, index=self.calibdb_list)

            else:
                master_calib_group = self.master_calib_df[master_mask].groupby('KEY')
                return master_calib_group.size()
        else:
            return pd.Series(0, index=self.calibdb_list)

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
        if not self.output_calibdb.empty:
            return self.output_calibdb.groupby('KEY').size()
        else:
            return pd.Series(0, index=self.calibdb_list)

    @property
    def calib_missing_mask(self) -> pd.Series:
        """calib_missing_mask.

        :rtype: pd.Series
        """
        return ~self.master_calib_df.FILE.isin(self.output_calibdb)

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
            files = files.explode()  # Pop lists in individual cells
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
                   fibers: Union[List[str], Tuple[str]] = ('AB', 'A', 'B', 'C'),
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

    def _load_log(self, logdirs: List[str]) -> Tuple[pd.DataFrame, List[str]]:
        """_load_log.

        Parse all log files and return a dataframe with only entries that have
        the current recipe in LOGFILE.

        :param logdirs:
        :type logdirs: List[str]
        :returns: dataframe and list of missing log files
        :rtype: Tuple[pd.DataFrame, list]
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

        # Use DIRECTORY as index and keep only relevant entries
        # whith self.recipe or extract recipes called
        log_df = log_df.set_index(['DIRECTORY'])

        # Important for extract recipes: keep only the recipe under test
        log_df = log_df[log_df.LOGFILE.str.contains(self.recipe)]

        # Paths without files
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
        try:
            master_calib_df = master_calib_df.loc[calib_inds]
        except KeyError:
            master_calib_df = master_calib_df[0:0]

        # Replace fibers by {FIBER} to match other tests
        # NOTE: Assumes fibers are at the end with \\b
        if self.fibers is not None:
            master_calib_df.index = master_calib_df.index.str.replace(
                '_(' + '|'.join(self.fibers) + ')\\b',
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
            files = files.explode()  # Pop lists in individual cells
        else:
            files = output_calib_paths.apply(CalibTest._fiberglob,
                                             fibers=self.fibers)
            for _ in self.fibers:
                files = files.explode()  # Nested lists so explode Nfiber times

        try:
            output_calibdb = files.apply(os.path.basename)
        except TypeError:
            output_calibdb = files[0:0]

        return output_calibdb

    @property
    def output_num_night(self) -> pd.Series:
        """Count outputs per night (and per output with multi-index).

        :rtype: pd.Series
        """
        ind_name = self.output_files.index.names
        return self.output_files.groupby(ind_name).size()

    # =========================================================================
    # Checks and stop common to all calibs
    # =========================================================================
    def check_qc(self, ncheck: int = 0) -> dict:
        """check_qc."""

        # Check passed QC and ENDED, can then access easily later
        num_logfits_qc_failed = len(self.log_all.qc_failed)

        # check: QC
        if num_logfits_qc_failed == 0:
            comments_check_qc = ''
            inspect_check_qc = ''
        else:
            comments_check_qc = 'One or more recipe have failed QC.'
            log_reset = self.log_all.qc_failed.reset_index()
            data_dict_check_qc = {
                'Night': log_reset.DIRECTORY.values,
                'QC_STRING': self.log_all.qc_failed.QC_STRING.values,
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
        num_logfits_ended_false = len(self.log_all.ended_false)

        # check: ENDED
        if num_logfits_ended_false == 0:
            comments_check_ended = ''
            inspect_check_ended = ''
        else:
            comments_check_ended = 'One or more recipe have failed to finish.'
            log_reset = self.log_all.ended_false.reset_index()
            data_dict_check_ended = {
                'Night': log_reset.DIRECTORY.values,
                'ERRORS': self.log_all.ended_false.ERRORS.values,
                'LOGFILE': self.log_all.ended_false.LOGFILE.values,
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
        qc_names = self.log_recipe.df.QC_NAMES.str.split(r'\|\|', expand=True)
        qc_names = qc_names.iloc[0]  # Keep only one row
        qc_values = self.log_recipe.df.QC_VALUES.str.split(r'\|\|', expand=True)
        qc_values.columns = qc_names
        # NOTE: .convert_dtypes will do in pd versions >= 1.0.0
        for col in qc_values.columns:
            try:
                qc_values[col] = qc_values[col].astype(float)
            except ValueError:
                del qc_values[col]

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

        # Define output duplicate mask
        (log_num_align,
         output_num_align) = self.log.get_align_count(self.output_num_night)
        output_dup_mask = log_num_align > output_num_align

        if output_dup_mask.any():
            # Check if '--master' in runstring of any file
            # Only for recipe log entries, not extract
            # These are a special case of duplicates
            if self.log_recipe.master_num > 0:
                comments_check_dup = ('An additional'
                                      f' {self.log_recipe.master_num}'
                                      ' recipe with --master=True was called'
                                      ' in the master directory'
                                      f' {self.log_recipe.master_nights}.'
                                      )
            else:
                comments_check_dup = ''

            # Get indices where outputs are duplicated in
            output_dup_ind = output_dup_mask[output_dup_mask]
            output_dup_ind = output_dup_ind.index.get_level_values('DIRECTORY')
            output_dup_ind = output_dup_ind.unique()

            # Check for non-master duplicates for each pattern
            dir_kwd = 'DIRECTORY'
            master_mask_group = self.log.master_mask.astype(int).groupby(dir_kwd)
            master_num_night = master_mask_group.sum()
            dup_num_night = master_mask_group.size()
            dup_not_master = dup_num_night - master_num_night
            log_dup_align, output_dup_align = dup_not_master.align(
                output_num_align[output_dup_mask]
            )
            true_dup = log_dup_align - output_dup_align  # number of dups
            true_dup = true_dup[true_dup > 0]  # Keep only non-zero
            true_dup.name = 'COUNT'
            true_dup = true_dup.reset_index()

            # TODO: implement something like this or just above w/o master/recipe
            #  dup = self.log_extract_num_align - self.output_num_align
            #  dup = dup[dup > 0]
            #  dup.name = 'COUNT'
            #  dup = dup.reset_index()

        else:
            comments_check_dup = ''
            true_dup = pd.DataFrame({'PATTERN': [],
                                     'DIRECTORY': [],
                                     'COUNT': []})

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
        # NOTE: Should we also check that the size is > 0 ?
        # Get total log count. For master get total else total - masters
        log_num = self.log.tot_num if self.ismaster else self.log.num

        # Condition for everything to pass
        passed = (self.output_num_unique == log_num).all()

        # Conditional pass
        cond = (self.output_num_unique < log_num).any()

        if passed:
            color = 'Lime'
            result = 'Yes'
            comment = ''
            inspect = ''
            data_dict = {}

        elif cond:
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
                    'Night': dup.DIRECTORY.values,
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
        if not self.output_files.empty:
            full_paths = (self.reduced_path
                          + os.path.sep
                          + self.output_files.index.get_level_values('DIRECTORY')
                          + os.path.sep
                          + self.output_files)
            headers = full_paths.loc[self.output_list[0]].apply(fits.getheader)
            header_df = pd.DataFrame(headers.tolist(), index=headers.index)

            # Keep only matching PIDs
            # NOTE: The output PID match the PID from cal_extract and not self.recipe
            # This will be corrected in 0.7
            log_pid_dir = self.log.df.reset_index().set_index('PID').DIRECTORY

            # make sure no duplicate PID per night
            log_pid_dir = log_pid_dir.reset_index().drop_duplicates().set_index(
                                                                        'PID'
                                                                        ).DIRECTORY
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
        else:
            missing_calib = pd.DataFrame([], columns=['LOC_DIR', 'FILE', 'CALIB_DIR'])

        return missing_calib

    def check_previous_calib(self, missing_calib, ncheck: int = 0) -> dict:
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
                                      ' used calibrations from'
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

    @property
    @abstractmethod
    def calls_extract(self) -> bool:
        """calls_extract.

        :rtype: bool
        """
