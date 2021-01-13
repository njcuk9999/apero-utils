"""
Check if preprocessing worked fine.

Tests performed:
    check1: how many raw files are there on disk?
    check2: how many pp files are there on disk?
    stop1: check2 == check1?
    check3: how many pp files are there in the index.fits?
    stop2: check3 == check2?
    check4: how many recipes were run? (cal_preprocess_{instrument} in
            tmp/*/log.fits)
    check5: how many unique odometers files were preprocessed according to the
            log.fits?
    stop3: check5 == check 2?
    check6: using the log.fits how many unique odometers failed one or more QC?
            Which odometers? Which QC?
    check7: using the log.fits how many unique odometers failed to finish?
            Which odometers? Why (using the ERRORS and LOGFILE columns)?

@author: charles
"""
import os
from datetime import datetime
from typing import Optional, List, Tuple

import numpy as np
import pandas as pd
from astropy.table import Table

from apero_tests import Test
import apero_tests_func as atf
from apero.core import constants


class PPTest(Test):
    """PPTEST."""

    # pylint: disable=too-many-public-methods
    # A few more (20+5=5) here to keep same logic as other tests even though
    # inherits directly from Test (no intermediate in between)


    def __init__(self, inst: str = 'SPIROU', setup: Optional[str] = None):
        """__init__.

        :param inst: instrument used
        :type inst: str
        :param setup: APERO setup
        :type setup: Optional[str]
        """
        super().__init__(inst=inst, setup=setup)

        # list raw night directories
        self._raw_nights = atf.list_nights(self.raw_path)

        # list PP data night directories
        self._pp_nights = atf.list_nights(self.pp_path)

        # Counts directly from glob (done here to avoid parsing each time)
        self._raw_num = atf.count_files_subdir(self.raw_path,
                                               subdir='all',
                                               files='*.fits')
        self._pp_num = atf.count_files_subdir(self.pp_path,
                                              subdir='all',
                                              files=self.output_list[0])


        # NOTE: Load index and log, discard missing with _ for now (not used)
        dirs = self.pp_nights.copy()
        self._log_df, _ = self._load_log_index('log', dirs)
        self._index_df, _ = self._load_log_index('index', dirs)


    # =========================================================================
    # Abstract properties from parent Test
    # =========================================================================
    @property
    def name(self) -> str:
        """Test full unique name."""
        return 'Preprocessing Recipe Test #1'


    @property
    def test_id(self) -> str:
        """Test short name (ID)."""
        return 'preprocessing_test1'

    @property
    def recipe(self) -> str:
        """recipe.

        :rtype: str
        """
        return f'cal_preprocess_{self.instrument.lower()}'

    @property
    def fibers(self):
        """Fibers to search in outputs
        No fibers for PP
        """

    @property
    def output_list(self) -> List[str]:
        """List of output string patterns

        :return: output_list
        :rtype: List[str]
        """
        return ['*_pp.fits']

    # =========================================================================
    # Properties specific to PP
    # =========================================================================
    @property
    def raw_path(self) -> str:
        """Path to raw directory

        :return: reduced_path
        :rtype: str
        """
        return self.params['DRS_DATA_RAW'].rstrip(os.path.sep)

    @property
    def pp_path(self) -> str:
        """Path to PP directory

        :return: reduced_path
        :rtype: str
        """
        return self.params['DRS_DATA_WORKING'].rstrip(os.path.sep)

    @property
    def raw_nights(self) -> List[str]:
        """raw_nights.

        :rtype: List[str]
        """
        return self._raw_nights

    @property
    def pp_nights(self) -> List[str]:
        """pp_nights.

        :rtype: List[str]
        """
        return self._pp_nights

    @property
    def log_df(self) -> pd.DataFrame:
        """Dataframe with log information

        :return: log_df
        :rtype: pd.DataFrame
        """
        return self._log_df

    @property
    def index_df(self) -> pd.DataFrame:
        """Dataframe with index information

        :return: index_df
        :rtype: pd.DataFrame
        """
        return self._index_df

    @property
    def log_df_unique(self) -> pd.DataFrame:
        """log_df_unique.

        :rtype: pd.DataFrame
        """
        return self.log_df.drop_duplicates(subset='ODOMETER')

    # =========================================================================
    # Properties derived from outputs
    # =========================================================================
    @property
    def raw_num(self) -> int:
        """raw_num.

        :rtype: int
        """
        return self._raw_num

    @property
    def pp_num(self) -> int:
        """pp_num.

        :rtype: int
        """
        return self._pp_num

    @property
    def pp_num_index(self) -> int:
        """pp_num_indexfits.

        :rtype: int
        """
        return self.index_df.shape[0]

    @property
    def pp_num_log(self) -> int:
        """pp_num_log.

        :rtype: int
        """
        return self.log_df.shape[0]

    @property
    def pp_num_log_unique(self) -> int:
        """pp_num_log_unique.

        :rtype: int
        """
        return self.log_df_unique.shape[0]

    @property
    def log_qc_failed(self) -> pd.DataFrame:
        """log_qc_failed.

        :rtype: pd.DataFrame
        """
        return self.log_df_unique[~self.log_df_unique.PASSED_ALL_QC]

    @property
    def log_ended_false(self) -> pd.DataFrame:
        """log_ended_false.

        :rtype: pd.DataFrame
        """
        return self.log_df_unique[~self.log_df_unique.ENDED]

    # =========================================================================
    # Utility methods
    # =========================================================================
    def _load_log_index(self,
                        ftype: str,
                        dirs: List[str]) -> Tuple[pd.DataFrame, list]:
        """_load_log_index.
        Generate a dataframe with log for each directory, with DIRECTORY as
        index. Only logs of the recipe being tested are kept.

        :param ftype: file type (log)
        :type ftype: str
        :param dirs:
        :type dirs: List[str]
        :rtype: Tuple[pd.DataFrame, list]
        """
        if ftype not in ['log', 'index']:
            raise ValueError('ftype should be log or index')

        # Get all logs in a dataframe
        allpaths = [os.path.join(self.pp_path, d, f'{ftype}.fits')
                    for d in dirs]
        paths = [p for p in allpaths if os.path.isfile(p)]
        dfs = [Table.read(f).to_pandas() for f in paths]
        df = pd.concat(dfs)

        # Decode strings (fits from astropy are encoded)
        for col in df.columns:
            # try/except for non-string columns (no .str attribute)
            try:
                df[col] = df[col].str.decode('utf-8')
            except AttributeError:
                pass

        # Use DIRECTORY/NIGHTNAME as index (all entries should be PP)
        night_label = 'DIRECTORY' if ftype == 'log' else 'NIGHTNAME'
        df = df.set_index([night_label])

        # Store missing logs in another list
        missing_logs = [p for p in allpaths if not os.path.isfile(p)]

        # Some extra formatting for logs
        if ftype == 'log':
            df['ARGS'] = df.ARGS.str.replace('persi_', '')
            df['ODOMETER'] = df.ARGS.str.split('.fits').str[0].str[-8:]

        return df, missing_logs

    # =========================================================================
    # Checks and stops
    # =========================================================================
    def stop_raw_pp(self) -> dict:
        """stop_raw_pp.

        :rtype: dict
        """
        # Stop 1
        if self.pp_num == self.raw_num:
            color = 'Lime'
            result = 'Yes'
            comment = ''
            inspect = ''
        elif self.pp_num < self.raw_num:
            color = 'Yellow'
            result = 'No'
            comment = 'Not all available raw files were reduced.'
            inspect = ''
        else:
            color = 'Red'
            result = 'No'
            comment = ('The number of pp files should always be '
                       'smaller than the number of raw files.')
            inspect = ''

        stop_dict = {
                'comment': comment,
                'inspect': inspect,
                'color': color,
                'result': result,
                }

        return stop_dict

    def stop_index_pp(self) -> dict:
        """stop_index_pp.

        :rtype: dict
        """
        # Stop2
        if self.pp_num_index == self.pp_num:
            color = 'Lime'
            result = 'Yes'
            comment = ''
            inspect = ''
        else:
            color = 'Red'
            result = 'No'
            comment = ''
            inspect = ''

        stop_dict = {
                'comment': comment,
                'inspect': inspect,
                'color': color,
                'result': result,
                }

        return stop_dict

    def stop_log_pp(self) -> dict:
        """stop_log_pp.

        :rtype: dict
        """
        # Stop 3
        if self.pp_num_log_unique == self.pp_num:
            color = 'Lime'
            result = 'Yes'
            comment = ''
            inspect = ''
        elif self.pp_num_log_unique > self.pp_num:
            color = 'Yellow'
            result = 'No'
            comment = 'Some files were processed more than once.'
            inspect = ''
        else:
            color = 'Red'
            result = 'No'
            comment = ''
            inspect = ''

        stop_dict = {
                'comment': comment,
                'inspect': inspect,
                'color': color,
                'result': result,
                }

        return stop_dict

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
                    'Odometer': self.log_qc_failed.ODOMETER.values,
                    'QC_STRING': self.log_qc_failed.QC_STRING.values,
                    }
            inspect_check_qc = atf.inspect_table(
                    'preprocessing_test1',
                    f'check{ncheck}',
                    data_dict_check_qc,
                    'Odometers that Failed One or More Quality Control'
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
                    'Odometer': self.log_ended_false.ODOMETER.values,
                    'ERRORS': self.log_ended_false.ERRORS.values,
                    'LOGFILE': self.log_ended_false.LOGFILE.values,
                                }
            inspect_check_ended = atf.inspect_table(
                    'badpixel_test1',
                    f'check{ncheck}',
                    data_dict_check_ended,
                    'Odometers that Failed to Finish'
                    )

        return comments_check_ended, inspect_check_ended

    # =========================================================================
    # Running the test
    # =========================================================================
    def runtest(self):
        """runtest."""

        # Stop 1: raw == pp?
        dict_stop1 = self.stop_raw_pp()

        # Stops for index and log == PP ?
        dict_stop2 = self.stop_index_pp()
        dict_stop3 = self.stop_log_pp()

        # QC/ENDED
        comments_check6, inspect_check6 = self.check_qc(ncheck=6)
        comments_check7, inspect_check7 = self.check_ended(ncheck=7)

        html_dict = {
                # Summary header info
                'setup:': self.setup,
                'instrument': self.instrument,
                'data': self.date,
                'output_list': self.output_list,

                # Check 1: raw
                'raw_path': self.raw_path,
                'raw_num': self.raw_num,

                # Check 2: PP
                'pp_path': self.pp_path,
                'pp_num': self.pp_num,

                # Stop 1: PP == raw
                'dict_stop1': dict_stop1,

                # Check 3: number of entries index.fits files
                'pp_num_index': self.pp_num_index,

                # Stop 2: PP == index
                'dict_stop2': dict_stop2,

                # Check 4: number of calls in logs
                'pp_num_log': self.pp_num_log,

                # Check 5: number of unique odometers in logs
                'pp_num_log_unique': self.pp_num_log_unique,

                # Stop 3: log_unique == PP
                'dict_stop3': dict_stop3,

                # Check 6: QC failed
                'pp_num_qc_failed': self.log_qc_failed.shape[0],
                'comments_check6': comments_check6,
                'inspect_check6': inspect_check6,

                # Check 6: QC failed
                'pp_num_ended_false': self.log_ended_false.shape[0],
                'comments_check7': comments_check7,
                'inspect_check7': inspect_check7,


                }

        self.gen_html(html_dict)


if __name__ == '__main__':
    test = PPTest()
    test.runtest()
