"""
Check if dark master calib worked fine.

Tests performed:
    check1: how many recipes were run? (cal_dark_master_{instrument} in
            other/log.fits)
    check2: how many of each output do we have?
            output1: {ODOMETER_CODE}_pp_dark_master.fits
    check3: how many of each unique output do we have?
            output1: {ODOMETER_CODE}_pp_dark_master.fits
    stop1: check3 == check1?
    check4: using the log.fits how many entry failed one or more QC?
            Which nights? Which QC?
    check5: using the log.fits how many entry failed to finish? Which nights?
            Why (using the ERRORS and LOGFILE columns)?
    check6: how many entry DARKM in master_calib_{INSTRUMENT}.txt?
    check7: for each calib entry how many are in the calibDB?
    stop2: check7 == check6?

@author: charles
"""
import os
from typing import Optional, Union, List

import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from apero_tests import Test, CalibTest
import apero_tests_func as atf


class DarkMTest(CalibTest):
    """DarkMTest."""

    def __init__(self,
                 inst: str = 'SPIROU',
                 setup: Optional[str] = None,
                 logdir: Union[str, list] = 'other'):
        """__init__.

        :param inst:
        :type inst: str
        :param setup:
        :type setup: Optional[str]
        :param logdir:
        :type logdir: Union[str, list]
        """
        super().__init__(inst=inst, setup=setup, logdir=logdir)

    # =========================================================================
    # Specific properties
    # =========================================================================
    @property
    def name(self) -> str:
        """name.

        :rtype: str
        """
        return 'Dark Master Recipe Test #1'

    @property
    def test_id(self) -> str:
        """test_id.

        :rtype: str
        """
        return 'darkmaster_test1'

    @property
    def output_list(self) -> List[str]:
        """List of output string patterns

        :return: output_list
        :rtype: list[str]
        """
        return ['*_pp_dark_master.fits']

    @property
    def calibdb_list(self) -> List[str]:
        """List of calibDB entries

        :return: calibdb_list
        :rtype: list[str]
        """
        return ['DARKM']

    @property
    def recipe(self) -> List[str]:
        """Recipe name

        :return: output_list
        :rtype: list[str]
        """
        return 'cal_dark_master_{self.instrument.lower()}'

    @property
    def fibers(self) -> None:
        """fibers.
        No fibers for darkmaster
        """

    # =========================================================================
    # Overwritten parent methods
    # =========================================================================
    @property
    def output_missing(self) -> pd.DataFrame:
        """Overwrite output_missing for specific case of dark master.

        Missing outputs defined as every night without a dark master calib

        :rtype: pd.DataFrame
        """
        # NOTE: should change in APERO v0.7
        output_frame = self.output_files.index.to_frame().reset_index(drop=True)
        night_comb = [(p, n)
                      for p in self.output_list
                      for n in self.reduced_nights]
        all_nights = pd.DataFrame(night_comb, columns=['PATTERN', 'DIRECTORY'])
        output_missing = all_nights[~all_nights.isin(output_frame).all(axis=1)]

        return output_missing

    # =========================================================================
    # Checks and stops
    # =========================================================================
    def check_duplicates(self) -> pd.DataFrame:
        """check_duplicates.

        Duplicate defined as every night with more than one dark master calib

        :rtype: pd.DataFrame
        """
        # NOTE: should change in APERO v0.7
        dup = self.output_night_num[self.output_night_num > 1]
        dup.name = 'COUNT'
        dup = dup.reset_index()

        return dup

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
            inspect_check_qc = atf.inspect_table(
                    'darkmaster_test1',
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
            inspect_check_ended = atf.inspect_table(
                    'darkmaster_test1',
                    f'check{ncheck}',
                    data_dict_check_ended,
                    'Nights that Failed to Finish'
                    )

        return comments_check_ended, inspect_check_ended

    def stop_output_log(self, dup: pd.DataFrame, nstop: int = 0) -> dict:
        """stop_output_log.

        :param dup: df of true duplicate output

        :rtype: dict
        """

        if (self.output_num_unique == self.log_tot_num).all():
            color = 'Lime'
            result = 'Yes'
            comment = ''
            inspect = ''
            data_dict = {}

        elif (self.output_num_unique < self.log_tot_num).any():

            color = 'Yellow'
            result = 'No'
            data_dict = {}

        else:
            color = 'Red'
            result = 'No'
            comment = ('The number of unique output files should always be '
                       'smaller than or equal to the number of recipe called.'
                       )
            inspect = ''
            data_dict = {}

        # if missing output
        # NOTE: moved out because can have missing or duplicates in a night
        # even if match overall
        # NOTE: could there be both missing output and duplicate?
        if self.output_missing.shape[0] > 0:
            comment = 'One or more output were not produced.'
            data_dict = {'Night': self.output_missing.DIRECTORY.values,
                         'File name': self.output_missing.PATTERN.values,
                         }
            inspect = atf.inspect_table('darkmaster_test1',
                                        f'stop{nstop}',
                                        data_dict,
                                        'Missing Outputs in {0}'.format(
                                            self.reduced_path
                                            ),
                                        )
        # if duplicates
        else:
            comment = ('Recipe called multiple times producing the same '
                       'outputs.')
            data_dict = {'Night': dup.DIRECTORY.values,
                         'File name': dup.PATTERN.values,
                         'Occurrence': dup.COUNT.values,
                         }
            inspect = atf.inspect_table('darkmaster_test1',
                            f'stop{nstop}',
                            data_dict,
                            ('Same Dark Master Recipe Called Twice or More '
                             'Producing the Same Outputs in {0}'
                             ).format(self.reduced_path)
                            )

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

        elif (self.output_num_calibdb < self.tot_num_entry).any():

            color = 'Yellow'
            result = 'No'


            # if missing output
            # NOTE: Could we have both missing and dup?
            if self.calib_missing_mask.any():

                comment = 'One or more outputs are not in the calibDB.'
                missing_mask = self.calib_missing_mask
                missing_calibdb_output = self.master_calib_df[missing_mask]
                data_dict = {
                        'Night': missing_calibdb_output.NIGHT.values,
                        'File name': missing_calibdb_output.FILE.values,
                        }
                inspect = atf.inspect_table(
                                'darkmaster_test1',
                                f'stop{nstop}',
                                data_dict,
                                f'Missing Output in {self.calibDB_path}'
                                )
            # if duplicates
            else:

                comment = ('Some entries in'
                           f'master_calib_{self.instrument}.txt'
                           'are identical.')
                data_dict = {
                        'Night': calib_dup.NIGHT.values,
                        'File name': calib_dup.FILE.values,
                        'Occurrence': calib_dup.COUNT.values,
                         }
                inspect = atf.inspect_table(
                        'darkmaster_test1',
                        f'stop{nstop}',
                        data_dict,
                        ('Duplicate Entries in '
                         'master_calib_{0}.txt').format(self.instrument)
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

    # =========================================================================
    # Running the test
    # =========================================================================
    def runtest(self):
        """runtest."""

        # =====================================================================
        # Inspect outputs
        # =====================================================================
        # Checking for duplicates
        dup = self.check_duplicates()

        # QC/ENDED
        comments_check4, inspect_check4 = self.check_qc()
        comments_check5, inspect_check5 = self.check_ended()

        dict_stop1 = self.stop_output_log(dup)

        # =====================================================================
        # Inspect calibDB
        # =====================================================================
        # Check if some files have duplicaltes in db (empty if none)
        # Use tot_num_entry for dark MASTER recipe (i.e. not remove masters)
        calib_dup_mask = self.master_calib_df.duplicated(keep=False)
        master_mask = self.master_calib_df.MASTER
        calib_dup = self.master_calib_df[calib_dup_mask & ~master_mask]
        calib_dup = calib_dup.set_index('FILE', append=True)
        ind_names = calib_dup.index.names  # get file and key to count
        calib_dup['COUNT'] = calib_dup.groupby(ind_names).size()
        calib_dup = calib_dup.reset_index('FILE')
        calib_dup = calib_dup.drop_duplicates()

        dict_stop2 = self.stop_calibdb(calib_dup)

        html_dict = {
                # Summary header info
                'setup': self.setup,
                'instrument': self.instrument,
                'date': self.date,
                'output_list': self.output_list,
                'calibdb_list': self.calibdb_list,

                # Check 1: number of calls in logfits
                'recipe_num_logfits': self.log_tot_num,  # Master recipe

                # Check 2 number of outputs
                'output_num_total': self.output_num_total,

                # Check 3 number of unique outputs
                'output_num_unique': self.output_num_unique,

                # Stop 1
                'dict_stop1': dict_stop1,

                # Check 4: QC failed
                'log_qc_failed': self.log_qc_failed,
                'comments_check4': comments_check4,
                'inspect_check4': inspect_check4,

                # check 5: not ended
                'log_ended_false': self.log_ended_false,
                'comments_check5': comments_check5,
                'inspect_check5': inspect_check5,

                # Check 6: calibdb entries
                'output_num_entry': self.tot_num_entry,  # Master recipe
                # 'comments_check6': comments_check6,

                # Check 7: calibdb outputs
                'output_num_calibdb': self.output_num_calibdb,
                'output_dict': self.calib_output_dict,

                # Stop 2: calib output == entries
                'dict_stop2': dict_stop2,
                }
        
        self.gen_html(html_dict)


if __name__ == '__main__':
    test = DarkMTest()
    test.runtest()
