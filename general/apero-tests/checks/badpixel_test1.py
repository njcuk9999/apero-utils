"""
Check if bad pixel calib worked fine.

Tests preformed
    check1: how many recipes were run (cal_badpix_{instrument} in log.fits)?
            how many in the master directory?
    check2: how many of each output do we have?
            #output1: {ODOMETER_CODE}_pp_badpixel.fits
            #output2: {ODOMETER_CODE}_pp_bmap.fits
    check3: how many of each unique output do we have?
            #output1: {ODOMETER_CODE}_pp_badpixel.fits
            #output2: {ODOMETER_CODE}_pp_bmap.fits
    stop1: check3 == check1?
    check4: using the log.fits how many entry failed one or more QC?
            Which nights? Which QC?
    check5: using the log.fits how many entry failed to finish? Which nights?
            Why (using the ERRORS and LOGFILE columns)?
    check6: how many entry BADPIX and BKGRDMAP in
            master_calib_{INSTRUMENT}.txt?
    check7: for each calib entry how many are in the calibDB?
    stop2: check7 == check6?

@author: charles
"""
from typing import List, Optional, Union

from apero_tests import CalibTest
import apero_tests_func as atf


class BadPixTest(CalibTest):
    """BadPixTest."""

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
        return 'Bad Pixel Correction Recipe Test #1'

    @property
    def test_id(self) -> str:
        """test_id.

        :rtype: str
        """
        return 'badpixel_test1'

    @property
    def output_list(self) -> List[str]:
        """List of output string patterns

        :return: output_list
        :rtype: list[str]
        """
        return ['*_pp_badpixel.fits', '*_pp_bmap.fits']

    @property
    def calibdb_list(self) -> List[str]:
        """List of calibDB entries

        :return: calibdb_list
        :rtype: list[str]
        """
        return ['BADPIX', 'BKGRDMAP']

    @property
    def recipe(self) -> List[str]:
        """Recipe name

        :return: output_list
        :rtype: list[str]
        """
        return 'cal_badpix_{}'.format(self.instrument.lower())


    @property
    def fibers(self) -> None:
        """fibers.
        No fibers for badpix
        """

    # =========================================================================
    # Checks and stops
    # =========================================================================
    def check_duplicates(self):
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
                comments_check1 = ('An additional {0} recipe with '
                        '--master=True was called in the master '
                        'directory {1}.'.format(self.master_recipe_num_logfits,
                                                self.master_nights)
                        )

            # Check for non-master duplicates for each pattern
            # NOTE: This replaces [night_]output{N}_dup and
            master_mask_group = self.master_mask.astype(int).groupby('DIRECTORY')
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
            comments_check1 = ''

        return comments_check1, true_dup

    def check_qc(self):
        """check_qc."""

        # Check passed QC and ENDED, can then access easily later
        num_logfits_qc_failed = self.log_qc_failed.shape[0]

        # check4: QC
        if num_logfits_qc_failed == 0:
            comments_check4 = ''
            inspect_check4 = ''
        else:
            comments_check4 = 'One or more recipe have failed QC.'
            log_reset = self.log_qc_failed.reset_index()
            data_dict_check4 = {
                    'Night': log_reset.DIRECTORY.tolist(),
                    'QC_STRING': self.log_qc_failed.QC_STRING,
                                }
            inspect_check4 = atf.inspect_table(
                    'badpixel_test1',
                    'check4',
                    data_dict_check4,
                    'Nights that Failed Quality Control'
                    )

        return comments_check4, inspect_check4

    def check_ended(self):
        """check_ended."""
        num_logfits_ended_false = self.log_ended_false.shape[0]

        # check5: ENDED
        if num_logfits_ended_false == 0:
            comments_check5 = ''
            inspect_check5 = ''
        else:
            comments_check5 = 'One or more recipe have failed to finish.'
            log_reset = self.log_ended_false.reset_index()
            data_dict_check5 = {
                    'Night': log_reset.DIRECTORY.tolist(),
                    'ERRORS': self.log_ended_false.ERRORS,
                    'LOGFILE': self.log_ended_false.LOGFILE,
                                }
            inspect_check5 = atf.inspect_table(
                    'badpixel_test1',
                    'check5',
                    data_dict_check5,
                    'Nights that Failed to Finish'
                    )

        return comments_check5, inspect_check5

    def stop_output_log(self, true_dup) -> dict:
        """stop_output_log.

        :param true_dup: df of true duplicate output

        :rtype: dict
        """

        # STOP 1: Unique output == log?
        if (self.output_unique_num == self.recipe_num_logfits).all():
            color = 'Lime'
            result = 'Yes'
            comment = ''
            inspect = ''
        elif (self.output_unique_num < self.recipe_num_logfits).any():

            color = 'Yellow'
            result = 'No'
            # if missing output
            # NOTE: could there be both missing output and duplicate?
            if self.output_missing.shape[0] > 0:

                comment = 'One or more output were not produced.'
                # NOTE: might need arrays instead of list
                data_dict = {
                        'Night': self.output_missing.DIRECTORY.tolist(),
                        'File name': self.output_missing.PATTERN.tolist(),
                         }
                inspect = atf.inspect_table(
                        'badpixel_test1',
                        'stop1',
                        data_dict,
                        'Missing Outputs in {0}'.format(self.reduced_path)
                        )
            # if duplicates
            else:
                comment = ('Recipe called multiple times producing the '
                                 'same outputs.')
                data_dict = {
                        'Night': self.true_dup.DIRECTORY.tolist(),
                        'File name': true_dup.PATTERN.tolist(),
                        'Occurrence': true_dup.COUNT.tolist(),
                                   }

                inspect_msg = ('Same Bad Pixel Correction Recipe Called Twice '
                               'or More Producing the Same Outputs in {0}'
                               ).format(self.reduced_path)
                inspect = atf.inspect_table(
                        'badpixel_test1',
                        'stop1',
                        data_dict,
                        inspect_msg
                        )

        else:
            color = 'Red'
            result = 'No'
            comment = ('The number of unique output files should always '
                       'be smaller or equal to the number of recipe '
                       'called.')
            inspect = ''

        stop_dict = {
                'data': data_dict,
                'comment': comment,
                'inspect': inspect,
                'color': color,
                'result': result,
                }

        return stop_dict


    def stop_calibdb(self, calib_dup) -> dict:
        """stop_calibdb.

        :param calib_dup: duplicates in calibdb
        :rtype: dict
        """
        # stop2
        if (self.output_num_calibdb == self.output_num_entry).all():
            color = 'Lime'
            result = 'Yes'
            comment = ''
            inspect = ''

        elif (self.output_num_calibdb < self.output_num_entry).any():

            color = 'Yellow'
            result = 'No'

            # if missing output
            # NOTE: Could we have both missing and dup?
            if self.calib_missing_mask.any():

                comment = 'One or more outputs are not in the calibDB.'
                missing_calibdb_output = self.master_calib_df[self.calib_missing_mask]
                data_dict = {
                        'Night': missing_calibdb_output.NIGHT.tolist(),
                        'File name': missing_calibdb_output.FILE.tolist(),
                        }
                inspect = atf.inspect_table(
                                'badpixel_test1',
                                'stop2',
                                data_dict,
                                f'Missing Output in {self.calibDB_path}'
                                )
            # if duplicates
            else:

                comment = ('Some entries in master_calib_{0}.txt are '
                                 'identical.').format(self.instrument)
                data_dict = {
                        'Night': calib_dup.NIGHT.tolist(),
                        'File name': calib_dup.FILE.tolist(),
                        'Occurrence': calib_dup.COUNT.tolist(),
                         }
                inspect = atf.inspect_table(
                        'badpixel_test1',
                        'stop2',
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

        stop_dict = {
                'data': data_dict,
                'comment': comment,
                'inspect': inspect,
                'color': color,
                'result': result,
                }

        return stop_dict



    def runtest(self):
        """runtest."""

        # Checking duplicates
        comments_check1, true_dup = self.check_duplicates()

        # QC/ENDED
        comments_check4, inspect_check4 = self.check_qc()
        comments_check5, inspect_check5 = self.check_ended()

        dict_stop1 = self.stop_output_log(true_dup)

        # =====================================================================
        # Inspect calibDB
        # =====================================================================

        # Generate comment on extra master entries
        comments_check6 = ('An additional {0} {1} and {2} {3} with master = 1 '
                           'are in the master_calibDB_{4}.'.format(
                               self.master_num_entry[self.calibdb_list[0]],
                               self.calibdb_list[0],
                               self.master_num_entry[self.calibdb_list[1]],
                               self.calibdb_list[1],
                               self.instrument)
                           )


        # Check if some files have duplicaltes in db
        if (self.output_num_calibdb < self.output_num_entry).any():

            # Get all duplicates
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

                # check 1 for logs
                'recipe_num_logfits': self.recipe_num_logfits,
                'comments_check1': comments_check1,

                # check 2 for outputsl
                'output_num_total': self.output_num_total,

                # check 3
                'output_num_unique': self.output_num_unique,

                # stop 1: output==log
                'dict_stop1': dict_stop1,

                # check 4: QC failed
                'log_qc_failed': self.log_qc_failed,
                'comments_check4': comments_check4,
                'inspect_check4': inspect_check4,

                # check 5: not ended
                'log_ended_false': self.log_ended_false,
                'comments_check5': comments_check5,
                'inspect_check5': inspect_check5,

                # Check 6: calibdb entries
                'output_num_entry': self.output_num_entry,

                # Check 7: calibdb outputs
                'output_num_calibdb': self.output_num_calibdb,

                # Stop 2: calib output == entries
                'dict_stop2': dict_stop2,
                }

        self.gen_html(html_dict)

if __name__ == '__main__':
    test = BadPixTest()
    test.runtest()
