"""
Check if localisation calib worked fine.

Tests preformed
check1: how many recipes were run (cal_loc_{instrument} in log.fits)?
        how many in the master directory?
check2: how many of each output do we have?
        output1: {ODOMETER_CODE}_pp_order_profile_C.fits
        output2: {ODOMETER_CODE}_pp_loco_C.fits
        output3: {ODOMETER_CODE}_pp_fwhm-order_C.fits
        output4: {ODOMETER_CODE}_pp_with-order_C.fits
check3: how many of each unique output do we have?
        output1: {ODOMETER_CODE}_pp_order_profile_C.fits
        output2: {ODOMETER_CODE}_pp_loco_C.fits
        output3: {ODOMETER_CODE}_pp_fwhm-order_C.fits
        output4: {ODOMETER_CODE}_pp_with-order_C.fits
stop1: check3 == check1?
check4: using the log.fits how many files failed one or more QC?
        Which odometers? Which nights? Which QC?
check5: plot the different QCs as a function of time.
check6: using the log.fits how many files failed to finish? Which odometers?
        Which nights? Why (using the ERRORS and LOGFILE columns)?
check7: how many entry ORDER_PROFILE_{FIBER} and LOC_{FIBER} in
        master_calib_{INSTRUMENT}.txt?
check8: for each calib entry how many are in the calibDB?
stop2: check7 == check8?
check9: which bad pixel calibrations were used for each file? Was it the one
        from this night? How far away is the calibration obs time from the loc
        input file obs time?

@author: charles
"""
from typing import List, Optional, Union

import pandas as pd

from apero_tests import CalibTest
import apero_tests_func as atf


class LocTest(CalibTest):
    """LocTest."""

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
        return 'Localisation Recipe Test #1'

    @property
    def test_id(self) -> str:
        """test_id.

        :rtype: str
        """
        return 'localisation_test1'

    @property
    def output_list(self) -> List[str]:
        """List of output string patterns

        :return: output_list
        :rtype: list[str]
        """
        return ['*_pp_order_profile_{FIBER}.fits',
                '*_pp_loco_{FIBER}.fits',
                '*_pp_fwhm-order_{FIBER}.fits',
                '*_pp_with-order_{FIBER}.fits']

    @property
    def calibdb_list(self) -> List[str]:
        """List of calibDB entries

        :return: calibdb_list
        :rtype: list[str]
        """
        return ['ORDER_PROFILE_{FIBER}', 'LOC_{FIBER}']

    @property
    def recipe(self) -> List[str]:
        """Recipe name

        :return: output_list
        :rtype: list[str]
        """
        return 'cal_loc_{}'.format(self.instrument.lower())

    @property
    def fibers(self) -> List[str]:
        """fibers.

        :rtype: List[str]
        """
        return ['AB', 'C']

    # =========================================================================
    # Checks and stops
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
                    'Night': log_reset.DIRECTORY.tolist(),
                    'QC_STRING': self.log_qc_failed.QC_STRING,
                                }
            inspect_check_qc = atf.inspect_table(
                    'badpixel_test1',
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
                    'Night': log_reset.DIRECTORY.tolist(),
                    'ERRORS': self.log_ended_false.ERRORS,
                    'LOGFILE': self.log_ended_false.LOGFILE,
                                }
            inspect_check_ended = atf.inspect_table(
                    'badpixel_test1',
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

        inspect_check_qc_plot = atf.inspect_plot(
                    'localisation_test1',
                    f'check{ncheck}',
                    data_dict_check_qc_plot,
                    'cal_localisation_{0}.py Quality Control'.format(
                    self.instrument.lower())
                    )

        return inspect_check_qc_plot

    def stop_output_log(self, true_dup: pd.DataFrame, nstop: int = 0) -> dict:
        """stop_output_log.

        :param true_dup: df of true duplicate output

        :rtype: dict
        """

        # STOP 1: Unique output == log?
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
                        'Night': self.output_missing.DIRECTORY.tolist(),
                        'File name': self.output_missing.PATTERN.tolist(),
                         }
                inspect = atf.inspect_table(
                        'badpixel_test1',
                        f'stop{nstop}',
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
                        f'stop{nstop}',
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

        :param calib_dup: duplicates in calibdb
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

            # if missing output
            # NOTE: Could we have both missing and dup?
            if self.calib_missing_mask.any():

                comment = 'One or more outputs are not in the calibDB.'
                missing_mask = self.calib_missing_mask
                missing_calibdb_output = self.master_calib_df[missing_mask]
                data_dict = {
                        'Night': missing_calibdb_output.NIGHT.tolist(),
                        'File name': missing_calibdb_output.FILE.tolist(),
                        }
                inspect = atf.inspect_table(
                                'badpixel_test1',
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
                        'Night': calib_dup.NIGHT.tolist(),
                        'File name': calib_dup.FILE.tolist(),
                        'Occurrence': calib_dup.COUNT.tolist(),
                         }
                inspect = atf.inspect_table(
                        'badpixel_test1',
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

    def check_previous_calib(self) -> dict:
        """check_previous_calib.

        :rtype: dict
        """

        # TODO: check9
        # check9
        # NOTE: can probably get paths from output files df.
        # Then need to open fits one by one still with apply or loop through
        # files

        # check9
        # PID = logfits.PID[index_recipe]
        # for k in range(len(output1_list_files)):
        #     hdul = fits.open('{0}/{1}/{2}'.format(reduced_path,
        #             reduced_nights[i], output1_list_files[k]))
        #     if hdul[0].header['DRSPID'] in PID:
        #         if 'CDBBAD' in hdul[0].header:
        #             if not (hdul[0].header['CDBBAD'] == 'None' or 
        #                     os.path.isfile('{0}/{1}/{2}'.format(reduced_path,
        #                     reduced_nights[i], hdul[0].header['CDBBAD']))):
        #                 calibration_night_missing.append(reduced_nights[i])
        #                 calibration_missing.append(hdul[0].header['CDBBAD'])
        #         if 'CDBBACK' in hdul[0].header:      
        #             if not (hdul[0].header['CDBBACK'] == 'None' or 
        #                     os.path.isfile('{0}/{1}/{2}'.format(reduced_path,
        #                     reduced_nights[i], hdul[0].header['CDBBACK']))):
        #                 calibration_night_missing.append(reduced_nights[i])
        #                 calibration_missing.append(hdul[0].header['CDBBACK'])
        #     else:
        #         continue

        # check9
        # for i in range(len(calibration_missing)):
        #     path = glob.glob('{0}/*/{1}'.format(reduced_path,
        #             calibration_missing[i]))
        #     path = path[0].replace('{0}/'.format(reduced_path), '')
        #     calibration_night.append(path[:path.index('/')])


        # check9
        # if len(calibration_missing) == 0:
        #     comments_check9 = ''
        #     inspect_check9 = ''
        # else:
        #     comments_check9 = ('One or more localisation recipe outputs used '
        #                        'the bad pixel calibrations from another night')
        #     data_dict_check9 = {'Localisation Night':calibration_night_missing,
        #                         'Calibration file name': calibration_missing,
        #                         'Calibration Night': calibration_night
        #                         }
        #     inspect_check9 = atf.inspect_table(
        #                 'localisation_test1',
        #                 'check9',
        #                 data_dict_check9,
        #                 ('Night where the Localisation Recipe Output '
        #                  'used Calibrations from Another Night')
        #                 )

    def runtest(self):
        """runtest."""

        # Checking duplicates
        comments_check1, true_dup = self.check_duplicates()

        # QC/ENDED
        comments_check4, inspect_check4 = self.check_qc()
        inspect_check5 = self.check_qc_plot()
        comments_check6, inspect_check6 = self.check_ended()

        dict_stop1 = self.stop_output_log(true_dup)

        # Generate comment on extra master entries in calib db
        comments_check7 = ('An additional {0} {1} and {2} {3} with master = 1 '
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

        # TODO: check9 returns
        # dict_check9 = self.check_previous_calib()

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

                # check 5: QC Plot
                'inspect_check5': inspect_check5,

                # check 6: not ended
                'log_ended_false': self.log_ended_false,
                'comments_check6': comments_check6,
                'inspect_check6': inspect_check6,

                # Check 7: calibdb entries
                'output_num_entry': self.output_num_entry,
                'comments_check7': comments_check7,

                # Check 8: calibdb outputs
                'output_num_calibdb': self.output_num_calibdb,
                'output_dict': self.calib_output_dict,

                # Stop 2: calib output == entries
                'dict_stop2': dict_stop2,

                # TODO: Check9: previous calibrations (here and in template)
                }

        self.gen_html(html_dict)


if __name__ == '__main__':
    test = LocTest()
    test.runtest()
