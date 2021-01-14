"""
Check if shape master calib worked fine.

Tests preformed
check1: how many recipes were run (cal_shape_master_{instrument} in log.fits)?
check2: how many of each output do we have?
        output1: {ODOMETER_CODE}_pp_shapex.fits
        output2: {ODOMETER_CODE}_pp_shapey.fits
        output3: {ODOMETER_CODE}_pp_fpmaster.fits
check3: how many of each unique output do we have?
        output1: {ODOMETER_CODE}_pp_shapex.fits
        output2: {ODOMETER_CODE}_pp_shapey.fits
        output3: {ODOMETER_CODE}_pp_fpmaster.fits
stop1: check3 == check1?
check4: using the log.fits how many files failed one or more QC?
        Which odometers? Which nights? Which QC?
check5: plot the different QCs as a function of time.
check6: using the log.fits how many files failed to finish? Which odometers?
        Which nights? Why (using the ERRORS and LOGFILE columns)?
check7: how many entry SHAPEX, SHAPEY and FPMASTER in
        master_calib_{INSTRUMENT}.txt?
check8: for each calib entry how many are in the calibDB?
stop2: check8 == check7?
check9: which bad pixel/localisation calibrations were used for each file? 
        Was it the one from this night? How far away is the calibration obs 
        time from the shape master input file obs time?

@author: charles
"""
import os
from typing import List, Optional, Union

import glob
import pandas as pd
from astropy.io import fits

from apero_tests import CalibTest
import apero_tests_func as atf


class ShapeMTest(CalibTest):
    """ShapeMTest."""

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

    @property
    def name(self) -> str:
        """name."""
        return 'Shape Master Recipe Test #1'

    @property
    def test_id(self) -> str:
        """test_id.

        :rtype: str
        """
        return 'shapemaster_test1'

    @property
    def output_list(self) -> List[str]:
        """List of output string patterns

        :return: output_list
        :rtype: list[str]
        """
        return ['*_pp_shapex.fits',
                '*_pp_shapey.fits',
                '*_pp_fpmaster.fits']

    @property
    def calibdb_list(self) -> List[str]:
        """List of calibDB entries

        :return: calibdb_list
        :rtype: list[str]
        """
        return ['SHAPEX', 'SHAPEY', 'FPMASTER']

    @property
    def previous_calibs(self) -> List[str]:
        """previous_calibs.

        :rtype: List[str]
        """
        return ['CDBBAD', 'CDBBACK', 'CDBORDP', 'CDBLOCO']

    @property
    def recipe(self) -> List[str]:
        """Recipe name

        :return: output_list
        :rtype: list[str]
        """
        return 'cal_shape_master_{}'.format(self.instrument.lower())

    @property
    def fibers(self) -> List[str]:
        """fibers.

        :rtype: List[str]
        """

    # =========================================================================
    # Overwritten parent methods
    # =========================================================================
    @property
    def output_missing(self) -> pd.DataFrame:
        """Overwrite output_missing for master recipe

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

        Duplicate defined as every night with more than one master_calib

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
                    self.test_id,
                    f'check{ncheck}',
                    data_dict_check_qc,
                    'Nights that Failed Quality Control'
                    )

        return comments_check_qc, inspect_check_qc

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
                    self.test_id,
                    f'check{ncheck}',
                    data_dict_check_qc_plot,
                    f'{self.recipe}.py Quality Control',
                    )

        return inspect_check_qc_plot

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
                    self.test_id,
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
                       'smaller than or equal to the number of recipe called.')
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
        if (self.output_num_calibdb == self.tot_num_entry).all():
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

        if missing_calib.shape[0]:
            comments_missing_calib = ''
            inspect_missing_calib = ''
        else:
            comments_missing_calib = ('One or more localisation recipe outputs'
                                      ' used the bad pixel calibrations from'
                                      ' another night')
            data_dict_missing_calib = {
                    'Localisation Night': missing_calib.LOC_DIR.values,
                    'Calibration file name': missing_calib.FILE.values,
                    'Calibration Night': missing_calib.CALIB_DIR.values,
                    }
            inspect_missing_calib = atf.inspect_table(
                        'localisation_test1',
                        f'check{ncheck}',
                        data_dict_missing_calib,
                        ('Night where the Localisation Recipe Output '
                         'used Calibrations from Another Night')
                        )

        return comments_missing_calib, inspect_missing_calib

    # =========================================================================
    # Running the test
    # =========================================================================
    def runtest(self):
        """runtest."""

        dup = self.check_duplicates()

        # QC/ENDED
        comments_check4, inspect_check4 = self.check_qc(ncheck=4)
        inspect_check5 = self.check_qc_plot(ncheck=5)
        comments_check6, inspect_check6 = self.check_ended(ncheck=6)


        dict_stop1 = self.stop_output_log(dup, nstop=1)

        # Check for duplicates in calibdb
        calib_dup_mask = self.master_calib_df.duplicated(keep=False)
        # master_mask = self.master_calib_df.MASTER  # master recipe
        calib_dup = self.master_calib_df[calib_dup_mask]
        calib_dup = calib_dup.set_index('FILE', append=True)
        ind_names = calib_dup.index.names  # get file and key to count
        calib_dup['COUNT'] = calib_dup.groupby(ind_names).size()
        calib_dup = calib_dup.reset_index('FILE')
        calib_dup = calib_dup.drop_duplicates()

        dict_stop2 = self.stop_calibdb(calib_dup, nstop=2)

        # Check previous calibs to see if missing any
        missing_previous = self.get_missing_previous_calib()
        comments_check9, inspect_check9 = ShapeMTest.check_previous_calib(
                                                            missing_previous,
                                                            ncheck=9)

        html_dict = {
                # Summary header info
                'name': self.name,
                'setup': self.setup,
                'instrument': self.instrument,
                'date': self.date,
                'output_list': self.output_list,
                'calibdb_list': self.calibdb_list,
                'recipe': self.recipe,

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

                # Check 5: QC Plot
                'inspect_check5': inspect_check5,

                # Check 6: not ended
                'log_ended_false': self.log_ended_false,
                'comments_check6': comments_check6,
                'inspect_check6': inspect_check6,

                # Check 7: calibdb entries
                'output_num_entry': self.tot_num_entry,  # Master recipe
                # 'comments_check6': comments_check6,

                # Check 8: calibdb outputs
                'output_num_calibdb': self.output_num_calibdb,
                'output_dict': self.calib_output_dict,

                # Stop 2: calib output == entries
                'dict_stop2': dict_stop2,

                # Check9: previous calibrations
                'missing_previous': missing_previous,
                'comments_check9': comments_check9,
                'inspect_check9': inspect_check9,
                }

        self.gen_html(html_dict)

if __name__ == '__main__':
    test = ShapeMTest()
    test.runtest()
