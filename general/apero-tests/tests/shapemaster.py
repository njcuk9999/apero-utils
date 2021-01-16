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
from typing import List, Optional, Union

import numpy as np
import pandas as pd

from .tests import CalibTest
from . import utils as ut


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

    @property
    def recipe_num_logfits(self) -> int:
        """recipe_num_logfits.

        For a master recipe, we don't need to subtract master entries.

        :rtype: int
        """
        return self.log_tot_num

    @property
    def output_num_entry(self) -> pd.Series:
        """Total number of entries per output in calibdb.

        For a master recipe, we don't need to subtract master entries.

        :rtype: pd.Series
        """
        return self.tot_num_entry

    def check_duplicates(self) -> pd.DataFrame:
        """check_duplicates.

        Duplicate defined as every night with more than one master_calib

        :rtype: pd.DataFrame
        """
        # NOTE: should change in APERO v0.7
        dup = self.output_num_align[self.output_num_align > 1]
        dup.name = 'COUNT'
        dup = dup.reset_index()

        return dup

    def check_qc_plot(self, ncheck: int = 0) -> dict:
        """check_qc_plot.

        Master recipe only has one night to plot so this will show QC per order

        :param ncheck:
        :type ncheck: int
        :rtype: dict
        """
        qc_names = self.log_df.QC_NAMES.str.split(r'\|\|', expand=True).iloc[0]
        qc_values = self.log_df.QC_VALUES.str.split(r'\|\|', expand=True)
        qc_values.columns = qc_names
        # NOTE: .convert_dtypes will do in pd versions >= 1.0.0
        float_mask = ~qc_values.isin(['True', 'False']).any()
        qc_values = qc_values.loc[:, float_mask].astype(float)

        # Only one value and one night here
        data_dict_check_qc_plot = {'Order': np.arange(1, 50),
                                   qc_names[0]: qc_values.values[0]}

        inspect_check_qc_plot = ut.inspect_plot(
                    self.test_id,
                    f'check{ncheck}',
                    data_dict_check_qc_plot,
                    f'{self.recipe}.py Quality Control',
                    )

        return inspect_check_qc_plot

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
                'reduced_path': self.reduced_path,
                'calibdb_path': self.calibdb_path,

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
