"""
Check if wavelength master calib worked fine.

Tests preformed
check1: how many recipes were run (cal_wave_master_{instrument} in log.fits)?
        how many in the master directory?
check2: how many recipes were run (cal_extract_{instrument} in log.fits)?
check3: how many of each output do we have?
        output1: {ODOMETER_CODE}_pp_e2ds_{FIBER}.fits
        output2: {ODOMETER_CODE}_pp_e2dsff_{FIBER}.fits
        output3: {ODOMETER_CODE}_pp_s1d_w_{FIBER}.fits
        output4: {ODOMETER_CODE}_pp_s1d_v_{FIBER}.fits
        output5: {ODOMETER_CODE}_pp_e2dsff_linelist_{FIBER}.dat
        output6: {ODOMETER_CODE}_pp_e2dsff_wavemres_{FIBER}.fits
        output7: {ODOMETER_CODE}_pp_e2dsff_wavem_hc_{FIBER}.fits
        output8: {ODOMETER_CODE}_pp_e2dsff_wavem_fp_{FIBER}.fits
        output9: cal_wave_results.tbl
        output10: {ODOMETER_CODE}_pp_e2dsff_mhc_lines_{FIBER}.tbl
        output11: {ODOMETER_CODE}_pp_wavem_hclines_{FIBER}.fits 
        output12: {ODOMETER_CODE}_pp_wavem_fplines_{FIBER}.fits
        output13: {ODOMETER_CODE}_pp_e2dsff_ccf_{FIBER}.fits
check4: how many of each unique output do we have?
        output1: {ODOMETER_CODE}_pp_e2ds_{FIBER}.fits
        output2: {ODOMETER_CODE}_pp_e2dsff_{FIBER}.fits
        output3: {ODOMETER_CODE}_pp_s1d_w_{FIBER}.fits
        output4: {ODOMETER_CODE}_pp_s1d_v_{FIBER}.fits
        output5: {ODOMETER_CODE}_pp_e2dsff_linelist_{FIBER}.dat
        output6: {ODOMETER_CODE}_pp_e2dsff_wavemres_{FIBER}.fits
        output7: {ODOMETER_CODE}_pp_e2dsff_wavem_hc_{FIBER}.fits
        output8: {ODOMETER_CODE}_pp_e2dsff_wavem_fp_{FIBER}.fits
        output9: cal_wave_results.tbl
        output10: {ODOMETER_CODE}_pp_e2dsff_mhc_lines_{FIBER}.tbl
        output11: {ODOMETER_CODE}_pp_wavem_hclines_{FIBER}.fits 
        output12: {ODOMETER_CODE}_pp_wavem_fplines_{FIBER}.fits
        output13: {ODOMETER_CODE}_pp_e2dsff_ccf_{FIBER}.fits
stop1:  check4 == check2?
check5: using the log.fits how many files failed one or more QC?
        Which odometers? Which nights? Which QC?
check6: plot the different QCs as a function of time.
check7: using the log.fits how many files failed to finish? Which odometers?
        Which nights? Why (using the ERRORS and LOGFILE columns)?
check8: how many entry WAVEM_{FIBER}, WAVEHCL_{FIBER}, WAVEFPL_{FIBER}
        in master_calib_{INSTRUMENT}.txt?
check9: for each calib entry how many are in the calibDB?
stop2:  check9 == check8?
check10: which previous calibrations (bad pixel, loc, shape master, shape,
        flat/blaze) were used? Are they from this night? How far away is the
        calibration obs time from the loc input file obs time?

@author: charles
"""
from typing import List, Optional, Union

import numpy as np
import pandas as pd
from . import utils as ut

from .tests import CalibTest


class WaveMTest(CalibTest):
    """WaveMTest."""
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
        return 'Master Wavelength Solution Recipe Test #1'

    @property
    def test_id(self) -> str:
        """test_id.

        :rtype: str
        """
        return 'wavelengthmaster_test1'

    @property
    def output_list(self) -> List[str]:
        """List of output string patterns

        :return: output_list
        :rtype: list[str]
        """

        # Which outputs should we look? Some are .dat or .tbl.

        return ['*_pp_e2dsff_{FIBER}_wavem_fp_{FIBER}.fits',
                '*_pp_e2dsff_{FIBER}_wavem_hclines_{FIBER}.fits',
                '*_pp_e2dsff_{FIBER}_wavem_fplines_{FIBER}.fits',
                '*_pp_e2dsff_{FIBER}_wavemres_{FIBER}.fits',
                '*_pp_e2dsff_{FIBER}_wavem_hc_{FIBER}.fits']

    @property
    def calibdb_list(self) -> List[str]:
        """List of calibDB entries

        :return: calibdb_list
        :rtype: list[str]
        """
        return ['WAVEM_{FIBER}', 'WAVEHCL_{FIBER}', 'WAVEFPL_{FIBER}']

    @property
    def previous_calibs(self) -> List[str]:
        """previous_calibs.

        :rtype: List[str]
        """
        return [
            'CDBBAD', 'CDBBACK', 'CDBORDP', 'CDBLOCO', 'CDBSHAPX', 'CDBSHAPY',
            'CDBSHAPL', 'CDBFLAT', 'CDBBLAZE', 'CDBTHERM'
        ]

    @property
    def recipe(self) -> List[str]:
        """Recipe name

        :return: output_list
        :rtype: list[str]
        """
        return 'cal_wave_master_{}'.format(self.instrument.lower())

    @property
    def ismaster(self) -> bool:
        """Is the test for a master recipe.

        :rtype: bool
        """
        return True

    @property
    def fibers(self) -> List[str]:
        """fibers.

        :rtype: List[str]
        """
        return ['AB', 'A', 'B', 'C']

    @property
    def calls_extract(self) -> bool:
        """Does this method call extract.

        :rtype: bool
        """
        return True

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
        output_frame = self.output_files.index.to_frame().reset_index(
            drop=True)
        night_comb = [(p, n) for p in self.output_list
                      for n in self.reduced_nights]
        all_nights = pd.DataFrame(night_comb, columns=['PATTERN', 'DIRECTORY'])
        output_missing = all_nights[~all_nights.isin(output_frame).all(axis=1)]

        return output_missing

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
        dup = self.output_num_night[self.output_num_night > 1]
        dup.name = 'COUNT'
        dup = dup.reset_index()

        return dup

    def check_qc_plot(self, ncheck: int = 0) -> dict:
        """check_qc_plot.

        Only one night (master night). A table is more appropriate.

        :rtype: dict
        """
        qc_names = self.log_recipe.df.QC_LOGIC.str.split(r'\|\|', expand=True)
        qc_names = qc_names.iloc[0]  # Keep only one row
        qc_values = self.log_recipe.df.QC_VALUES.str.split(r'\|\|', expand=True)
        qc_values.columns = qc_names

        # NOTE: Duplicate QC_NAMES, why? For now, it overwrites.
        data_dict_check_qc_plot = {'Night': qc_values.index.tolist()}
        for key, series in qc_values.iteritems():
            data_dict_check_qc_plot[key] = series.tolist()

        inspect_check_qc_plot = ut.inspect_table(
            self.test_id,
            f'check{ncheck}',
            data_dict_check_qc_plot,
            f'{self.recipe}.py Quality Control'
        )

        return inspect_check_qc_plot

    # =========================================================================
    # Run the full test
    # =========================================================================
    def runtest(self):
        """runtest."""

        dup = self.check_duplicates()

        # QC/ENDED
        comments_check5, inspect_check5 = self.check_qc(ncheck=5)
        inspect_check6 = self.check_qc_plot(ncheck=6)
        comments_check7, inspect_check7 = self.check_ended(ncheck=7)

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

        # If the number of outputs is 0, comment on the default wave master files
        # used (MASTER_WAVE_{FIBER}.fits)
        if ((self.output_num_entry[0] + self.output_num_entry[1] +
                self.output_num_entry[2]) == 0):
            comments_check8 = ('4 WAVEM_D_{FIBER} entries not DRS processed. '
                    'Default MASTER_WAVE_2400416c_{FIBER}.fits files used.')
        else:
            comments_check8 = ''

        dict_stop2 = self.stop_calibdb(calib_dup, nstop=2)

        # Check previous calibs to see if missing any

        # No DRSPID in the header. Fixed by Neil, it will be corrected
        # in the next mini_data run.
        #missing_previous = self.get_missing_previous_calib()
        #comments_check10, inspect_check9 = self.check_previous_calib(
        #                                                    missing_previous,
        #                                                    ncheck=10)

        html_dict = {
            # Summary header info
            'name': self.name,
            'setup': self.setup,
            'instrument': self.instrument,
            'recipe': self.recipe,
            'date': self.date,
            'reduced_path': self.reduced_path,
            'output_list': self.output_list,
            'calibdb_list': self.calibdb_list,
            'calibdb_path': self.calibdb_path,

            # check 1 for logs
            'recipe_num_logfits': self.log_recipe.tot_num,

            # check 2 for logs
            'recipe_extract_num_logfits': self.log_extract.tot_num,

            # check 3 for outputs
            'output_num_total': self.output_num_total,

            # check 4
            'output_num_unique': self.output_num_unique,

            # stop 1: output==log
            'dict_stop1': dict_stop1,

            # check 5: QC failed
            'log_qc_failed': self.log_all.qc_failed,
            'comments_check5': comments_check5,
            'inspect_check5': inspect_check5,

            # check 6: QC Plot
            'inspect_check6': inspect_check6,

            # check 7: not ended
            'log_ended_false': self.log_all.ended_false,
            'comments_check7': comments_check7,
            'inspect_check7': inspect_check7,

            # Check 8: calibdb entries
            'output_num_entry': self.output_num_entry,
            'comments_check8': comments_check8,

            # Check 9: calibdb outputs
            'output_num_calibdb': self.output_num_calibdb,
            'output_dict': self.calib_output_dict,

            # Stop 2: calib output == entries
            'dict_stop2': dict_stop2,

            # Check10: previous calibrations
            'missing_previous': '',#missing_previous,
            'comments_check10': '',#comments_check10,
            'inspect_check10': '',#inspect_check10,
        }

        self.gen_html(html_dict)


if __name__ == '__main__':
    test = WaveMTest()
    test.runtest()
