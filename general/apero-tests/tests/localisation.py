"""
Check if localisation calib worked fine.

Tests preformed
check1: how many recipes were run (cal_loc_{instrument} in log.fits)?
        how many in the master directory?
check2: how many of each output do we have?
        output1: {ODOMETER_CODE}_pp_order_profile_{FIBER}.fits
        output2: {ODOMETER_CODE}_pp_loco_{FIBER}.fits
        output3: {ODOMETER_CODE}_pp_fwhm-order_{FIBER}.fits
        output4: {ODOMETER_CODE}_pp_with-order_{FIBER}.fits
check3: how many of each unique output do we have?
        output1: {ODOMETER_CODE}_pp_order_profile_{FIBER}.fits
        output2: {ODOMETER_CODE}_pp_loco_{FIBER}.fits
        output3: {ODOMETER_CODE}_pp_fwhm-order_{FIBER}.fits
        output4: {ODOMETER_CODE}_pp_with-order_{FIBER}.fits
stop1: check3 == check1?
check4: using the log.fits how many files failed one or more QC?
        Which odometers? Which nights? Which QC?
check5: plot the different QCs as a function of time.
check6: using the log.fits how many files failed to finish? Which odometers?
        Which nights? Why (using the ERRORS and LOGFILE columns)?
check7: how many entry ORDER_PROFILE_{FIBER} and LOC_{FIBER} in
        master_calib_{INSTRUMENT}.txt?
check8: for each calib entry how many are in the calibDB?
stop2: check8 == check7?
check9: which bad pixel calibrations were used for each file? Was it the one
        from this night? How far away is the calibration obs time from the loc
        input file obs time?

@author: charles
"""
from typing import List, Optional, Union

from .tests import CalibTest


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
    def previous_calibs(self) -> List[str]:
        """previous_calibs.

        :rtype: List[str]
        """
        return ['CDBBAD', 'CDBBACK']

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
    # Run the full test
    # =========================================================================
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

        # Check previous calibs to see if missing any
        missing_previous = self.get_missing_previous_calib()
        comments_check9, inspect_check9 = LocTest.check_previous_calib(
                                                            missing_previous,
                                                            ncheck=9)

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

                # Check9: previous calibrations
                'missing_previous': missing_previous,
                'comments_check9': comments_check9,
                'inspect_check9': inspect_check9,
                }

        self.gen_html(html_dict)


if __name__ == '__main__':
    test = LocTest()
    test.runtest()
