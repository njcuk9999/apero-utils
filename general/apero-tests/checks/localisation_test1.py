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
stop2: check7 == check6?
check9: which bad pixel calibrations were used for each file? Was it the one
        from this night? How far away is the calibration obs time from the loc
        input file obs time?

@author: charles
"""
from typing import List, Optional, Union

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

    def runtest(self):
        """runtest."""

        # Checking duplicates
        comments_check1, true_dup = self.check_duplicates()


        comments_check6 = ('An additional {0} {1} and {2} {3} with master = 1 '
                           'are in the master_calibDB_{4}.'.format(
                               master_output1_num_entry,
                               calibDB_entry_list[0], master_output2_num_entry,
                               calibDB_entry_list[1], instrument)
                           )

        #check5
        # TODO: implement this
        split_QCnames = np.array(QCnames[0].split('||'))
        split_QCvalues = np.array(
                [[x for x in e.split('||')] for e in QCvalues])
        index_float = ~np.logical_or(split_QCvalues[0] == 'True', 
                split_QCvalues[0] == 'False')
    
        split_QCnames = split_QCnames[index_float]
        split_QCvalues = split_QCvalues[:,index_float]
        split_QCvalues = split_QCvalues.astype(dtype=float)

        data_dict_check5 = {'Night': nights}
        for i in range(len(split_QCnames)):
            data_dict_check5[split_QCnames[i]] = split_QCvalues[:,i]

        inspect_check5 = atf.inspect_plot('localisation_test1',
                    'check5',
                    data_dict_check5,
                    'cal_localisation_{0}.py Quality Control'.format(
                    instrument.lower())
                    )

        # checking for duplicates
        if (len(output1_calibDB) < output1_num_entry
                or len(output2_calibDB) < output2_num_entry):

            (file_col_output1_unique,
             index_output1_dup,
             return_counts_output1) = np.unique(
                         file_col[index_key_output1][~master_col[index_key_output1]],
                         return_index=True,
                         return_counts=True
                         )
            (file_col_output2_unique,
             index_output2_dup,
             return_counts_output2) = np.unique(
                         file_col[index_key_output2][~master_col[index_key_output2]],
                         return_index=True,
                         return_counts=True
                         )

        # stop2
        if (output1_num_calibDB == output1_num_entry
                and output2_num_calibDB == output2_num_entry):
            color_stop2 = 'Lime'
            result_stop2 = 'Yes'
            comment_stop2 = ''
            inspect_stop2 = ''

        elif (output1_num_calibDB < output1_num_entry
                or output2_num_calibDB < output2_num_entry):

            color_stop2 = 'Yellow'
            result_stop2 = 'No'

            # if missing output
            if len(output1_missing_calibDB) > 0 or len(output2_missing_calibDB) > 0:

                comment_stop2 = 'One or more outputs are not in the calibDB.'
                data_dict_stop2 = {'Night': np.concatenate(
                                                    (night_output1_missing_calibDB,
                                                     night_output2_missing_calibDB)
                                            ),
                                   'File name': np.concatenate(
                                                        (output1_missing_calibDB,
                                                         output2_missing_calibDB)
                                                ),
                                   }
                inspect_stop2 = atf.inspect_table('badpixel_test1',
                                                  'stop2',
                                                  data_dict_stop2,
                                                  'Missing Output in {0}'.format(
                                                      calibDB_path)
                                                  )
            # if duplicates
            else:

                comment_stop2 = ('Some entries in master_calib_{0}.txt are '
                                 'identical.').format(instrument)
                data_dict_stop2 = {'Night': np.concatenate(
                                                    (night_output1_dup_calibDB,
                                                     night_output2_dup_calibDB)
                                            ),
                                   'File name': np.concatenate(
                                                        (output1_dup_calibDB,
                                                         output2_dup_calibDB)
                                                        ),
                                   'Occurrence': np.concatenate(
                                                (return_counts_output1[count_mask1],
                                                 return_counts_output2[count_mask2])
                                                 )
                                   }
                inspect_stop2 = atf.inspect_table(
                        'badpixel_test1',
                        'stop2',
                        data_dict_stop2,
                        ('Duplicate Entries in '
                         'master_calib_{0}.txt').format(instrument)
                        )

        else:
            color_stop2 = 'Red'
            result_stop2 = 'No'
            comment_stop2 = ('The calibDB should not have more output files than what '
                             'was produced.')
            inspect_stop2 = ''

        # TODO: use jinja template
        html_text = ''
        with open('localisation_test1/localisation_test1.html', 'w') as f:
            f.write(html_text)


if __name__ == '__main__':
    test = LocTest()
    test.runtest()
