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
        return 'cal_dark_master_{}'.format(self.instrument.lower())

    def runtest(self):
        """runtest."""

        # TODO: overwrite counting to include master (is a master recipe)
        # Either ismaster property in calib parent or just reimplement

        # TODO: properly re-align log (other dir) with output (nights)

        # checking for missing output
        if (recipe_num_logfits > output1_num_unique and len(np.unique(logfits.nights[index_recipe])) > len(np.unique(nights_output1))):
            night_output1_missing = logfits.nights[index_recipe][np.in1d(logfits.nights[index_recipe], ~np.unique(nights_output1))]
            output1_missing = output_list * len(night_output1_missing)

        # checking for duplicates
        if recipe_num_logfits > output1_num_unique and len(np.unique(logfits.nights[index_recipe])) == len(np.unique(nights_output1)):

            nights_logfits_unique, return_counts_output1 = np.unique(logfits.nights[index_recipe], return_counts=True)

            night_output1_dup = nights_logfits_unique[return_counts_output1 > 1]
            output1_dup = output_list * len(night_output1_dup)

        indexQCfalse = np.array(logfits.indexQCfalse) * np.array(index_recipe)  # QC false + correct recipe
        indexENDEDfalse = np.array(logfits.indexENDEDfalse) * np.array(index_recipe)  # ENDED false + correct recipe

        num_logfits_QCfalse = sum(indexQCfalse)  # check4
        num_logfits_ENDEDfalse = sum(indexENDEDfalse)  # check5

        nights_logfits_QCfalse = logfits.nights[indexQCfalse]  # check4
        QCstr_logfits_QCfalse = logfits.QCstr[indexQCfalse]    # check4

        nights_logfits_ENDEDfalse = logfits.nights[indexENDEDfalse]    # check5
        ERRORS_logfits_ENDEDfalse = logfits.ERRORS[indexENDEDfalse]    # check5
        LOGFILE_logfits_ENDEDfalse = logfits.LOGFILE[indexENDEDfalse]  # check5


        # check4
        if num_logfits_QCfalse == 0:
            comments_check4 = ''
            inspect_check4 = ''
        else:
            comments_check4 = 'One or more recipe have failed QC.'
            data_dict_check4 = {'Night': nights_logfits_QCfalse,
                                'QC_STRING': QCstr_logfits_QCfalse,
                                }
            inspect_check4 = atf.inspect_table('darkmaster_test1',
                                               'check4',
                                               data_dict_check4,
                                               'Nights that Failed Quality Control'
                                               )

        # check5
        if num_logfits_ENDEDfalse == 0:
            comments_check5 = ''
            inspect_check5 = ''
        else:
            comments_check5 = 'One or more recipe have failed to finish.'
            data_dict_check5 = {'Night': nights_logfits_ENDEDfalse,
                                'ERRORS': ERRORS_logfits_ENDEDfalse,
                                'LOGFILE': LOGFILE_logfits_ENDEDfalse,
                                }
            inspect_check5 = atf.inspect_table('darkmaster_test1',
                                               'check5',
                                               data_dict_check5,
                                               'Nights that Failed to Finish')


        # stop1
        if output1_num_unique == recipe_num_logfits:
            color_stop1 = 'Lime'
            result_stop1 = 'Yes'
            comment_stop1 = ''
            inspect_stop1 = ''
        elif output1_num_unique < recipe_num_logfits:

            color_stop1 = 'Yellow'
            result_stop1 = 'No'

            # if missing output
            if len(output1_missing) > 0:
                comment_stop1 = 'One or more output were not produced.'
                data_dict_stop1 = {'Night': night_output1_missing,
                                   'File name': output1_missing,
                                   }
                inspect_stop1 = atf.inspect_table('darkmaster_test1',
                                                  'stop1',
                                                  data_dict_stop1,
                                                  'Missing Outputs in {0}'.format(
                                                      reduced_path
                                                      ),
                                                  )
            # if duplicates
            else:
                comment_stop1 = ('Recipe called multiple times producing the same '
                                 'outputs.')
                data_dict_stop1 = {'Night': night_output1_dup,
                                   'File name': output1_dup,
                                   }
                inspect_stop1 = atf.inspect_table(
                        'darkmaster_test1',
                        'stop1',
                        data_dict_stop1,
                        ('Same Dark Master Recipe Called Twice or More '
                         'Producing the Same Outputs in {0}'
                         ).format(reduced_path)
                         )
        else:
            color_stop1 = 'Red'
            result_stop1 = 'No'
            comment_stop1 = ('The number of unique output files should always be '
                             'smaller or equal to the number of recipe called.'
                             )
            inspect_stop1 = ''

        # Inspect calibDB
        f = open("{0}/master_calib_{1}.txt".format(calibDB_path, instrument), "r")
        master_calib_txt = f.read()
        nprocessed = master_calib_txt.index('# DRS processed')
        index_start = master_calib_txt[:nprocessed].count('\n')
        key_col, master_col, night_col, file_col = np.genfromtxt(
                "{0}/master_calib_{1}.txt".format(calibDB_path, instrument),
                delimiter=' ',
                unpack=True,
                usecols=(0, 1, 2, 3),
                skip_header=index_start,
                dtype=str)
        master_col = master_col.astype(dtype = bool) # str to bool

        index_key_output1 = key_col == calibDB_entry_list[0]

        output1_num_entry = len(key_col[index_key_output1]) # check6

        output1_calibDB = atf.list_files("{0}".format(calibDB_path),
                                         files=output_list[0][1:])
        output1_num_calibDB = len(output1_calibDB)  # check7


        # checking for missing output

        output1_missing_calibDB = []

        if sum(np.in1d(file_col[index_key_output1], output1_calibDB)) < len(file_col[index_key_output1]): 

            index_output1_missing = ~np.in1d(
                    file_col[index_key_output1],
                    output1_calibDB)

            night_output1_missing_calibDB = night_col[index_key_output1][index_output1_missing]
            output1_missing_calibDB = file_col[index_key_output1][index_output1_missing]


        # checking for duplicates

        if len(output1_calibDB) < output1_num_entry:

            (file_col_output1_unique,
             index_output1_dup,
             return_counts_output1) = np.unique(
                         file_col[index_key_output1],
                         return_index=True,
                         return_counts=True
                         )

            count_mask1 = return_counts_output1 > 1

            night_output1_dup_calibDB = night_col[index_key_output1][index_output1_dup][count_mask1]
            output1_dup_calibDB = file_col_output1_unique[count_mask1]


        # stop2
        if output1_num_calibDB == output1_num_entry:
            color_stop2 = 'Lime'
            result_stop2 = 'Yes'
            comment_stop2 = ''
            inspect_stop2 = ''

        elif output1_num_calibDB < output1_num_entry:

            color_stop2 = 'Yellow'
            result_stop2 = 'No'

            # if missing output
            if len(output1_missing_calibDB) > 0:

                comment_stop2 = 'One or more outputs are not in the calibDB.'
                data_dict_stop2 = {'Night': night_output1_missing_calibDB,
                                   'File name': output1_missing_calibDB,
                                   }
                inspect_stop2 = atf.inspect_table('darkmaster_test1',
                                                  'stop2',
                                                  data_dict_stop2,
                                                  'Missing Output in {0}'.format(
                                                      calibDB_path)
                                                  )
            # if duplicates
            else:

                comment_stop2 = ('Some entries in master_calib_{0}.txt are '
                                 'identical.').format(instrument)
                data_dict_stop2 = {'Night': night_output1_dup_calibDB,
                                   'File name': output1_dup_calibDB,
                                   'Occurrence': return_counts_output1[count_mask1]
                                   }
                inspect_stop2 = atf.inspect_table(
                        'darkmaster',
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


        # TODO: fill template
        html_text = ''
        with open('darkmaster_test1/darkmaster_test1.html', 'w') as f:
            f.write(html_text)


if __name__ == '__main__':
    test = DarkMTest()
    test.runtest()