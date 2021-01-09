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
import numpy as np

from checks.apero_tests import Test
import apero_tests_func as atf
from apero.core import constants


class PPTest(Test):
    """PPTEST."""

    def __init__(self):
        """__init__."""
        
        # TODO: Call parent init

        self._name = 'Preprocessing Recipe Test #1'

        raw_path = params['DRS_DATA_RAW']
        if raw_path[-1] == '/':
            raw_path = raw_path[:-1]            # raw path without / at the end
        raw_nights = atf.list_nights(raw_path)  # list raw night directories

        pp_path = params['DRS_DATA_WORKING']
        if pp_path[-1] == '/':
            pp_path = pp_path[:-1]  # pp path without / at the end
        pp_nights = atf.list_nights(pp_path)  # list PP data night directories

        # output list
        output_list = ['*_pp.fits']


    @property
    def name(self):
        """name."""
        return self._name

    def runtest(self):
        """runtest."""

        print(self.name)

        # =============================================================================
        # TESTS
        # =============================================================================
        # Check 1
        raw_num = atf.count_files_subdir(raw_path,
                                         subdir='all',
                                         files='*.fits')

        # Check 2
        pp_num = atf.count_files_subdir(pp_path,
                                        subdir='all',
                                        files=output_list[0])

        # Stop 1
        if pp_num == raw_num:
            color_stop1 = 'Lime'
            result_stop1 = 'Yes'
            comment_stop1 = ''
            inspect_stop1 = ''
        elif pp_num < raw_num:
            color_stop1 = 'Yellow'
            result_stop1 = 'No'
            comment_stop1 = 'Not all available raw files were reduced.'
            inspect_stop1 = ''
        else:
            color_stop1 = 'Red'
            result_stop1 = 'No'
            comment_stop1 = ('The number of pp files should always be '
                             'smaller than the number of raw files.')
            inspect_stop1 = ''


        # Inspect all pp_nights index.fits and log.fits
        pp_num_indexfits = 0                  # Check 3
        pp_num_logfits = 0                    # Check 4
        pp_num_logfits_unique = 0             # Check 5
        pp_num_logfits_unique_QCfalse = 0     # Check 6
        pp_num_logfits_unique_ENDEDfalse = 0  # Check 7

        odometers_logfits_QCfalse = []        # Check 6
        nights_logfits_QCfalse = []           # Check 6
        QCstr_logfits_QCfalse = []            # Check 6

        odometers_logfits_ENDEDfalse = []     # Check 7
        nights_logfits_ENDEDfalse = []        # Check 7
        ERRORS_logfits_ENDEDfalse = []        # Check 7
        LOGFILE_logfits_ENDEDfalse = []       # Check 7

        missing_indexfits = []
        missing_logfits = []

        for i in range(len(pp_nights)):

            # Inspect index.fits if the file exists
            if os.path.isfile('{0}/{1}/index.fits'.format(pp_path, pp_nights[i])):

                indexfits = atf.index_fits('{0}/{1}/index.fits'.format(pp_path,
                                                                       pp_nights[i])
                                           )
                pp_num_indexfits += indexfits.len  # Check 3

            # Missing index.fits
            else:
                missing_indexfits.append('{0}/{1}/index.fits'.format(pp_path,
                                                                     pp_nights[i])
                                         )

            # Inspect log.fits if the file exists
            if os.path.isfile('{0}/{1}/log.fits'.format(pp_path, pp_nights[i])):

                logfits = atf.log_fits('{0}/{1}/log.fits'.format(pp_path,
                                                                 pp_nights[i])
                                       )
                pp_num_logfits += logfits.len  # Check 4

                # Don't consider duplicates pp files
                args = logfits.args

                odometers = []

                for j in range(len(args)):
                    if 'persi_' in args[j]:
                        args[j].replace('persi_', '')
                    index = args[j].index('.fits')
                    odometers.append(args[j][index-8:index])

                odometers, index_unique = np.unique(odometers, return_index=True)

                tbl_unique = logfits.tbl[index_unique]
                pp_num_logfits_unique += len(tbl_unique)  # Check 5

                indexQCfalse = ~tbl_unique['PASSED_ALL_QC']
                indexENDEDfalse = ~tbl_unique['ENDED']

                # Check 6
                pp_num_logfits_unique_QCfalse += sum(indexQCfalse)
                odometers_logfits_QCfalse.extend(odometers[indexQCfalse])
                nights_logfits_QCfalse.extend(tbl_unique['DIRECTORY'][indexQCfalse])
                QCstr_logfits_QCfalse.extend(tbl_unique['QC_STRING'][indexQCfalse])

                # Check 7
                pp_num_logfits_unique_ENDEDfalse += sum(indexENDEDfalse)
                odometers_logfits_ENDEDfalse.extend(odometers[indexENDEDfalse])
                nights_logfits_ENDEDfalse.extend(
                        tbl_unique['DIRECTORY'][indexENDEDfalse]
                        )
                ERRORS_logfits_ENDEDfalse.extend(tbl_unique['ERRORS'][indexENDEDfalse])
                LOGFILE_logfits_ENDEDfalse.extend(
                        tbl_unique['LOGFILE'][indexENDEDfalse]
                        )

            # Missing log.fits
            else:
                missing_logfits.append('{0}/{1}/log.fits'.format(pp_path, pp_nights[i]))


        # Stop2
        if pp_num_indexfits == pp_num:
            color_stop2 = 'Lime'
            result_stop2 = 'Yes'
            comment_stop2 = ''
            inspect_stop2 = ''
        else:
            color_stop2 = 'Red'
            result_stop2 = 'No'
            comment_stop2 = ''
            inspect_stop2 = ''

        # Stop 3
        if pp_num_logfits_unique == pp_num:
            color_stop3 = 'Lime'
            result_stop3 = 'Yes'
            comment_stop3 = ''
            inspect_stop3 = ''
        elif pp_num_logfits_unique > pp_num:
            color_stop3 = 'Yellow'
            result_stop3 = 'No'
            comment_stop3 = 'Some files were processed more than once.'
            inspect_stop3 = ''
        else:
            color_stop3 = 'Red'
            result_stop3 = 'No'
            comment_stop3 = ''
            inspect_stop3 = ''


        data_dict_check6 = {'Night': nights_logfits_QCfalse,
                            'Odometer': odometers_logfits_QCfalse,
                            'QC_STRING': QCstr_logfits_QCfalse,
                            }
        inspect_check6 = atf.inspect_table('preprocessing_test1',
                                           'check6',
                                           data_dict_check6,
                                           ('Odometers that Failed One or More'
                                            'Quality Control')
                                           )

        data_dict_check7 = {'Night': nights_logfits_ENDEDfalse,
                            'Odometer': odometers_logfits_ENDEDfalse,
                            'ERRORS': ERRORS_logfits_ENDEDfalse,
                            'LOGFILE': LOGFILE_logfits_ENDEDfalse,
                            }
        inspect_check7 = atf.inspect_table('preprocessing_test1',
                                           'check7',
                                           data_dict_check7,
                                           'Odometers that Failed to Finish'
                                           )


        # TODO: use jinja template
        html_text = ''
        with open('preprocessing_test1/preprocessing_test1.html', 'w') as f:
            f.write(html_text)


if __name__ == '__main__':
    test = PPTest()
    test.runtest()
