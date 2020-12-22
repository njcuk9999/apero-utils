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
import os
import glob
from datetime import datetime
import numpy as np
from astropy.io import fits

from apero_tests import Test
import apero_tests_func as atf
from apero.core import constants


class LocTest(Test):
    """LocTest."""

    def __init__(self):
        """__init__."""

        self._name = 'Localisation Recipe Test #1'

    @property
    def name(self):
        """name."""
        return self._name

    def runtest(self):
        """runtest."""

        print(self.name)
        # =============================================================================
        # Define Constants
        # =============================================================================
        params = constants.load('SPIROU')

        setup = os.environ['DRS_UCONFIG']  # setup
        instrument = params['INSTRUMENT']  # instrument
        date = datetime.now()
        date = date.strftime("%Y-%m-%d %H:%M:%S")  # date

        reduced_path = params['DRS_DATA_REDUC']
        if reduced_path[-1] == '/':
            reduced_path = reduced_path[:-1]        # reduced path without / at the end
        reduced_nights = atf.list_nights(reduced_path)  # list reduced night dirs

        calibDB_path = params['DRS_CALIB_DB']  # calibDB path
        if calibDB_path[-1] == '/':
            calibDB_path = calibDB_path[:-1]   # calibDB path without / at the end

        # output list
        output_list = ['*_pp_order_profile_{FIBER}.fits',
                       '*_pp_loco_{FIBER}.fits',
                       '*_pp_fwhm-order_{FIBER}.fits',
                       '*_pp_with-order_{FIBER}.fits']


        # calibDB entries list
        calibDB_entry_list = ['ORDER_PROFILE_{FIBER}', 'LOC_{FIBER}']


        # =============================================================================
        # TESTS
        # =============================================================================
        # inspect all reduced_nights log.fits

        tot_recipe_num_logfits = 0       # check1
        master_recipe_num_logfits = 0    # check1
        num_logfits_QCfalse = 0          # check4
        num_logfits_ENDEDfalse = 0       # check6

        nights_logfits_QCfalse = []      # check4
        QCstr_logfits_QCfalse = []       # check4

        nights = []                      # check5
        QCnames = []                     # check5
        QCvalues = []                    # check5

        nights_logfits_ENDEDfalse = []   # check6
        ERRORS_logfits_ENDEDfalse = []   # check6
        LOGFILE_logfits_ENDEDfalse = []  # check6

        calibration_night_missing = []   # check9
        calibration_missing = []         # check9
        calibration_night = []           # check9

        output1 = []
        output2 = []
        output3 = []
        output4 = []

        night_output1_missing = []
        output1_missing = []
        night_output2_missing = []
        output2_missing = []
        night_output3_missing = []
        output3_missing = []
        night_output4_missing = []
        output4_missing = []

        night_output1_dup = []
        output1_dup = []
        night_output2_dup = []
        output2_dup = []
        night_output3_dup = []
        output3_dup = []
        night_output4_dup = []
        output4_dup = []

        missing_logfits = []

        for i in range(len(reduced_nights)):

            output1_list_files = atf.list_files(
                    '{0}/{1}'.format(reduced_path, reduced_nights[i]),
                    files=(output_list[0][1:].replace('{FIBER}', 'AB'),
                           output_list[0][1:].replace('{FIBER}', 'C'))
                    )
            output2_list_files = atf.list_files(
                    '{0}/{1}'.format(reduced_path, reduced_nights[i]),
                    files=(output_list[1][1:].replace('{FIBER}', 'AB'),
                           output_list[1][1:].replace('{FIBER}', 'C'))
                    )
            output3_list_files = atf.list_files(
                    '{0}/{1}'.format(reduced_path, reduced_nights[i]),
                    files=(output_list[2][1:].replace('{FIBER}', 'AB'),
                           output_list[2][1:].replace('{FIBER}', 'C'))
                    )
            output4_list_files = atf.list_files(
                    '{0}/{1}'.format(reduced_path, reduced_nights[i]),
                    files=(output_list[3][1:].replace('{FIBER}', 'AB'),
                           output_list[3][1:].replace('{FIBER}', 'C'))
                    )

            output1.extend(output1_list_files)
            output2.extend(output2_list_files)
            output3.extend(output3_list_files)
            output4.extend(output4_list_files)

            # inspect log.fits if the file exists
            if os.path.isfile('{0}/{1}/log.fits'.format(reduced_path,
                                                        reduced_nights[i])
                              ):
                logfits = atf.log_fits('{0}/{1}/log.fits'.format(reduced_path,
                                       reduced_nights[i])
                                       )

                index_recipe = (
                        logfits.recipe == 'cal_loc_{0}'.format(instrument.lower())
                        )
                tot_recipe_num_logfits += sum(index_recipe)  # check1

                # checking for missing output
                if sum(index_recipe) > 0 and len(output1_list_files) == 0:
                    night_output1_missing.append(reduced_nights[i])
                    output1_missing.append(output_list[0])
                if sum(index_recipe) > 0 and len(output2_list_files) == 0:
                    night_output2_missing.append(reduced_nights[i])
                    output2_missing.append(output_list[1])
                if sum(index_recipe) > 0 and len(output3_list_files) == 0:
                    night_output3_missing.append(reduced_nights[i])
                    output3_missing.append(output_list[2])
                if sum(index_recipe) > 0 and len(output4_list_files) == 0:
                    night_output4_missing.append(reduced_nights[i])
                    output4_missing.append(output_list[3])

                # checking for duplicates
                if sum(index_recipe) > len(output1_list_files):

                    # checking if it is in the master directory
                    for j in range(sum(index_recipe)):
                        if '--master=True' in logfits.runstr[index_recipe][j]:
                            master_recipe_num_logfits += 1
                            comments_check1 = (
                                    'An additional {0} recipes with '
                                    '--master=True were called in the master '
                                    'directory {1}.'.format(master_recipe_num_logfits,
                                                            reduced_nights[i])
                                    )

                    non_dup_num = sum(index_recipe) - master_recipe_num_logfits

                    if non_dup_num == len(output1_list_files):
                        pass
                    else:
                        night_output1_dup.append(reduced_nights[i])
                        output1_dup.append(output_list[0])
                    if non_dup_num == len(output2_list_files):
                        pass
                    else:
                        night_output2_dup.append(reduced_nights[i])
                        output2_dup.append(output_list[1])
                    if non_dup_num == len(output3_list_files):
                        pass
                    else:
                        night_output3_dup.append(reduced_nights[i])
                        output3_dup.append(output_list[2])
                    if non_dup_num == len(output4_list_files):
                        pass
                    else:
                        night_output4_dup.append(reduced_nights[i])
                        output4_dup.append(output_list[3])


                indexQCfalse = np.array(logfits.indexQCfalse) * np.array(index_recipe)  # QC false + correct recipe
                indexENDEDfalse = np.array(logfits.indexENDEDfalse) * np.array(index_recipe)  # ENDED false + correct recipe

                num_logfits_QCfalse += sum(indexQCfalse)  # check4        
                num_logfits_ENDEDfalse += sum(indexENDEDfalse)  # check6

                # check4
                nights_logfits_QCfalse.extend(logfits.nights[indexQCfalse])
                QCstr_logfits_QCfalse.extend(logfits.QCstr[indexQCfalse])

                # check5    
                nights.extend(logfits.nights[index_recipe]) 
                QCnames.extend(logfits.QCnames[index_recipe]) 
                QCvalues.extend(logfits.QCvalues[index_recipe])

                # check6
                nights_logfits_ENDEDfalse.extend(logfits.nights[indexENDEDfalse])
                ERRORS_logfits_ENDEDfalse.extend(logfits.ERRORS[indexENDEDfalse])
                LOGFILE_logfits_ENDEDfalse.extend(logfits.LOGFILE[indexENDEDfalse])

                # check9
                PID = logfits.PID[index_recipe]
                for k in range(len(output1_list_files)):
                    hdul = fits.open('{0}/{1}/{2}'.format(reduced_path,
                            reduced_nights[i], output1_list_files[k]))
                    if hdul[0].header['DRSPID'] in PID:
                        if 'CDBBAD' in hdul[0].header:
                            if not (hdul[0].header['CDBBAD'] == 'None' or 
                                    os.path.isfile('{0}/{1}/{2}'.format(reduced_path,
                                    reduced_nights[i], hdul[0].header['CDBBAD']))):
                                calibration_night_missing.append(reduced_nights[i])
                                calibration_missing.append(hdul[0].header['CDBBAD'])
                        if 'CDBBACK' in hdul[0].header:      
                            if not (hdul[0].header['CDBBACK'] == 'None' or 
                                    os.path.isfile('{0}/{1}/{2}'.format(reduced_path,
                                    reduced_nights[i], hdul[0].header['CDBBACK']))):
                                calibration_night_missing.append(reduced_nights[i])
                                calibration_missing.append(hdul[0].header['CDBBACK'])
                    else:
                        continue


            # missing log.fits
            else:
                missing_logfits.append('{0}/{1}/log.fits'.format(reduced_path,
                                                                 reduced_nights[i])
                                       )
        recipe_num_logfits = tot_recipe_num_logfits-master_recipe_num_logfits  # check1
        output1_num = len(output1)                    # check2
        output1_num_unique = len(np.unique(output1))  # check3
        output2_num = len(output2)                    # check2
        output2_num_unique = len(np.unique(output2))  # check3
        output3_num = len(output3)                    # check2
        output3_num_unique = len(np.unique(output3))  # check3
        output4_num = len(output3)                    # check2
        output4_num_unique = len(np.unique(output4))  # check3


        # check4
        if num_logfits_QCfalse == 0:
            comments_check4 = ''
            inspect_check4 = ''
        else:
            comments_check4 = 'One or more recipe have failed QC.'
            data_dict_check4 = {'Night': nights_logfits_QCfalse,
                                'QC_STRING': QCstr_logfits_QCfalse,
                                }
            inspect_check4 = atf.inspect_table(
                    'localisation_test1',
                    'check4',
                    data_dict_check4,
                    'Nights that Failed Quality Control'
                    )

        #check5
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

        # check6
        if num_logfits_ENDEDfalse == 0:
            comments_check6 = ''
            inspect_check6 = ''
        else:
            comments_check6 = 'One or more recipe have failed to finish.'
            data_dict_check6 = {'Night': nights_logfits_ENDEDfalse,
                                'ERRORS': ERRORS_logfits_ENDEDfalse,
                                'LOGFILE': LOGFILE_logfits_ENDEDfalse,
                                }
            inspect_check6 = atf.inspect_table(
                    'localisation_test1',
                    'check6',
                    data_dict_check6,
                    'Nights that Failed to Finish'
                    )


        # stop1
        if (output1_num_unique == recipe_num_logfits
                and output2_num_unique == recipe_num_logfits
                and output3_num_unique == recipe_num_logfits
                and output4_num_unique == recipe_num_logfits):
            color_stop1 = 'Lime'
            result_stop1 = 'Yes'
            comment_stop1 = ''
            inspect_stop1 = ''

        elif (output1_num_unique < recipe_num_logfits
                or output2_num_unique < recipe_num_logfits
                or output3_num_unique < recipe_num_logfits
                or output4_num_unique < recipe_num_logfits):

            color_stop1 = 'Yellow'
            result_stop1 = 'No'

            # if missing outputs
            if (len(output1_missing) > 0
                    or len(output2_missing) > 0
                    or len(output3_missing) > 0
                    or len(output4_missing) > 0):

                comment_stop1 = 'One or more outputs were not produced.'
                data_dict_stop1 = {'Night': np.concatenate((night_output1_missing,
                                                            night_output2_missing,
                                                            night_output3_missing,
                                                            night_output4_missing)),
                                   'File name': np.concatenate((output1_missing,
                                                                output2_missing,
                                                                output3_missing,
                                                                output4_missing)),
                                   }
                inspect_stop1 = atf.inspect_table(
                        'localisation_test1',
                        'stop1',
                        data_dict_stop1,
                        'Missing Outputs in {0}'.format(reduced_path)
                        )

            # if duplicates
            else:
                comment_stop1 = ('Recipe called multiple times producing the same '
                                 'outputs.')
                data_dict_stop1 = {'Night': np.concatenate((night_output1_dup,
                                                            night_output2_dup,
                                                            night_output3_dup,
                                                            night_output4_dup)),
                                   'File name': np.concatenate((output1_dup,
                                                                output2_dup,
                                                                output3_dup,
                                                                output4_dup)),
                                   }
                inspect_stop1 = atf.inspect_table(
                        'localisation_test1',
                        'stop1',
                        data_dict_stop1,
                        ('Same Localisation Recipe Called Twice or More Producing '
                         'the Same Outputs in {0}'
                         ).format(reduced_path)
                        )

        else:
            color_stop1 = 'Red'
            result_stop1 = 'No'
            comment_stop1 = ('The number of unique output files should always be '
                             'smaller than or equal to the number of recipe called.')
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
        master_col = master_col.astype(dtype=bool)  # str to bool

        index_key_output1 = np.logical_or(
                key_col == calibDB_entry_list[0].replace('{FIBER}', 'AB'),
                key_col == calibDB_entry_list[0].replace('{FIBER}', 'C')
                )
        index_key_output2 = np.logical_or(
                key_col == calibDB_entry_list[1].replace('{FIBER}', 'AB'),
                key_col == calibDB_entry_list[1].replace('{FIBER}', 'C')
                )

        tot_output1_num_entry = len(key_col[index_key_output1])  # check7
        master_output1_num_entry = np.sum(master_col[index_key_output1])
        output1_num_entry = tot_output1_num_entry - master_output1_num_entry
        output1_calibDB = atf.list_files(
                "{0}".format(calibDB_path),
                files=(output_list[0][1:].replace('{FIBER}', 'AB'),
                       output_list[0][1:].replace('{FIBER}', 'C'))
                )
        output1_num_calibDB = len(output1_calibDB)  # check8

        tot_output2_num_entry = len(key_col[index_key_output2])  # check7
        master_output2_num_entry = np.sum(master_col[index_key_output2])
        output2_num_entry = tot_output2_num_entry - master_output2_num_entry
        output2_calibDB = atf.list_files(
                "{0}".format(calibDB_path),
                files=(output_list[1][1:].replace('{FIBER}', 'AB'),
                       output_list[1][1:].replace('{FIBER}', 'C'))
                )
        output2_num_calibDB = len(output2_calibDB)  # check8

        comments_check7 = ('An additional {0} {1} and {2} {3} with master = 1 '
                           'are in the master_calibDB_{4}.'.format(
                               master_output1_num_entry,
                               calibDB_entry_list[0], master_output2_num_entry,
                               calibDB_entry_list[1], instrument)
                           )


        # checking for missing output
        output1_missing_calibDB = []
        output2_missing_calibDB = []
        if (sum(np.in1d(file_col[index_key_output1], output1_calibDB)) < len(file_col[index_key_output1])
                or sum(np.in1d(file_col[index_key_output2], output2_calibDB)) < len(file_col[index_key_output2])):

            index_output1_missing = ~np.in1d(
                    file_col[index_key_output1],
                    output1_calibDB)
            index_output2_missing = ~np.in1d(
                    file_col[index_key_output2],
                    output2_calibDB)

            night_output1_missing_calibDB = night_col[index_key_output1][index_output1_missing]
            output1_missing_calibDB = file_col[index_key_output1][index_output1_missing]
            night_output2_missing_calibDB = night_col[index_key_output2][index_output2_missing]
            output2_missing_calibDB = file_col[index_key_output2][index_output2_missing]


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

            # Keeping long lines for now, will split later
            # (indices on next line were hard to read)
            count_mask1 = return_counts_output1 > 1
            count_mask2 = return_counts_output2 > 1
            night_output1_dup_calibDB = night_col[index_key_output1][~master_col[index_key_output1]][index_output1_dup][count_mask1]
            output1_dup_calibDB = file_col_output1_unique[count_mask1]
            night_output2_dup_calibDB = night_col[index_key_output2][~master_col[index_key_output2]][index_output2_dup][count_mask2]
            output2_dup_calibDB = file_col_output2_unique[count_mask2]


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
                inspect_stop2 = atf.inspect_table('localisation_test1',
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
                        'localisation_test1',
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


        # check9
        for i in range(len(calibration_missing)):
            path = glob.glob('{0}/*/{1}'.format(reduced_path,
                    calibration_missing[i]))
            path = path[0].replace('{0}/'.format(reduced_path), '')
            calibration_night.append(path[:path.index('/')])
           

        if len(calibration_missing) == 0:
            comments_check9 = ''
            inspect_check9 = ''
        else:
            comments_check9 = ('One or more localisation recipe outputs used '
                               'the bad pixel calibrations from another night')
            data_dict_check9 = {'Localisation Night':calibration_night_missing,
                                'Calibration file name': calibration_missing,
                                'Calibration Night': calibration_night
                                }
            inspect_check9 = atf.inspect_table(
                        'localisation_test1',
                        'check9',
                        data_dict_check9,
                        ('Night where the Localisation Recipe Output '
                         'used Calibrations from Another Night')
                        )    

        # Build localisation_test1.html

        html_text = f"""
        <html>


        <head>
        <title>APERO Tests</title>
        <style>
        table {{
          width:75%;
        }}
        table, th, td {{
          border: 1px solid black;
          border-collapse: collapse;
        }}
        th, td {{
          padding: 15px;
          text-align: left;
        }}
        #t01 tr:nth-child(even) {{
          background-color: #eee;
        }}
        #t01 tr:nth-child(odd) {{
         background-color: #fff;
        }}
        #t01 th {{
          background-color: white;
          color: black;
        }}
        </style>
        </head>

        <body>

        <h3>Localisation Recipe Test #1</h3>
        <p><b>Setup: {setup}</b><br>
        <p><b>Instrument: {instrument}</b><br>
        <p><b>Date: {date}</b><br>
        <br>
        <p>Script: cal_loc_{instrument.lower()}.py<br>
        <p>Output files: {', '.join(output_list)}<br>
        <p>Calibration database entry: {', '.join(calibDB_entry_list)}<br>
        <p><a href='https://github.com/njcuk9999/apero-drs#84-localisation-recipe'>Link</a> to Localisation Recipe description</p>
        <br></br>

        <table id="t01">

          <colgroup>
             <col span="1" style="width: 5%;">
             <col span="1" style="width: 45%;">
             <col span="1" style="width: 21%;">
             <col span="1" style="width: 24%;">
             <col span="1" style="width: 5%;">
          </colgroup>
          <tr>
            <th>Check</th>
            <th>Description</th>
            <th>Result</th>
            <th>Comments</th>
            <th>Details</th>
          </tr>
          <tr>
            <td>1</td>
            <td># of time cal_loc_{instrument.lower()}.py was called</td>
            <td>{recipe_num_logfits}</td>
            <td>{comments_check1}</td>
            <td></td>
          </tr>
          <tr>
            <td>2</td>
            <td># of outputs in {reduced_path}</td>
            <td>{output_list[0]}: {output1_num}<br>{output_list[1]}: {output2_num}<br>{output_list[2]}: {output3_num}<br>{output_list[3]}: {output4_num}</td>
            <td></td>
            <td></td>
          </tr>
          <tr>
            <td>3</td>
            <td># of unique outputs in {reduced_path}</td>
            <td>{output_list[0]}: {output1_num_unique}<br>{output_list[1]}: {output2_num_unique}<br>{output_list[2]}: {output3_num_unique}<br>{output_list[3]}: {output4_num_unique}</td>
            <td></td>
            <td></td>
          </tr>
          <tr>
            <td></td>
            <td>Check 3 == Check 1?</td>
            <td bgcolor={color_stop1}>{result_stop1}</td>
            <td>{comment_stop1}</td>
            <td>{inspect_stop1}</td>
          </tr>
          <tr>
            <td>4</td>
            <td># of entry in {reduced_path}/log.fits that failed one or more QC</td>
            <td>{num_logfits_QCfalse}</td>
            <td>{comments_check4}</td>
            <td></td>
          </tr>
          <tr>
            <td>5</td>
            <td>Plot the different QCs as a function of time</td>
            <td>{inspect_check5}</td>
            <td></td>
            <td></td>
          </tr>
          <tr>
            <td>6</td>
            <td># of entry in {reduced_path}/log.fits that failed to finish</td>
            <td>{num_logfits_ENDEDfalse}</td>
            <td>{comments_check6}</td>
            <td></td>
          </tr>
          <tr>
            <td>7</td>
            <td># of {', '.join(calibDB_entry_list)} entry in {calibDB_path}/master_calib_{instrument}.txt</td>
            <td>{calibDB_entry_list[0]}: {output1_num_entry}<br>{calibDB_entry_list[1]}: {output2_num_entry}</td>
            <td>{comments_check7}</td>
            <td></td>
          </tr>
          <tr>
            <td>8</td>
            <td># of outputs in {calibDB_path}</td>
            <td>{output_list[0]}: {output1_num_calibDB}<br>{output_list[1]}: {output2_num_calibDB}</td>
            <td></td>
            <td></td>
          </tr>
          <tr>
            <td></td>
            <td>Check 8 == Check 7?</td>
            <td bgcolor={color_stop2}>{result_stop2}</td>
            <td>{comment_stop2}</td>
            <td>{inspect_stop2}</td>
          </tr>
          <tr>
            <td>9</td>
            <td># of outputs in {reduced_path} that used the bad pixel calibration files from another night</td>
            <td>{len(calibration_missing)}</td>
            <td>{comments_check9}</td>
            <td>{inspect_check9}</td>
          </tr>
        </table>


        </body>
        </html>
        """


        with open('localisation_test1/localisation_test1.html', 'w') as f:
            f.write(html_text)


if __name__ == '__main__':
    test = LocTest()
    test.runtest()
