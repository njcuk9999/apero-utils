"""
Check if localisation calib worked fine.

Tests preformed
check1: how many recipes were run? (cal_loc_{instrument} in log.fits)
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
import os
from datetime import datetime
import numpy as np

import apero_tests_func as atf
from apero.core import constants

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
output_list = ['*_pp_order_profile_AB.fits',
               '*_pp_order_profile_C.fits',
               '*_pp_loco_AB.fits',
               '*_pp_loco_C.fits',
               '*_pp_fwhm-order_AB.fits',
               '*_pp_fwhm-order_C.fits',
               '*_pp_with-order_AB.fits',
               '*_pp_with-order_C.fits']

# calibDB entries list
calibDB_entry_list = ['ORDER_PROFILE_AB', 'ORDER_PROFILE_C', 'LOC_AB', 'LOC_C']


# =============================================================================
# TESTS
# =============================================================================
# inspect all reduced_nights log.fits

recipe_num_logfits = 0      # check1
num_logfits_QCfalse = 0     # check4
num_logfits_ENDEDfalse = 0  # check5

nights_logfits_QCfalse = []  # check4
QCstr_logfits_QCfalse = []   # check4

# check5
nights_logfits_ENDEDfalse = []
ERRORS_logfits_ENDEDfalse = []
LOGFILE_logfits_ENDEDfalse = []

output1 = []
output2 = []
output3 = []
output4 = []
output5 = []
output6 = []
output7 = []
output8 = []

night_output1_missing = []
output1_missing = []
night_output2_missing = []
output2_missing = []
night_output3_missing = []
output3_missing = []
night_output4_missing = []
output4_missing = []
night_output5_missing = []
output5_missing = []
night_output6_missing = []
output6_missing = []
night_output7_missing = []
output7_missing = []
night_output8_missing = []
output8_missing = []

night_output1_dup = []
output1_dup = []
night_output2_dup = []
output2_dup = []
night_output3_dup = []
output3_dup = []
night_output4_dup = []
output4_dup = []
night_output5_dup = []
output5_dup = []
night_output6_dup = []
output6_dup = []
night_output7_dup = []
output7_dup = []
night_output8_dup = []
output8_dup = []

missing_logfits = []

for i in range(len(reduced_nights)):

    output1_list_files = atf.list_files(
            '{0}/{1}'.format(reduced_path, reduced_nights[i]),
            files=output_list[0][1:]
            )
    output2_list_files = atf.list_files(
            '{0}/{1}'.format(reduced_path, reduced_nights[i]),
            files=output_list[1][1:]
            )
    output3_list_files = atf.list_files(
            '{0}/{1}'.format(reduced_path, reduced_nights[i]),
            files=output_list[2][1:]
            )
    output4_list_files = atf.list_files(
            '{0}/{1}'.format(reduced_path, reduced_nights[i]),
            files=output_list[3][1:]
            )
    output5_list_files = atf.list_files(
            '{0}/{1}'.format(reduced_path, reduced_nights[i]),
            files=output_list[4][1:]
            )
    output6_list_files = atf.list_files(
            '{0}/{1}'.format(reduced_path, reduced_nights[i]),
            files=output_list[5][1:]
            )
    output7_list_files = atf.list_files(
            '{0}/{1}'.format(reduced_path, reduced_nights[i]),
            files=output_list[6][1:]
            )
    output8_list_files = atf.list_files(
            '{0}/{1}'.format(reduced_path, reduced_nights[i]),
            files=output_list[7][1:]
            )

    output1.extend(output1_list_files)
    output2.extend(output2_list_files)
    output3.extend(output3_list_files)
    output4.extend(output4_list_files)
    output5.extend(output5_list_files)
    output6.extend(output6_list_files)
    output7.extend(output7_list_files)
    output8.extend(output8_list_files)

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
        recipe_num_logfits += sum(index_recipe)  # check1

        # checking for missing output
        if sum(index_recipe) > 0 and len(output1_list_files) == 0:
            night_output1_missing.extend(reduced_nights[i])
            output1_missing.extend(output_list[0])
        if sum(index_recipe) > 0 and len(output2_list_files) == 0:
            night_output2_missing.extend(reduced_nights[i])
            output2_missing.extend(output_list[1])
        if sum(index_recipe) > 0 and len(output3_list_files) == 0:
            night_output3_missing.extend(reduced_nights[i])
            output3_missing.extend(output_list[2])
        if sum(index_recipe) > 0 and len(output4_list_files) == 0:
            night_output4_missing.extend(reduced_nights[i])
            output4_missing.extend(output_list[3])
        if sum(index_recipe) > 0 and len(output5_list_files) == 0:
            night_output5_missing.extend(reduced_nights[i])
            output5_missing.extend(output_list[4])
        if sum(index_recipe) > 0 and len(output6_list_files) == 0:
            night_output6_missing.extend(reduced_nights[i])
            output6_missing.extend(output_list[5])
        if sum(index_recipe) > 0 and len(output7_list_files) == 0:
            night_output7_missing.extend(reduced_nights[i])
            output7_missing.extend(output_list[6])
        if sum(index_recipe) > 0 and len(output8_list_files) == 0:
            night_output8_missing.extend(reduced_nights[i])
            output8_missing.extend(output_list[7])

        # checking for duplicates
        if sum(index_recipe) > len(output1_list_files):
            night_output1_dup.extend(
                    [reduced_nights[i]] * len(output1_list_files)
                    )
            output1_dup.extend(output1_list_files)
        if sum(index_recipe) > len(output2_list_files):
            night_output2_dup.extend(
                    [reduced_nights[i]] * len(output2_list_files)
                    )
            output2_dup.extend(output2_list_files)
        if sum(index_recipe) > len(output3_list_files):
            night_output3_dup.extend(
                    [reduced_nights[i]] * len(output3_list_files)
                    )
            output3_dup.extend(output3_list_files)
        if sum(index_recipe) > len(output4_list_files):
            night_output4_dup.extend(
                    [reduced_nights[i]] * len(output4_list_files)
                    )
            output4_dup.extend(output4_list_files)
        if sum(index_recipe) > len(output5_list_files):
            night_output5_dup.extend(
                    [reduced_nights[i]] * len(output5_list_files)
                    )
            output5_dup.extend(output5_list_files)
        if sum(index_recipe) > len(output6_list_files):
            night_output6_dup.extend(
                    [reduced_nights[i]] * len(output6_list_files)
                    )
            output6_dup.extend(output6_list_files)
        if sum(index_recipe) > len(output7_list_files):
            night_output7_dup.extend(
                    [reduced_nights[i]] * len(output7_list_files)
                    )
            output7_dup.extend(output7_list_files)
        if sum(index_recipe) > len(output8_list_files):
            night_output8_dup.extend(
                    [reduced_nights[i]] * len(output8_list_files)
                    )
            output8_dup.extend(output8_list_files)


        indexQCfalse = np.array(logfits.indexQCfalse) * np.array(index_recipe)  # QC false + correct recipe
        indexENDEDfalse = np.array(logfits.indexENDEDfalse) * np.array(index_recipe)  # ENDED false + correct recipe

        num_logfits_QCfalse += sum(indexQCfalse)  # check4        
        num_logfits_ENDEDfalse += sum(indexENDEDfalse)  # check5

        # Check 4
        nights_logfits_QCfalse.extend(logfits.nights[indexQCfalse])
        QCstr_logfits_QCfalse.extend(logfits.QCstr[indexQCfalse])

        # Check 5
        nights_logfits_ENDEDfalse.extend(logfits.nights[indexENDEDfalse])
        ERRORS_logfits_ENDEDfalse.extend(logfits.ERRORS[indexENDEDfalse])
        LOGFILE_logfits_ENDEDfalse.extend(logfits.LOGFILE[indexENDEDfalse])

    # missing log.fits
    else:
        missing_logfits.append('{0}/{1}/log.fits'.format(reduced_path,
                                                         reduced_nights[i])
                               )

output1_num = len(output1)                    # check2
output1_num_unique = len(np.unique(output1))  # check3
output2_num = len(output2)                    # check2
output2_num_unique = len(np.unique(output2))  # check3
output3_num = len(output3)                    # check2
output3_num_unique = len(np.unique(output3))  # check3
output4_num = len(output3)                    # check2
output4_num_unique = len(np.unique(output4))  # check3
output5_num = len(output5)                    # check2
output5_num_unique = len(np.unique(output5))  # check3
output6_num = len(output6)                    # check2
output6_num_unique = len(np.unique(output6))  # check3
output7_num = len(output7)                    # check2
output7_num_unique = len(np.unique(output7))  # check3
output8_num = len(output8)                    # check2
output8_num_unique = len(np.unique(output8))  # check3

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
    inspect_check5 = atf.inspect_table(
            'localisation_test1',
            'check5',
            data_dict_check5,
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
        or output4_num_unique < recipe_num_logfits
        or output5_num_unique < recipe_num_logfits
        or output6_num_unique < recipe_num_logfits
        or output7_num_unique < recipe_num_logfits
        or output8_num_unique < recipe_num_logfits):

    color_stop1 = 'Yellow'
    result_stop1 = 'No'

    # if missing output
    if (len(output1_missing) > 0
            or len(output2_missing) > 0
            or len(output3_missing) > 0
            or len(output4_missing) > 0
            or len(output5_missing) > 0
            or len(output6_missing) > 0
            or len(output6_missing) > 0
            or len(output7_missing) > 0
            or len(output8_missing) > 0):
        comment_stop1 = 'One or more output were not produced.'
        data_dict_stop1 = {'Night': np.concatenate((night_output1_missing,
                                                    night_output2_missing,
                                                    night_output3_missing,
                                                    night_output4_missing,
                                                    night_output5_missing,
                                                    night_output6_missing,
                                                    night_output7_missing,
                                                    night_output8_missing)),
                           'File name': np.concatenate((output1_missing,
                                                        output2_missing,
                                                        output3_missing,
                                                        output4_missing,
                                                        output5_missing,
                                                        output6_missing,
                                                        output7_missing,
                                                        output8_missing)),
                           }
        inspect_stop1 = atf.inspect_table(
                'localisation_test1',
                'stop1',
                data_dict_stop1,
                'Missing Outputs in {0}'.format(reduced_path)
                )
    # if duplicates
    else:
        comment_stop1 = ('Recipe called 3 times or more in the same night '
                         'directory.')
        data_dict_stop1 = {'Night': np.concatenate((night_output1_dup,
                                                    night_output2_dup,
                                                    night_output3_dup,
                                                    night_output4_dup,
                                                    night_output5_dup,
                                                    night_output6_dup,
                                                    night_output7_dup,
                                                    night_output8_dup)),
                           'File name': np.concatenate((output1_dup,
                                                        output2_dup,
                                                        output3_dup,
                                                        output4_dup,
                                                        output5_dup,
                                                        output6_dup,
                                                        output7_dup,
                                                        output8_dup)),
                           }
        inspect_stop1 = atf.inspect_table(
                'localisation_test1',
                'stop1',
                data_dict_stop1,
                ('Localisation Recipe Called 5 Times or More While Producing '
                 'the Same Four Outputs in {0}'
                 ).format(reduced_path)
                )

else:
    color_stop1 = 'Red'
    result_stop1 = 'No'
    comment_stop1 = ('The number of unique output files should always be '
                     'equal or smaller than the number of recipe called.')
    inspect_stop1 = ''


# Build badpixel_test1.html

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
     <col span="1" style="width: 50%;">
     <col span="1" style="width: 20%;">
     <col span="1" style="width: 20%;">
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
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>2</td>
    <td># of outputs in {reduced_path}</td>
    <td>{output_list[0]}: {output1_num}<br>{output_list[1]}: {output2_num}<br>{output_list[2]}: {output3_num}<br>{output_list[3]}: {output4_num}<br>{output_list[4]}: {output5_num}<br>{output_list[5]}: {output6_num}<br>{output_list[6]}: {output7_num}<br>{output_list[7]}: {output8_num}</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>3</td>
    <td># of unique outputs in {reduced_path}</td>
    <td>{output_list[0]}: {output1_num_unique}<br>{output_list[1]}: {output2_num_unique}<br>{output_list[2]}: {output3_num_unique}<br>{output_list[3]}: {output4_num_unique}<br>{output_list[4]}: {output5_num_unique}<br>{output_list[5]}: {output6_num_unique}<br>{output_list[6]}: {output7_num_unique}<br>{output_list[7]}: {output8_num_unique}</td>
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
    <td># of entry in {reduced_path}/log.fits that failed to finish</td>
    <td>{num_logfits_ENDEDfalse}</td>
    <td>{comments_check5}</td>
    <td></td>
  </tr>
</table>


</body>
</html>
"""


with open('localisation_test1/localisation_test1.html', 'w') as f:
    f.write(html_text)
