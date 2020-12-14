import os
from apero.core import constants
from datetime import datetime
from apero_tests_func import *
import numpy as np

#ALL TESTS CONDUCT BY 'badpixel_test1.py'

#check1: how many recipes were run? (cal_badpix_{instrument} in log.fits)
#check2: how many of each output do we have?
         #output1: {ODOMETER_CODE}_pp_badpixel.fits
         #output2: {ODOMETER_CODE}_pp_bmap.fits
         #output3: DEBUG_{ODOMETER_CODE}_pp_background.fits
#check3: how many of each unique output do we have?
         #output1: {ODOMETER_CODE}_pp_badpixel.fits
         #output2: {ODOMETER_CODE}_pp_bmap.fits
         #output3: DEBUG_{ODOMETER_CODE}_pp_background.fits
#stop1: check3 == check1?
#check4: using the log.fits how many entry failed one or more QC?
#check5: using the log.fits how many entry failed to finish?
#check6: how many entry BADPIX and BKGRDMAP in master_calib_{INSTRUMENT}.txt? 
#check7: for each calib entry how many are in the calibDB?
#stop2: check7 == check6?


#constants

params = constants.load('SPIROU')

setup = os.environ['DRS_UCONFIG'] #setup
instrument = params['INSTRUMENT'] #instrument
date = datetime.now()
date = date.strftime("%Y-%m-%d %H:%M:%S") #date

reduced_path = params['DRS_DATA_REDUC']
if reduced_path[-1] == '/' :
    reduced_path = reduced_path[:-1] #reduced path without / at the end
reduced_nights = list_nights(reduced_path) #list reduced night directories

calibDB_path = params['DRS_CALIB_DB'] #calibDB path
if calibDB_path[-1] == '/' :
    calibDB_path = calibDB_path[:-1] #calibDB path without / at the end

#output list
output_list = ['*_pp_badpixel.fits', '*_pp_bmap.fits']

#calibDB entries list

calibDB_entry_list = ['BADPIX', 'BKGRDMAP']


#TESTS


#inspect all reduced_nights log.fits

recipe_num_logfits = 0 #check1

output1 = []
output2 = []

missing_logfits = []

for i in range(len(reduced_nights)):

    output1.extend(list_files('{0}/{1}'.format(reduced_path, reduced_nights[i]), files = '_pp_badpixel.fits'))
    output2.extend(list_files('{0}/{1}'.format(reduced_path, reduced_nights[i]), files = '_pp_bmap.fits'))

    #inspect log.fits if the file exists
    if os.path.isfile('{0}/{1}/log.fits'.format(reduced_path, reduced_nights[i])):
        logfits = log_fits('{0}/{1}/log.fits'.format(reduced_path, reduced_nights[i]))

        index_recipe = logfits.recipe == 'cal_badpix_{0}'.format(instrument.lower())
        recipe_num_logfits += sum(index_recipe)  #check1

        indexQCfalse = np.array(logfits.indexQCfalse) * np.array(index_recipe) #QC false + correct recipe
        indexENDEDfalse = np.array(logfits.indexENDEDfalse) * np.array(index_recipe) #ENDED false + correct recipe
        
        num_logfits_QCfalse = sum(indexQCfalse) #check4
        if num_logfits_QCfalse == 0:
            comments_check4 = ''
        else :
            comments_check4 = 'Inspect the log.fits of {0}/{1}. One or more entry have somehow failed QC.'.format(reduced_path, logfits.nights[0])

        num_logfits_ENDEDfalse = sum(indexENDEDfalse) #check5
        if num_logfits_ENDEDfalse == 0:
            comments_check5 = ''
        else :
            comments_check5 = 'Inspect the log.fits of {0}/{1}. One or more entry have somehow failed to finish.'.format(reduced_path, logfits.nights[0])

    #missing log.fits
    else:
        missing_logfits.append('{0}/{1}/log.fits'.format(reduced_path, reduced_nights[i]))

output1_num = len(output1) #check2
output1_num_unique = len(np.unique(output1)) #check3
output2_num = len(output2) #check2
output2_num_unique = len(np.unique(output2)) #check3

#stop1

if output1_num_unique == recipe_num_logfits and output2_num_unique == recipe_num_logfits:
    color1 = 'Lime'
    stop1 = 'Yes'
    comment1 = ''
    inspect1 = ''
elif output1_num_unique < recipe_num_logfits or output2_num_unique < recipe_num_logfits:
    color1 = 'Yellow'
    stop1 = 'No'
    comment1 = 'One or more output were not produced OR recipe called multiple times in the same directory.'
    inspect1 = ''
else :
    color1 = 'Red'
    stop1 = 'No'
    comment1 = 'The number of unique output files should always be equal or smaller than the number of recipe called.'
    inspect1 = ''

#inspect calibDB

f = open("{0}/master_calib_{1}.txt".format(calibDB_path, instrument), "r")
master_calib_txt = f.read()
index_start = master_calib_txt[:master_calib_txt.index('# DRS processed')].count('\n')
entry_col, night_col, file_col = np.genfromtxt("{0}/master_calib_{1}.txt".format(calibDB_path, instrument), delimiter = ' ', unpack = True, usecols = (0,2,3), skip_header = index_start, dtype=str)

output1_num_entry = master_calib_txt.count(calibDB_entry_list[0]) #check6
output1_calibDB = list_files("{0}".format(calibDB_path), files = '_pp_badpixel.fits')
output1_num_calibDB = len(output1_calibDB) #check7
output2_num_entry = master_calib_txt.count(calibDB_entry_list[1]) #check6
output2_calibDB = list_files("{0}".format(calibDB_path), files = '_pp_bmap.fits')
output2_num_calibDB = len(output1_calibDB) #check7

test, test1 = np.unique(file_col[entry_col == 'BADPIX'], return_counts=True)
dup = test[test1 > 1]
print(dup)
#stop2

if output1_num_calibDB == output1_num_entry and output2_num_calibDB == output2_num_entry:
    color2 = 'Lime'
    stop2 = 'Yes'
    comment2 = ''
    inspect2 = ''
elif output1_num_calibDB < output1_num_entry or output2_num_calibDB < output2_num_entry:
    color2 = 'Yellow'
    stop2 = 'No'
    comment2 = 'One or more output are not in the calibDB OR some entries are identical.'
    inspect2 = ''
else :
    color2 = 'Red'
    stop2 = 'No'
    comment2 = 'The calibDB should not have more output files than what was produced.'
    inspect2 = ''


#Build badpixel_test1.html
 
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

<h3>Bad Pixel Correction Recipe Test #1</h3>
<p><b>Setup: {setup}</b><br>
<p><b>Instrument: {instrument}</b><br>
<p><b>Date: {date}</b><br>
<br>
<p>Script: cal_badpix_{instrument.lower()}.py<br>
<p>Output files: {', '.join(output_list)}<br>
<p>Calibration database entry: {', '.join(calibDB_entry_list)}<br>
<p><a href='https://github.com/njcuk9999/apero-drs#83-bad-pixel-correction-recipe'>Link</a> to Bad Pixel Correction Recipe description</p>
<br></br>

<table id="t01">

  <colgroup>
     <col span="1" style="width: 5%;">
     <col span="1" style="width: 55%;">
     <col span="1" style="width: 15%;">
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
    <td># of time cal_badpix_{instrument.lower()}.py was called</td>
    <td>{recipe_num_logfits}</td>
    <td></td> 
    <td></td> 
  </tr>
  <tr>
    <td>2</td>
    <td># of {', '.join(output_list)} in {reduced_path}</td>
    <td>{output_list[0]}: {output1_num}<br>{output_list[1]}: {output2_num}</td>
    <td></td> 
    <td></td> 
  </tr>
  <tr>
    <td>3</td>
    <td># of unique {', '.join(output_list)} in {reduced_path}</td>
    <td>{output_list[0]}: {output1_num_unique}<br>{output_list[1]}: {output2_num_unique}</td>
    <td></td> 
    <td></td> 
  </tr>
  <tr>
    <td></td>
    <td>Check 3 == Check 1?</td>
    <td bgcolor={color1}>{stop1}</td>
    <td>{comment1}</td>
    <td>{inspect1}</td>
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
  <tr>
    <td>6</td>
    <td># of {', '.join(calibDB_entry_list)} entry in {calibDB_path}/master_calib_{instrument}.txt</td>
    <td>{calibDB_entry_list[0]}: {output1_num_entry}<br>{calibDB_entry_list[1]}: {output2_num_entry}</td>
    <td></td> 
    <td></td> 
  </tr>
  <tr>
    <td>7</td>
    <td># of {', '.join(output_list)} in {calibDB_path}</td>
    <td>{output_list[0]}: {output1_num_calibDB}<br>{output_list[1]}: {output2_num_calibDB}</td>
    <td></td> 
    <td></td> 
  </tr>
  <tr>
    <td></td>
    <td>Check 7== Check 6?</td>
    <td bgcolor={color2}>{stop2}</td>
    <td>{comment2}</td>
    <td>{inspect2}</td>
  </tr>
</table>


</body>
</html>
"""


with open('badpixel_test1/badpixel_test1.html', 'w') as f:
    f.write(html_text)
