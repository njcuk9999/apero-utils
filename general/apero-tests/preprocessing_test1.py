import os
from apero.core import constants
from datetime import datetime
from apero_tests_func import *
import numpy as np

#ALL TESTS CONDUCT BY 'preprocessing_test1.py'

#check1: how many raw files are there on disk?
#check2: how many pp files are there on disk?
#stop1: check2 == check1?
#check3: how many pp files are there in the index.fits?
#stop2: check3 == check2?
#check4: how many recipes were run? (cal_preprocess_{instrument} in tmp/*/log.fits)
#check5: how many unique odometers files were preprocessed according to the log.fits?
#stop3: check5 == check 2?
#check6: using the log.fits how many unique odometers failed one or more QC? Which odometer? Which QC?
#check7: using the log.fits how many unique odometers failed to finish? Which odometers? Why (using the ERRORS and LOGFILE columns)?


#constants

params = constants.load('SPIROU')

setup = os.environ['DRS_UCONFIG'] #setup
instrument = params['INSTRUMENT'] #instrument
date = datetime.now()
date = date.strftime("%Y-%m-%d %H:%M:%S") #date

raw_path = params['DRS_DATA_RAW']
if raw_path[-1] == '/' :
    raw_path = raw_path[:-1] #raw path without / at the end
raw_nights = list_nights(raw_path) #list raw night directories

pp_path = params['DRS_DATA_WORKING']
if pp_path[-1] == '/' :
    pp_path = pp_path[:-1] #pp path without / at the end
pp_nights = list_nights(pp_path) #list preprocessed data night directories

#output list
output_list = ['*_pp.fits']


#TESTS

raw_num = count_files_subdir(raw_path, subdir = 'all', files = '*.fits')  #check1
pp_num = count_files_subdir(pp_path, subdir = 'all', files = '*_pp.fits')    #check2

#stop1

if pp_num == raw_num:
    color1 = 'Lime'
    stop1 = 'Yes'
    comment1 = ''
    inspect1 = ''
elif pp_num < raw_num:
    color1 = 'Yellow'
    stop1 = 'No'
    comment1 = 'Not all available raw files were reduced.'
    inspect1 = ''
else :
    color1 = 'Red'
    stop1 = 'No'
    comment1 = 'The number of pp files should always be smaller than the number of raw files.'
    inspect1 = ''


#inspect all pp_nights index.fits and log.fits

pp_num_indexfits = 0 #check3
pp_num_logfits = 0 #check4
pp_num_logfits_unique = 0 #check5
pp_num_logfits_unique_QCfalse = 0 #check6
pp_num_logfits_unique_ENDEDfalse = 0 #check7

odometers_logfits_QCfalse = [] #check6
nights_logfits_QCfalse = [] #check6
QCstr_logfits_QCfalse = [] #check6

odometers_logfits_ENDEDfalse = [] #check7
nights_logfits_ENDEDfalse = [] #check7
ERRORS_logfits_ENDEDfalse = [] #check7
LOGFILE_logfits_ENDEDfalse = [] #check7

missing_indexfits = []
missing_logfits = []

#odometers_all = []
#test = []

for i in range(len(pp_nights)):

    #inspect index.fits if the file exists
    if os.path.isfile('{0}/{1}/index.fits'.format(pp_path, pp_nights[i])):
        
        indexfits = index_fits('{0}/{1}/index.fits'.format(pp_path, pp_nights[i]))
        pp_num_indexfits += indexfits.len #check3
    
    #missing index.fits
    else:
        missing_indexfits.append('{0}/{1}/index.fits'.format(pp_path, pp_nights[i]))


    #inspect log.fits if the file exists
    if os.path.isfile('{0}/{1}/log.fits'.format(pp_path, pp_nights[i])):
        
        logfits = log_fits('{0}/{1}/log.fits'.format(pp_path, pp_nights[i]))
        pp_num_logfits += logfits.len #check4

        #don't consider duplicates pp files
        args = logfits.args
        
        odometers = []

        for j in range(len(args)):
            if 'persi_' in args[j]: args[j].replace('persi_','')
            index = args[j].index('.fits')
            odometers.append(args[j][index-8:index])

        odometers, index_unique = np.unique(odometers, return_index = True)

        tbl_unique = logfits.tbl[index_unique]
        pp_num_logfits_unique += len(tbl_unique) #check5

        indexQCfalse = tbl_unique['PASSED_ALL_QC'] == False
        indexENDEDfalse = tbl_unique['ENDED'] == False

        pp_num_logfits_unique_QCfalse += sum(indexQCfalse) #check6
        odometers_logfits_QCfalse.extend(odometers[indexQCfalse]) #check6
        nights_logfits_QCfalse.extend(tbl_unique['DIRECTORY'][indexQCfalse]) #check6
        QCstr_logfits_QCfalse.extend(tbl_unique['QC_STRING'][indexQCfalse]) #check6

        pp_num_logfits_unique_ENDEDfalse += sum(indexENDEDfalse) #check7
        odometers_logfits_ENDEDfalse.extend(odometers[indexENDEDfalse]) #check7
        nights_logfits_ENDEDfalse.extend(tbl_unique['DIRECTORY'][indexENDEDfalse]) #check7
        ERRORS_logfits_ENDEDfalse.extend(tbl_unique['ERRORS'][indexENDEDfalse]) #check7
        LOGFILE_logfits_ENDEDfalse.extend(tbl_unique['LOGFILE'][indexENDEDfalse]) #check7

    #missing log.fits
    else:
        missing_logfits.append('{0}/{1}/log.fits'.format(pp_path, pp_nights[i]))
  

    #odometers_all.extend(odometers)
    #file_name = list_files('{0}/{1}'.format(pp_path, pp_nights[i]), files = '_pp.fits')
    #test.extend(file_name[:8])
    

#print(len(odometers_all))
#print(len(test))
#print(len(np.intersect1d(odometers_all, test)))
#print(np.setdiff1d(odometers_all, test))

#stop2

if pp_num_indexfits == pp_num:
    color2 = 'Lime'
    stop2 = 'Yes'
    comment2 = ''
    inspect2 = ''
    
else:
    color2 = 'Red'
    stop2 = 'No'
    comment2 = ''
    inspect2 = ''

#stop3

if pp_num_logfits_unique == pp_num:
    color3 = 'Lime'
    stop3 = 'Yes'
    comment3 = ''
    inspect3 = ''

elif pp_num_logfits_unique > pp_num:
    color3 = 'Yellow'
    stop3 = 'No'
    comment3 = 'Some files were processed more than once.'
    inspect3 = ''
    
else:
    color3 = 'Red'
    stop3 = 'No'
    comment3 = ''
    inspect3 = ''


data_dict_check6 = {'Night': nights_logfits_QCfalse,
             'Odometer': odometers_logfits_QCfalse,
             'QC_STRING': QCstr_logfits_QCfalse,
}
inspect_check6 = inspect('preprocessing_test1', 'check6', data_dict_check6, 'Odometers that Failed One or More Quality Control')

data_dict_check7 = {'Night': nights_logfits_ENDEDfalse,
             'Odometer': odometers_logfits_ENDEDfalse,
             'ERRORS': ERRORS_logfits_ENDEDfalse,
             'LOGFILE': LOGFILE_logfits_ENDEDfalse,
}
inspect_check7 = inspect('preprocessing_test1', 'check7', data_dict_check7, 'Odometers that Failed to Finish')


#Build preprocessing_test1.html
 
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

<h3>Preprocessing Recipe Test #1</h3>
<p><b>Setup: {setup}</b><br>
<p><b>Instrument: {instrument}</b><br>
<p><b>Date: {date}</b><br>
<br>
<p>Script: cal_preprocessing_{instrument.lower()}.py<br>
<p>Output files: {output_list}<br>
<p><a href='https://github.com/njcuk9999/apero-drs#81-preprocessing-recipe'>Link</a> to Preprocessing Recipe description</p>
<br></br>

<table id="t01">

  <colgroup>
     <col span="1" style="width: 5%;">
     <col span="1" style="width: 55%;">
     <col span="1" style="width: 5%;">
     <col span="1" style="width: 30%;">
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
    <td># of raw files in {raw_path}</td>
    <td>{raw_num}</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>2</td>
    <td># of pp files in {pp_path}</td>
    <td>{pp_num}</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td> </td>
    <td>Check 2 == Check 1?</td>
    <td bgcolor={color1}>{stop1}</td>
    <td>{comment1}</td>
    <td>{inspect1}</td>
  </tr>
  <tr>
    <td>3</td>
    <td># of pp files in {pp_path}/*/index.fits </td>
    <td>{pp_num_indexfits}</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td> </td>
    <td>Check 3 == Check 2?</td>
    <td bgcolor={color2}>{stop2}</td>
    <td>{comment2}</td>
    <td>{inspect2}</td>
  </tr>
  <tr>
    <td>4</td>
    <td># of time cal_preprocess_{instrument.lower()}.py was called</td>
    <td>{pp_num_logfits}</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>5</td>
    <td># of unique odometers preprocessed according to {pp_path}/*/log.fits</td>
    <td>{pp_num_logfits_unique}</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td> </td>
    <td>Check 5 == Check 2?</td>
    <td bgcolor={color3}>{stop3}</td>
    <td>{comment3}</td>
    <td>{inspect3}</td>
  </tr>
  <tr>
    <td>6</td>
    <td># of unique odometers in {pp_path}/*/log.fits that failed one or more QC</td>
    <td>{pp_num_logfits_unique_QCfalse}</td>
    <td></td>
    <td>{inspect_check6}</td>
  </tr>
  <tr>
    <td>7</td>
    <td># of unique odometers in {pp_path}/*/log.fits that failed to finish</td>
    <td>{pp_num_logfits_unique_ENDEDfalse}</td>
    <td></td>
    <td>{inspect_check7}</td>
  </tr>
</table>


</body>
</html>
"""


with open('preprocessing_test1/preprocessing_test1.html', 'w') as f:
    f.write(html_text)
