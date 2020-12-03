import os
from apero.core import constants
from datetime import datetime
from apero_tests_func import *

#ALL TESTS CONDUCT BY 'preprocessing_test1.py'

#check1: how many raw files are there on disk?
#check2: how many pp files are there on disk?
#stop1: check2 == check1?
#check3: how many pp files are there in the index.fits?
#stop2: check3 == check2?
#check4: how many pp files are there in the log.fits?
#check5: how many unique pp files are there in the log.fits?
#stop3: check5 == check 2?
#check6: using the log.fits how many failed one or more QC? Which odometer? Which QC?
#check7: using the log.fits how many failed to finish? Which odometers? Why (using the ERRORS and LOGFILE columns)?


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


#TESTS

raw_num = count_files_subdir(raw_path, subdir = 'all', files = '*.fits')  #check1
pp_num = count_files_subdir(pp_path, subdir = 'all', files = '*pp.fits')    #check2

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
pp_num_unique_logfits = 0 #check5
pp_num_logfits_QC = 0
pp_num_logfits_QCtrue = 0
pp_num_logfits_ENDED = 0
pp_num_logfits_ENDEDtrue = 0 

odometers_logfits_QCfalse = []
nights_logfits_QCfalse = []
QCstr_logfits_QCfalse = []

odometers_logfits_ENDEDfalse = []
nights_logfits_ENDEDfalse = []
ERRORS_logfits_ENDEDfalse = []
LOGFILE_logfits_ENDEDfalse = []

missing_indexfits = []
missing_logfits = []

for i in range(len(pp_nights)):

    if os.path.isfile('{0}/{1}/index.fits'.format(pp_path, pp_nights[i])):
        indexfits = index_fits('{0}/{1}/index.fits'.format(pp_path, pp_nights[i]))
        pp_num_indexfits += indexfits.len #check3

    else:
        missing_indexfits.append('{0}/{1}/index.fits'.format(pp_path, pp_nights[i]))


    if os.path.isfile('{0}/{1}/log.fits'.format(pp_path, pp_nights[i])):
        logfits = log_fits('{0}/{1}/log.fits'.format(pp_path, pp_nights[i]))

        pp_num_logfits += logfits.len #check4
        pp_num_unique_logfits += logfits.len_unique #check5
        pp_num_logfits_QC += len(logfits.QC)
        pp_num_logfits_QCtrue += sum(logfits.QC)
        pp_num_logfits_ENDED += len(logfits.ENDED)
        pp_num_logfits_ENDEDtrue += sum(logfits.ENDED)

        indexQCfalse = logfits.indexQCfalse
        indexENDEDfalse = logfits.indexENDEDfalse

        odometers_logfits_QCfalse.extend(logfits.odometers[indexQCfalse]) #check6
        nights_logfits_QCfalse.extend(logfits.nights[indexQCfalse]) #check6
        QCstr_logfits_QCfalse.extend(logfits.QCstr[indexQCfalse]) #check6

        odometers_logfits_ENDEDfalse.extend(logfits.odometers[indexENDEDfalse]) #check7
        nights_logfits_ENDEDfalse.extend(logfits.nights[indexENDEDfalse]) #check7
        ERRORS_logfits_ENDEDfalse.extend(logfits.ERRORS[indexENDEDfalse]) #check7
        LOGFILE_logfits_ENDEDfalse.extend(logfits.LOGFILE[indexENDEDfalse]) #check7


    else:
        missing_logfits.append('{0}/{1}/log.fits'.format(pp_path, pp_nights[i]))

pp_num_logfits_QCfalse = pp_num_logfits_QC - pp_num_logfits_QCtrue #check6
pp_num_logfits_ENDEDfalse = pp_num_logfits_ENDED - pp_num_logfits_ENDEDtrue #check7


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

if pp_num_unique_logfits == pp_num:
    color3 = 'Lime'
    stop3 = 'Yes'
    comment3 = ''
    inspect3 = ''

elif pp_num_unique_logfits > pp_num:
    color3 = 'Yellow'
    stop3 = 'No'
    comment3 = 'Some files were processed more than once.'
    inspect3 = ''
    
else:
    color3 = 'Red'
    stop3 = 'No'
    comment3 = ''
    inspect3 = ''


data_dict6 = {'Night': nights_logfits_QCfalse,
             'Odometer': odometers_logfits_QCfalse,
             'QC_STRING': QCstr_logfits_QCfalse,
}
inspect6 = inspect('preprocessing_test1', 'check7', data_dict6, 'Odometers that failed one or more Quality Control')

data_dict7 = {'Night': nights_logfits_ENDEDfalse,
             'Odometer': odometers_logfits_ENDEDfalse,
             'ERRORS': ERRORS_logfits_ENDEDfalse,
             'LOGFILE': LOGFILE_logfits_ENDEDfalse,
}
inspect7 = inspect('preprocessing_test1', 'check8', data_dict7, 'Odometers that failed to finish')


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

<h3>Preprocessing Test #1</h3>
<p><b>Setup: {setup}</b></p>
<p><b>Instrument: {instrument}</b></p>
<p><b>Date: {date}</b></p>
<p>   </p>

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
    <td></td>
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
    <td></td>
  </tr>
  <tr>
    <td>4</td>
    <td># of pp files in {pp_path}/*/log.fits</td>
    <td>{pp_num_logfits}</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>5</td>
    <td># of unique pp files in {pp_path}/*/log.fits</td>
    <td>{pp_num_unique_logfits}</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td> </td>
    <td>Check 5 == Check 2?</td>
    <td bgcolor={color3}>{stop3}</td>
    <td>{comment3}</td>
    <td></td>
  </tr>
  <tr>
    <td>6</td>
    <td># of unique pp files in {pp_path}/*/log.fits that failed one or more QC</td>
    <td>{pp_num_logfits_QCfalse}</td>
    <td></td>
    <td>{inspect6}</td>
  </tr>
  <tr>
    <td>7</td>
    <td># of unique pp files in {pp_path}/*/log.fits that failed to finish</td>
    <td>{pp_num_logfits_ENDEDfalse}</td>
    <td></td>
    <td>{inspect7}</td>
  </tr>
</table>


</body>
</html>
"""


with open('preprocessing_test1/preprocessing_test1.html', 'w') as f:
    f.write(html_text)
