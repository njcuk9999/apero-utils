import os
from apero.core import constants
from datetime import datetime
from apero_tests_func import *

#test1: how many raw files are there on disk
#test2: how many pp files are there on disk (compared to raw files)
#test3: how many pp files are there in the index.fits
#test4: using the log.fits how many passed all QC
#test5: using the log.fits how many failed one or more QC
#test6: which odometer failed one or more QC
#test7: which QC did they fail
#test8: using the log.fits how many failed to finish
#test9: which odometer failed to finish
#test10: why did they failed to finish using the 'ERRORS' columns
#test11: why did they failed to finish using the 'LOGFILE' columns

params = constants.load('SPIROU')

setup = os.environ['DRS_UCONFIG'] #setup
instrument = params['INSTRUMENT'] #instrument
date = datetime.now()
date = date.strftime("%Y-%m-%d %H:%M:%S") #date

raw_path = params['DRS_DATA_RAW']
raw_nights = list_nights(raw_path)

pp_path = params['DRS_DATA_WORKING']
pp_nights = list_nights(pp_path)

raw_num = count_files_subdir(raw_path, subdir = 'all', files = '*.fits')  #test1
pp_num = count_files_subdir(pp_path, subdir = 'all', files = '*.fits')    #test2

pp_indexfits_num = 0
pp_logfits_QCpassed_num = 0
pp_logfits_QCfailed_num = 0
pp_logfits_endedfailed_num = 0

odometer_QC = []
NIGHTNAME_QC = []
QCstr = []
odometer_endedfailed = []
NIGHTNAME_endedfailed = []
ERRORS = []
LOGFILE = []

for i in range(len(pp_nights)):

    indexfits = index_fits('{0}/{1}/index.fits'.format(pp_path, pp_nights[i]))
    pp_indexfits_num += indexfits.len #test3
    
    logfits = log_fits('{0}/{1}/log.fits'.format(pp_path, pp_nights[i]))

    pp_logfits_QCpassed_num += logfits.lenQCpass #test4
    pp_logfits_QCfailed_num += logfits.lenQCfail #test5

    odometer_QC.extend(logfits.odometerQCfail) #test6
    NIGHTNAME_QC.extend(logfits.nightQCfail)
    QCstr.extend(logfits.QCstrfail) #test7

    pp_logfits_endedfailed_num += logfits.lenEndfail #test8
    odometer_endedfailed.extend(logfits.odometerEndfail) #test9
    NIGHTNAME_endedfailed.extend(logfits.nightEndfail)
    ERRORS.extend(logfits.errorEndfail) #test10
    LOGFILE.extend(logfits.logfileEndfail) #test11



html_hover_QC = []
for i in range(len(odometer_QC)):
    html_hover_QC.append("""<a href=" " title="NIGHTNAME: {0}\nQC_STRING: {1}" style="color:#000000;text-decoration:none">{2}</a>""".format(NIGHTNAME_QC[i], QCstr[i], odometer_QC[i]))
html_hover_QC = ", ".join(html_hover_QC)


html_hover_endedfailed = []
for i in range(len(odometer_endedfailed)):
    html_hover_endedfailed.append("""<a href=" " title="NIGHTNAME: {0}\nERRORS: {1}\nLOGFILE: {2}" style="color:#000000;text-decoration:none">{3}</a>""".format(NIGHTNAME_endedfailed[i], ERRORS[i], LOGFILE[i], odometer_endedfailed[i]))
html_hover_endedfailed = ", ".join(html_hover_endedfailed)


 
html_text = f"""
<html>



<head>
<title>APERO Tests</title>
<style>
table {{
  width:100%;
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
     <col span="1" style="width: 50%;">
     <col span="1" style="width: 50%;">
  </colgroup>
  <tr>
    <th>Test</th>
    <th>Result</th> 
  </tr>
  <tr>
    <td># of raw files in {raw_path}</td>
    <td>{raw_num}</td>
  </tr>
  <tr>
    <td># of pp files in {pp_path}</td>
    <td>{pp_num}</td>
  </tr>
  <tr>
    <td># of pp files in {pp_path}/*/index.fits</td>
    <td>{pp_indexfits_num}</td>
  </tr>
  <tr>
    <td># of pp files in {pp_path}/*/log.fits that passed all QC</td>
    <td>{pp_logfits_QCpassed_num}</td>
  </tr>
  <tr>
    <td># of pp files in {pp_path}/*/log.fits that failed one or more QC</td>
    <td>{pp_logfits_QCfailed_num}</td>
  </tr>
  <tr>
    <td>pp files that failed one or more QC, and which QC?</td>
    <td><small>{html_hover_QC}</small></td>
  </tr>
  <tr>
    <td># of pp files in {pp_path}/*/log.fits that failed to finish</td>
    <td>{pp_logfits_endedfailed_num}</td>
  </tr>
  <tr>
    <td>pp files that failed to finish, and why?</td>
    <td><small>{html_hover_endedfailed}</small></td>
  </tr>
</table>


</body>
</html>
"""


with open('preprocessing_test1/preprocessing_test1.html', 'w') as f:
    f.write(html_text)
