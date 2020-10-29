import os
from apero.core import constants
from datetime import datetime
from apero_tests_func import *

#check1: how many raw files are there on disk
#check2: how many pp files are there on disk (compared to raw files)

#stop1: check2 == check1?

#check3: how many pp files are there in the index.fits

#stop2: check3 == check2?

#check4: using the log.fits how many passed all QC
#check5: using the log.fits how many failed one or more QC
#check6: which odometer failed one or more QC
#check7: which QC did they fail
#check8: using the log.fits how many failed to finish
#check9: which odometer failed to finish
#check10: why did they failed to finish using the 'ERRORS' columns
#check11: why did they failed to finish using the 'LOGFILE' columns

params = constants.load('SPIROU')

setup = os.environ['DRS_UCONFIG'] #setup
instrument = params['INSTRUMENT'] #instrument
date = datetime.now()
date = date.strftime("%Y-%m-%d %H:%M:%S") #date

raw_path = params['DRS_DATA_RAW']
raw_nights = list_nights(raw_path)

pp_path = params['DRS_DATA_WORKING']
pp_nights = list_nights(pp_path)

raw_num = count_files_subdir(raw_path, subdir = 'all', files = '*.fits')  #check1
pp_num = count_files_subdir(pp_path, subdir = 'all', files = '*pp.fits')    #check2

#stop1

if pp_num == raw_num:
    color1 = 'Lime'
    stop1 = 'Yes'
    comment1 = ''
elif pp_num < raw_num:
    color1 = 'Yellow'
    stop1 = 'No'
    comment1 = 'Not all available raw files were reduced'
else :
    color1 = 'Red'
    stop1 = 'No'
    comment1 = 'The number of pp files should always be smaller than the number of raw files'


pp_indexfits_num = 0
pp_logfits_QCpassed_num = 0
pp_logfits_QCfailed_num = 0
pp_logfits_endedfailed_num = 0

odometer_index = []
odometer_QC = []
NIGHTNAME_QC = []
QCstr = []
odometer_endedfailed = []
NIGHTNAME_endedfailed = []
ERRORS = []
LOGFILE = []

for i in range(len(pp_nights)):

    indexfits = index_fits('{0}/{1}/index.fits'.format(pp_path, pp_nights[i]))
    odometer_index.extend(indexfits.odometer)

    pp_indexfits_num += indexfits.len #check3
    
    logfits = log_fits('{0}/{1}/log.fits'.format(pp_path, pp_nights[i]))

    pp_logfits_QCpassed_num += logfits.lenQCpass #check4
    pp_logfits_QCfailed_num += logfits.lenQCfail #check5

    odometer_QC.extend(logfits.odometerQCfail) #check6
    NIGHTNAME_QC.extend(logfits.nightQCfail)
    QCstr.extend(logfits.QCstrfail) #check7

    pp_logfits_endedfailed_num += logfits.lenEndfail #check8
    odometer_endedfailed.extend(logfits.odometerEndfail) #check9
    NIGHTNAME_endedfailed.extend(logfits.nightEndfail)
    ERRORS.extend(logfits.errorEndfail) #check10
    LOGFILE.extend(logfits.logfileEndfail) #check11

#stop2

if pp_indexfits_num == pp_num:
    color2 = 'Lime'
    stop2 = 'Yes'
    comment2 = """<a href=" ">Inspect</a>"""
    
else:
    color2 = 'Red'
    stop2 = 'No'
    comment2 = """<a href=" ">Inspect</a>"""



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
  width:70%;
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
     <col span="1" style="width: 60%;">
     <col span="1" style="width: 5%;">
     <col span="1" style="width: 30%;">
  </colgroup>
  <tr>
    <th>Check</th>
    <th>Description</th>
    <th>Result</th>
    <th>Comment</th>  
  </tr>
  <tr>
    <td>1</td>
    <td># of raw files in {raw_path}</td>
    <td>{raw_num}</td>
    <td></td>
  </tr>
  <tr>
    <td>2</td>
    <td># of pp files in {pp_path}</td>
    <td>{pp_num}</td>
    <td></td>
  </tr>
  <tr>
    <td> </td>
    <td>Check 2 == Check 1?</td>
    <td bgcolor={color1}>{stop1}</td>
    <td>{comment1}</td>
  </tr>
  <tr>
    <td>3</td>
    <td># of pp files in {pp_path} index.fits </td>
    <td>{pp_indexfits_num}</td>
    <td></td>
  </tr>
  <tr>
    <td> </td>
    <td>Check 3 == Check 2?</td>
    <td bgcolor={color2}>{stop2}</td>
    <td>{comment2}</td>
  </tr>
  <tr>
    <td>4</td>
    <td># of pp files in {pp_path}/*/log.fits that passed all QC</td>
    <td>{pp_logfits_QCpassed_num}</td>
    <td></td>
  </tr>
  <tr>
    <td>5</td>
    <td># of pp files in {pp_path}/*/log.fits that failed one or more QC</td>
    <td>{pp_logfits_QCfailed_num}</td>
    <td></td>
  </tr>
  <tr>
    <td>6</td>
    <td>pp files that failed one or more QC, and which QC?</td>
    <td><small>{html_hover_QC}</small></td>
    <td></td>
  </tr>
  <tr>
    <td>7</td>
    <td># of pp files in {pp_path}/*/log.fits that failed to finish</td>
    <td>{pp_logfits_endedfailed_num}</td>
    <td></td>
  </tr>
  <tr>
    <td>8</td>
    <td>pp files that failed to finish, and why?</td>
    <td><small>{html_hover_endedfailed}</small></td>
    <td></td>
  </tr>
</table>


</body>
</html>
"""


with open('preprocessing_test1/preprocessing_test1.html', 'w') as f:
    f.write(html_text)
