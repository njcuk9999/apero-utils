import os
from apero.core import constants
from datetime import datetime
from apero_tests_func import *
import numpy as np

#ALL TESTS CONDUCT BY 'badpixel_test1.py'

#check1: how many recipes were run?
#check2: how many of each output do we have?
#check3: using the log.fits how many failed one or more QC? Which odometer? Which QC?
#check4: for each calib entry how many are in the calibDB?
#check5: plot the different QCs as a function of time
#check6: using the log.fits how many failed to finish? Which odometers? Why (using the ERRORS and LOGFILE columns)?


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




#inspect all reduced_nights log.fits

recipe_num_logfits = 0 #check1

for i in range(len(reduced_nights)):

    #inspect index.fits if the file exists
    if os.path.isfile('{0}/{1}/index.fits'.format(reduced_path, reduced_nights[i])):
        
    
    #missing index.fits
    else:
        missing_indexfits.append('{0}/{1}/index.fits'.format(reduced_path, reduced_nights[i]))


    if os.path.isfile('{0}/{1}/log.fits'.format(reduced_path, reduced_nights[i])):
        logfits = log_fits('{0}/{1}/log.fits'.format(reduced_path, reduced_nights[i]))

        recipe_num_logfits += sum(logfits.recipe == 'cal_badpix_spirou')  #check1





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
    <td>How many recipes were run</td>
    <td>{recipe_num_logfits}</td>
    <td></td> 
    <td></td> 
  </tr>
  <tr>
    <td>2</td>
    <td>How many of each output do we have</td>
    <td></td>
    <td></td> 
    <td></td> 
  </tr>
  <tr>
    <td>3</td>
    <td># of unique files in {reduced_path}/*/log.fits that failed one or more QC</td>
    <td></td>
    <td></td> 
    <td></td> 
  </tr>
  <tr>
    <td>4</td>
    <td>For each calib entry how many are in the calibDB</td>
    <td></td>
    <td></td> 
    <td></td> 
  </tr>
  <tr>
    <td>5</td>
    <td>Plot the different QCs as a function of time</td>
    <td></td>
    <td></td> 
    <td></td> 
  </tr>
  <tr>
    <td>6</td>
    <td># of unique files in {reduced_path}/*/log.fits that failed to finish</td>
    <td></td>
    <td></td> 
    <td></td>
  </tr>
</table>


</body>
</html>
"""


with open('badpixel_test1/badpixel_test1.html', 'w') as f:
    f.write(html_text)
