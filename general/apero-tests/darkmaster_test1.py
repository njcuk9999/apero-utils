import os
from apero.core import constants
from datetime import datetime
from apero_tests_func import *

#ALL TESTS CONDUCT BY 'darkmaster_test1.py'

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
