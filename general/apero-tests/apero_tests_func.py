import os
import glob
from datetime import datetime
from astropy.io import fits


def list_nights(path):

    list = [x for x in os.listdir(path) if os.path.isdir(os.path.join(path, x))]
    list.sort(key=lambda date: datetime.strptime(date[:10], "%Y-%m-%d"))
    return list


def count_files(path, files = 'all'):

    if files == 'all':
        f = '*'
    else:
        f = files

    return len(glob.glob('{0}/{1}'.format(path, f)))


def count_files_subdir(path, subdir='all', files='all'):

    if subdir == 'all':
        s = '*'
    else:
        s = subdir

    if files == 'all':
        f = '*'
    else:
        f = files

    return len(glob.glob('{0}/{1}/{2}'.format(path, s, f)))


class index_fits:

    def __init__(self, path):

        self.tbl = fits.getdata(path)
        self.len = len(fits.getdata(path))


class log_fits:

    def __init__(self, path):

        tbl = fits.getdata(path)
        self.tbl = tbl
        self.len = len(tbl)
        self.lenQCpass = sum(tbl['PASSED_ALL_QC'])
        self.lenQCfail = len(tbl) - sum(tbl['PASSED_ALL_QC'])
        self.lenEndpass = sum(tbl['ENDED'])
        self.lenEndfail = len(tbl) - sum(tbl['ENDED'])

        indexQCpass = tbl['PASSED_ALL_QC'] == True
        indexQCfail = tbl['PASSED_ALL_QC'] == False
        self.indexQCpass = indexQCpass
        self.indexQCfail = indexQCfail
        indexEndpass = tbl['ENDED'] == True
        indexEndfail = tbl['ENDED'] == False          
        self.indexEndpass = indexEndpass  
        self.indexEndfail = indexEndfail

        self.nightQCpass = tbl['DIRECTORY'][indexQCpass]
        self.nightQCfail = tbl['DIRECTORY'][indexQCfail]
        self.nightEndpass = tbl['DIRECTORY'][indexEndpass]
        self.nightEndfail = tbl['DIRECTORY'][indexEndfail]
        self.QCstrpass = tbl['QC_STRING'][indexQCpass]
        self.QCstrfail = tbl['QC_STRING'][indexQCfail]
        self.errorEndpass = tbl['ERRORS'][indexEndpass]
        self.errorEndfail = tbl['ERRORS'][indexEndfail]
        self.logfileEndpass = tbl['LOGFILE'][indexEndpass]
        self.logfileEndfail = tbl['LOGFILE'][indexEndfail]

        args = tbl['ARGS']
        argsQCpass = tbl['ARGS'][indexQCpass]
        argsQCfail = tbl['ARGS'][indexQCfail]
        argsEndpass = tbl['ARGS'][indexEndpass]
        argsEndfail = tbl['ARGS'][indexEndfail]

        self.args = args      
        self.argsQCpass = argsQCpass
        self.argsQCfail = argsQCfail
        self.argsEndpass = argsEndpass
        self.argsEndfail = argsEndfail

        odometer = []

        for i in range(len(args)):
            index = args[i].index('.fits')
            odometer.append(args[i][index-8:index])

        self.odometer = odometer
        self.odometerQCpass = [x for x, y in zip(odometer, indexQCpass) if y == True]
        self.odometerQCfail = [x for x, y in zip(odometer, indexQCfail) if y == True]
        self.odometerEndpass = [x for x, y in zip(odometer, indexEndpass) if y == True]
        self.odometerEndfail = [x for x, y in zip(odometer, indexEndfail) if y == True]
