import os
import glob
from datetime import datetime
from astropy.io import fits
import numpy as np

def intersection(lst1, lst2): 

    lst3 = [value for value in lst1 if value in lst2] 
    return lst3 


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


def list_odometers(path, files = 'all'):

    #files = 'all' => a, c, d, f, and o odometers are listed
    #files = 'a' => only a files
    #files = 'c' => only c files
    #files = 'd' => only d files
    #files = 'f' => only f files
    #files = 'o' => only o files

    files_list = os.listdir('{0}'.format(path))
    
    odometers = []

    if files == 'all':

        for i in range(len(files_list)):
            if 'a' in files_list[i][7]: odometers.append(files_list[i][:7])
            elif 'c' in files_list[i][7]: odometers.append(files_list[i][:7])
            elif 'd' in files_list[i][7]: odometers.append(files_list[i][:7])
            elif 'f' in files_list[i][7]: odometers.append(files_list[i][:7])
            elif 'o' in files_list[i][7]: odometers.append(files_list[i][:7])

    if files == 'a':

        for i in range(len(files_list)):
            if 'a' in files_list[i][7]: odometers.append(files_list[i][:7])

    if files == 'c':

        for i in range(len(files_list)):
            if 'c' in files_list[i][7]: odometers.append(files_list[i][:7])

    if files == 'd':

        for i in range(len(files_list)):
            if 'd' in files_list[i][7]: odometers.append(files_list[i][:7])

    if files == 'f':

        for i in range(len(files_list)):
            if 'f' in files_list[i][7]: odometers.append(files_list[i][:7])

    if files == 'o':

        for i in range(len(files_list)):
            if files_list[i][7] == 'o': odometers.append(files_list[i][:7])


    return odometers


class index_fits:

    def __init__(self, path):
        

        tbl = fits.getdata(path)
        self.tbl = tbl
        self.len = len(tbl)

        filename = tbl['FILENAME']
        self.filename= filename
        self.nights = tbl['NIGHTNAME']

        odometers = []

        for i in range(len(tbl)):
            index = filename[i].index('_pp.fits')
            odometers.append(filename[i][:index])

        self.odometers = np.array(odometers)


class log_fits:

    def __init__(self, path):


        tbl = fits.getdata(path)
        self.tbl = tbl
        self.len = len(tbl)

        self.QC = tbl['PASSED_ALL_QC']
        self.ENDED = tbl['ENDED']

        indexQCtrue = tbl['PASSED_ALL_QC'] == True
        indexQCfalse = tbl['PASSED_ALL_QC'] == False
        self.indexQCtrue = indexQCtrue
        self.indexQCfalse = indexQCfalse
        indexENDEDtrue = tbl['ENDED'] == True
        indexENDEDfalse = tbl['ENDED'] == False          
        self.indexENDEDtrue = indexENDEDtrue  
        self.indexENDEDfalse = indexENDEDfalse

        self.nights = tbl['DIRECTORY']

        args = tbl['ARGS']
        self.args = args
        self.QCstr = tbl['QC_STRING']
        self.ERRORS = tbl['ERRORS']
        self.LOGFILE = tbl['LOGFILE']

        odometers = []

        for i in range(len(args)):
            index = args[i].index('.fits')
            odometers.append(args[i][index-8:index])

        self.odometers = np.array(odometers)

from bokeh.io.saving import save
from bokeh.io.output import output_file
from bokeh.models import ColumnDataSource, DataTable, DateFormatter, TableColumn
from bokeh.models.widgets import Div
from bokeh.layouts import column

def inspect(test, check, data_dict, title):
    
    if not os.path.isdir('{0}/{1}'.format(test, check)):
        os.system('mkdir {0}/{1}'.format(test, check))

    source = ColumnDataSource(data_dict)
    
    keys_list = list(data_dict.keys())
    columns = []
    for i in range(len(keys_list)):
        columns.append(TableColumn(field = keys_list[i], title = keys_list[i]))
    
    table_title = Div(text="""<font size="+1"> <b>{0}</b> </font>""".format(title), height = 30)
    data_table = DataTable(source=source, columns=columns, autosize_mode = 'fit_columns', width=1500, height=500)
    layout = column(table_title, data_table)

    output_file("{0}/{1}/{1}.html".format(test, check), title="{0}".format(check))
    save(layout)

    html_str = """<a href='{0}/{0}.html'>Inspect""".format(check)

    return html_str

def summary(test):
    f = open("{0}/{0}.html".format(test), "r")
    html = f.read()

    passed = html.count('Lime')
    conditional = html.count('Yellow')
    failed = html.count('Red')

    n = passed + conditional + failed

    if passed == n :
        color = 'Lime'
    elif conditional >= 1 and failed == 0 :
        color = 'Yellow'
    else : 
        color = 'Red'

    html_str = """Passed: {0}/{3}<br>Passed with conditions: {1}/{3}<br>Failed: {2}/{3} """.format(passed, conditional, failed, n)

    return html_str, color
