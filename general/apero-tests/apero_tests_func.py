"""
Utility functions for APERO tests

@author: charles
"""
import os
import glob
from datetime import datetime
from astropy.io import fits
from bokeh.io.saving import save
from bokeh.io.output import output_file
from bokeh.models import ColumnDataSource, DataTable, DateFormatter, TableColumn
from bokeh.models.widgets import Div
from bokeh.layouts import column


def intersection(lst1, lst2):
    """
    Returns interesection of list
    """

    # TODO: can probably use numpy for this
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def list_nights(path):
    """
    Return a list of datetime objects corresponding to each night
    """

    nlist = [x
             for x in os.listdir(path)
             if os.path.isdir(os.path.join(path, x))
             ]
    if 'other' in nlist:
        nlist.remove('other')
    nlist.sort(key=lambda date: datetime.strptime(date[:10], "%Y-%m-%d"))

    return list_nights


def list_files(path, files='all'):
    """
    Get list of files satisfying a pattern.
    """

    if files == 'all':
        f = ''
    else:
        f = files

    flist = [x for x in os.listdir(path) if x.endswith(f)]

    if len(flist) == 0:
        print(path)
    return flist


def count_files(path, files='all'):
    """
    Count number of files of a give pattern in a directory.
    """

    if files == 'all':
        f = '*'
    else:
        f = files

    return len(glob.glob('{0}/{1}'.format(path, f)))


def count_files_subdir(path, subdir='all', files='all'):
    """
    Count files in selected subdirectories of a given dir.
    """

    if subdir == 'all':
        s = '*'
    else:
        s = subdir

    if files == 'all':
        f = '*'
    else:
        f = files

    return len(glob.glob('{0}/{1}/{2}'.format(path, s, f)))


def list_raw_odometers(path, files='all'):
    """
    Get a list of raw data odometers
    """

    # TODO: Simplify structure
    # files = 'all' => a, c, d, f, and o odometers are listed
    # files = 'a' => only a files
    # files = 'c' => only c files
    # files = 'd' => only d files
    # files = 'f' => only f files
    # files = 'o' => only o files

    files_list = os.listdir('{0}'.format(path))

    odometers = []

    if files == 'all':

        for i in range(len(files_list)):

            if 'persi_' in files_list[i]:
                files_list[i].replace('persi_', '')

            if 'a' in files_list[i][7]:
                odometers.append(files_list[i][:8])
            elif 'c' in files_list[i][7]:
                odometers.append(files_list[i][:8])
            elif 'd' in files_list[i][7]:
                odometers.append(files_list[i][:8])
            elif 'f' in files_list[i][7]:
                odometers.append(files_list[i][:8])
            elif 'o' in files_list[i][7]:
                odometers.append(files_list[i][:8])

    if files == 'a':

        for i in range(len(files_list)):

            if 'persi_' in files_list[i]:
                files_list[i].replace('persi_', '')

            if 'a' in files_list[i][7]:
                odometers.append(files_list[i][:8])

    if files == 'c':

        for i in range(len(files_list)):

            if 'persi_' in files_list[i]:
                files_list[i].replace('persi_', '')
            if 'c' in files_list[i][7]:
                odometers.append(files_list[i][:8])

    if files == 'd':

        for i in range(len(files_list)):

            if 'persi_' in files_list[i]:
                files_list[i].replace('persi_', '')

            if 'd' in files_list[i][7]:
                odometers.append(files_list[i][:8])

    if files == 'f':

        for i in range(len(files_list)):

            if 'persi_' in files_list[i]:
                files_list[i].replace('persi_', '')

            if 'f' in files_list[i][7]:
                odometers.append(files_list[i][:8])

    if files == 'o':

        for i in range(len(files_list)):

            if 'persi_' in files_list[i]:
                files_list[i].replace('persi_', '')

            if files_list[i][7] == 'o':
                odometers.append(files_list[i][:8])


    return odometers


class index_fits:

    def __init__(self, path):

        tbl = fits.getdata(path)
        self.tbl = tbl
        self.len = len(tbl)
        self.filename = tbl['FILENAME']
        self.nights = tbl['NIGHTNAME']
        self.object = tbl['KW_OBJNAME']


class log_fits:

    def __init__(self, path):

        tbl = fits.getdata(path)
        self.tbl = tbl
        self.len = len(tbl)
        self.recipe = tbl['RECIPE']
        self.QC = tbl['PASSED_ALL_QC']
        self.ENDED = tbl['ENDED']

        self.indexQCtrue = tbl['PASSED_ALL_QC']
        self.indexQCfalse = ~tbl['PASSED_ALL_QC']
        self.indexENDEDtrue = tbl['ENDED']
        self.indexENDEDfalse = ~tbl['ENDED']

        self.nights = tbl['DIRECTORY']
        self.args = tbl['ARGS']

        self.QCstr = tbl['QC_STRING']
        self.ERRORS = tbl['ERRORS']
        self.LOGFILE = tbl['LOGFILE']


def inspect(test, check, data_dict, title):
    """
    Write an 'inspect' html element.
    """

    if not os.path.isdir('{0}/{1}'.format(test, check)):
        os.system('mkdir {0}/{1}'.format(test, check))

    source = ColumnDataSource(data_dict)

    keys_list = list(data_dict.keys())
    columns = []
    for i in range(len(keys_list)):
        columns.append(TableColumn(field=keys_list[i], title=keys_list[i]))

    table_title = Div(
            text="""<font size="+1"> <b>{0}</b> </font>""".format(title),
            height=50
            )
    data_table = DataTable(
            source=source,
            columns=columns,
            autosize_mode='fit_columns',
            min_width=1000,
            max_width=1500,
            width_policy='min',
            height=600,
            editable=True)
    layout = column(table_title, data_table)

    output_file(
            "{0}/{1}/{1}.html".format(test, check),
            title="{0}".format(check))
    save(layout)

    html_str = """<a href='{0}/{0}.html'>Inspect</a>""".format(check)

    return html_str


def summary(test):
    """
    Write the html summary.
    """

    f = open("{0}/{0}.html".format(test), "r")
    html = f.read()

    passed = html.count('Lime')
    conditional = html.count('Yellow')
    failed = html.count('Red')

    n = passed + conditional + failed

    if passed == n:
        color = 'Lime'
    elif conditional >= 1 and failed == 0:
        color = 'Yellow'
    else:
        color = 'Red'

    html_str = """Passed: {0}/{3}<br>Passed with conditions: {1}/{3}<br>Failed: {2}/{3} """.format(passed, conditional, failed, n)

    return html_str, color
