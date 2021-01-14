"""
Utility functions for APERO tests

@author: charles
"""
import os
from pathlib import Path
import glob
import pandas as pd
from datetime import datetime
from astropy.io import fits
from bokeh.plotting import figure
from bokeh.io.saving import save
from bokeh.io.output import output_file
from bokeh.models import ColumnDataSource, DataTable, DateFormatter, TableColumn
from bokeh.models import LinearAxis, BasicTicker, CustomJS
from bokeh.models.widgets import Div, Select
from bokeh.layouts import row, column


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
    # TODO: Make this safer and keep only strings with a date-compatible fmt
    # or use a reject_dir list to allow various date formats, maybe this can
    # be checked by datetime
    if 'other' in nlist:
        nlist.remove('other')
    nlist.sort(key=lambda date: datetime.strptime(date[:10], "%Y-%m-%d"))

    return nlist


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
        print('No ', f, 'in', path)
    return flist


def count_files(path, files='all'):
    """
    Count number of files of a give pattern in a directory.
    """

    if files == 'all':
        fpattern = '*'
    else:
        fpattern = files

    return len(glob.glob(os.path.join(path, fpattern)))


def count_files_subdir(path, subdir='all', files='all'):
    """
    Count files in selected subdirectories of a given dir.
    """

    if subdir == 'all':
        spattern = '*'
    else:
        spattern = subdir

    if files == 'all':
        fpattern = '*'
    else:
        fpattern = files

    return len(glob.glob(os.path.join(path, spattern, fpattern)))


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
        self.PID = tbl['PID']
        self.QC = tbl['PASSED_ALL_QC']
        self.ENDED = tbl['ENDED']

        self.indexQCtrue = tbl['PASSED_ALL_QC']
        self.indexQCfalse = ~tbl['PASSED_ALL_QC']
        self.indexENDEDtrue = tbl['ENDED']
        self.indexENDEDfalse = ~tbl['ENDED']

        self.nights = tbl['DIRECTORY']
        self.runstr = tbl['RUNSTRING']
        self.args = tbl['ARGS']

        self.QCstr = tbl['QC_STRING']
        self.QCnames = tbl['QC_NAMES']
        self.QCvalues = tbl['QC_VALUES']
        self.QClogic = tbl['QC_LOGIC']
        self.ERRORS = tbl['ERRORS']
        self.LOGFILE = tbl['LOGFILE']


def inspect_table(test, subtest, data_dict, title):
    """
    Write an html table from a data set in a dictionary.
    """

    p = Path('..', 'out', test, subtest)
    if not p.is_dir():
        p.mkdir(exist_ok=True)

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
            max_width=1600,
            width_policy='min',
            height=600,
            editable=True)
    layout = column(table_title, data_table)

    output_file(
            os.path.join('..', 'out', test, subtest, subtest+'.html'),
            title=subtest)
    save(layout)

    html_str = """<a href='{0}/{0}.html'>Inspect</a>""".format(subtest)

    return html_str


def inspect_plot(test, subtest, data_dict, title, order = False):
    """
    Write an html interactive plot from a data set in a dictionary.
    """

    p = Path('..', 'out', test, subtest)
    if not p.is_dir():
        p.mkdir(exist_ok=True)

    # bokeh tools
    TOOLS = ["crosshair", "hover", "pan", "box_zoom", "undo", "redo", "reset",
             "save", "tap"]

    if order == True:

        # y variable list
        axis_map = data_dict.copy()
        axis_map.pop('Order')

        axis_map_list = list(axis_map.keys())

        # create widget
        y_axis_widget = Select(title="Quality Control",
                        options=axis_map_list,
                        value = axis_map_list[0],
                        width=260)

        # data set
        data_dict['x'] = data_dict['Order']
        data_dict['y'] = data_dict[axis_map_list[0]]
        source_visible = ColumnDataSource(data_dict)

        # plot
        p = figure(plot_width=1200,
                   plot_height=700,
                   tools=TOOLS,
                   toolbar_location = "left",
                   x_axis_label = 'Order',
                   title = title)
        p.title.text_font_size = '12pt'
        p.xaxis.axis_label_text_font_size = '12pt'
        p.yaxis.visible = False

    else:

        # Night to datetime
        for i in range(len(data_dict['Night'])):
            if '_persi' in data_dict['Night'][i]:
                data_dict['Night'][i] = data_dict['Night'][i][:10]
        data_dict['Night'] = pd.to_datetime(data_dict['Night'])

        # y variable list
        axis_map = data_dict.copy()
        axis_map.pop('Night')

        axis_map_list = list(axis_map.keys())

        # create widget
        y_axis_widget = Select(title="Quality Control",
                        options=axis_map_list,
                        value = axis_map_list[0],
                        width=260)

        # data set
        data_dict['x'] = data_dict['Night']
        data_dict['y'] = data_dict[axis_map_list[0]]
        source_visible = ColumnDataSource(data_dict)

        # plot
        p = figure(plot_width=1200,
                   plot_height=700,
                   tools=TOOLS,
                   toolbar_location = "left",
                   x_axis_label = 'Night',
                   x_axis_type="datetime",
                   title = title)
        p.title.text_font_size = '12pt'
        p.xaxis.axis_label_text_font_size = '12pt'
        p.yaxis.visible = False


    p.circle('x',
           'y',
           source = source_visible,
           line_width=2
           )

    y_axis = LinearAxis(axis_label = y_axis_widget.value,
                  axis_label_text_font_size = '12pt')
    p.add_layout(y_axis, 'left')

    # javascript callback
    callback_y_axis = CustomJS(args=dict(source_visible=source_visible,
                                         y_axis=y_axis,
                                         p=p),
                      code="""
                      var selected_y_axis = cb_obj.value
                      var data_visible = source_visible.data

                      data_visible.y = data_visible[selected_y_axis]
                      source_visible.change.emit()
                      y_axis.axis_label = selected_y_axis
                      y_axis.change.emit()

                      p.reset.emit()
                           """)

    y_axis_widget.js_on_change('value', callback_y_axis)

    #html doc
    layout = row(y_axis_widget, p)

    output_file(
            os.path.join('..', 'out', test, subtest, subtest+'.html'),
            title=subtest)
    save(layout)

    html_str = """<a href='{0}/{0}.html'>Inspect</a>""".format(subtest)

    return html_str


def summary(test):
    """
    Write the html summary.
    """

    f = open("../out/{0}/{0}.html".format(test), "r")
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
