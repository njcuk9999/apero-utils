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
    # TODO: Make this safer and keep only strings with a date-compatible fmt
    # or use a reject_dir list to allow various darte formats, maybe this can
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
        print(path)
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


def inspect_table(test, check, data_dict, title):
    """
    Write an 'inspect' html element.
    """

    # TODO: handle path better and use os.path.join
    if not os.path.isdir('../out/{0}/{1}'.format(test, check)):
        os.system('mkdir ../out/{0}/{1}'.format(test, check))

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
            "../out/{0}/{1}/{1}.html".format(test, check),
            title="{0}".format(check))
    save(layout)

    html_str = """<a href='{0}/{0}.html'>Inspect</a>""".format(check)

    return html_str


def inspect_plot(test, check, data_dict, title):
    """
    WIP
    """

    if not os.path.isdir('../out/{0}/{1}'.format(test, check)):
        os.system('mkdir ./out/{0}/{1}'.format(test, check))

    # y variable list
    axis_map = data_dict
    axis_map_list = list(axis_map.keys())

    # create widget
    y_axis = Select(title="Quality Control",
                    options=axis_map_list,
                    width=260)

    # data sets
    source = ColumnDataSource(data_dict)

    # bokeh tools
    TOOLS = ["crosshair", "hover", "pan", "box_zoom", "undo", "redo", "reset",
             "save", "tap"]

    # titles
    titlestr = title

    p = figure(plot_width=1200, plot_height=700, tools=TOOLS,
               tooltips=TOOLTIPS, title = titlestr)
    p.title.text_font_size = '12pt'

    plot = p.circle('Night',
                    axis_map_list[0],
                    source=source
                    )

    table_title = Div(text="""<font size="+1"> <b>{0}</b> </font>""".format(title), height=50)
    data_table = DataTable(source=source, columns=columns,
                           autosize_mode='fit_columns', min_width=1000,
                           max_width=1500, width_policy='min', height=600,
                           editable=True)
    layout = column(table_title, data_table)

    output_file(
            "../out/{0}/{1}/{1}.html".format(test, check), title="{0}".format(check)
            )
    save(layout)

    html_str = """<a href='{0}/{0}.html'>Inspect</a>""".format(check)

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
