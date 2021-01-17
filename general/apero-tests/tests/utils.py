"""
Utility functions for APERO tests

@author: charles
"""
import os
from pathlib import Path
from datetime import datetime

import glob
import pandas as pd
from bokeh.plotting import figure
from bokeh.io.saving import save
from bokeh.io.output import output_file
from bokeh.models import ColumnDataSource, DataTable, TableColumn
from bokeh.models import LinearAxis, CustomJS
from bokeh.models.widgets import Div, Select
from bokeh.layouts import column, layout

from . import OUTDIR


def list_nights(path):
    """
    Return a list of datetime objects corresponding to each night
    """
    nlist = [x
             for x in os.listdir(path)
             if os.path.isdir(os.path.join(path, x))
             ]
    # NOTE: may need to change this if statement in v0.7
    if 'other' in nlist:
        nlist.remove('other')
    nlist.sort(key=lambda date: datetime.strptime(date[:10], "%Y-%m-%d"))

    return nlist


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


def inspect_table(test, subtest, data_dict, title):
    """
    Write an html table from a data set in a dictionary.
    """

    save_path = Path(OUTDIR, test, subtest)
    save_path.mkdir(exist_ok=True)

    source = ColumnDataSource(data_dict)

    keys_list = list(data_dict.keys())
    columns = []
    for k in keys_list:
        columns.append(TableColumn(field=k, title=k))

    parent_link = Div(
            text=f'<a href="../{test}.html">Go back</a>',
            height=25,
            )
    table_title = Div(
            text=f'<font size="+1"> <b>{title}</b> </font>',
            height=50,
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
    grid_layout = column(parent_link, table_title, data_table)

    output_file(os.path.join(str(save_path), subtest+'.html'), title=subtest)
    save(grid_layout)

    # Keep only subtest dir and file to put in parent html
    html_path = '/'.join([save_path.parts[-1], subtest+'.html'])

    return html_path


def inspect_plot(test, subtest, data_dict, title):
    """
    Write an html interactive plot from a data set in a dictionary.
    """

    save_path = Path(OUTDIR, test, subtest)
    save_path.mkdir(exist_ok=True)

    # bokeh tools
    TOOLS = ["crosshair", "hover", "pan", "box_zoom", "undo", "redo", "reset",
             "save", "tap"]

    if 'Order' in data_dict:

        # y variable list
        axis_map = data_dict.copy()
        axis_map.pop('Order')

        axis_map_list = list(axis_map.keys())

        # create widget
        y_axis_widget = Select(title="Quality Control",
                        options=axis_map_list,
                        value=axis_map_list[0],
                        width=260)

        # data set
        data_dict['x'] = data_dict['Order']
        data_dict['y'] = data_dict[axis_map_list[0]]
        source_visible = ColumnDataSource(data_dict)

        # plot
        p = figure(plot_width=1200,
                   plot_height=700,
                   tools=TOOLS,
                   toolbar_location="left",
                   x_axis_label='Order',
                   title=title)
        p.title.text_font_size = '12pt'
        p.xaxis.axis_label_text_font_size = '12pt'
        p.yaxis.visible = False

    elif 'Night' in data_dict:

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
    else:
        KeyError('Expected Order or Night key for x axis')

    p.circle('x', 'y', source=source_visible, line_width=2)

    y_axis = LinearAxis(axis_label=y_axis_widget.value,
                        axis_label_text_font_size='12pt')
    p.add_layout(y_axis, 'left')

    # javascript callback
    js_code = """
              var selected_y_axis = cb_obj.value
              var data_visible = source_visible.data

              data_visible.y = data_visible[selected_y_axis]
              source_visible.change.emit()
              y_axis.axis_label = selected_y_axis
              y_axis.change.emit()

              p.reset.emit()
              """
    callback_y_axis = CustomJS(args=dict(source_visible=source_visible,
                                         y_axis=y_axis,
                                         p=p),
                               code=js_code)

    y_axis_widget.js_on_change('value', callback_y_axis)

    #html doc
    parent_link = Div(
            text=f'<a href="../{test}.html">Go back</a>',
            height=25,
            )
    grid_layout = layout([
        [parent_link],
        [y_axis_widget, p]
        ])

    output_file(os.path.join(str(save_path), subtest+'.html'), title=subtest)
    save(grid_layout)

    # Keep only subtest dir and file to put in parent html
    html_path = '/'.join([save_path.parts[-1], subtest+'.html'])

    return html_path


def summary(test):
    """
    Write the html summary.
    """

    f = open(os.path.join(OUTDIR, test, test+'.html'), 'r')
    html = f.read()

    npassed = html.count('Lime')
    ncond = html.count('Yellow')
    nfailed = html.count('Red')

    ntotal = npassed + ncond + nfailed

    if npassed == ntotal:
        color = 'Lime'
    elif ncond >= 1 and nfailed == 0:
        color = 'Yellow'
    else:
        color = 'Red'

    return ntotal, npassed, ncond, nfailed, color
