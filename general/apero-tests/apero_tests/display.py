import os
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from bokeh.io.output import output_file
from bokeh.io.saving import save
from bokeh.layouts import layout
from bokeh.models import (
    Button,
    ColumnDataSource,
    CustomJS,
    DataTable,
    HoverTool,
    LinearAxis,
    TableColumn,
)
from bokeh.models.widgets import Div, Select
from bokeh.plotting import figure


def inspect_table(
    test_html_path: Path, subtest: str, data_dict: Dict, title: str
) -> str:
    """
    Write an html table from a data set in a dictionary.
    """

    test_dir = test_html_path.parent
    test_file = test_html_path.name
    save_path = Path(test_dir, subtest, subtest + ".html")
    save_dir = save_path.parent

    save_dir.mkdir(exist_ok=True)

    source = ColumnDataSource(data_dict)

    keys_list = list(data_dict.keys())
    columns = []
    for k in keys_list:
        columns.append(TableColumn(field=k, title=k))

    parent_link = Div(
        text=f'<a href="../{test_file}">Go back</a>',
        height=25,
    )
    table_title = Div(
        text=f'<font size="+1"> <b>{title}</b> </font>',
        width=800,
        height=50,
    )
    data_table = DataTable(
        source=source,
        columns=columns,
        index_header="",
        autosize_mode="fit_columns",
        width=800,
        height=400,
        editable=True,
    )

    download = Button(label="Download to CSV", button_type="success", width=80)

    download.js_on_click(
        CustomJS(
            args=dict(source=source),
            code=open(
                os.path.join(os.path.dirname(__file__), "download.js")
            ).read(),
        )
    )

    grid_layout = layout(
        [[parent_link], [table_title], [data_table, download]]
    )

    output_file(save_path, title=subtest)
    save(grid_layout)

    # keep only subtest dir and file to put in parent html as link
    html_path = "/".join(save_path.parts[-2:])

    return html_path


def inspect_plot(test_html_path, subtest, data_dict, title):
    """
    Write an html interactive plot from a data set in a dictionary.
    """

    test_dir = test_html_path.parent
    test_file = test_html_path.name
    save_path = Path(test_dir, subtest, subtest + ".html")
    save_dir = save_path.parent

    save_dir.mkdir(exist_ok=True)

    # bokeh tools
    TOOLS = [
        "crosshair",
        "hover",
        "pan",
        "box_zoom",
        "undo",
        "redo",
        "reset",
        "save",
    ]

    if "Odometer" in data_dict:

        # night to datetime
        for i in range(len(data_dict["Night"])):
            if "_persi" in data_dict["Night"][i]:
                data_dict["Night"][i] = data_dict["Night"][i][:10]
        data_dict["PLOTDATE"] = pd.to_datetime(data_dict["Night"])

        # y variable list
        axis_map = data_dict.copy()
        axis_map.pop("Night")
        axis_map.pop("PLOTDATE")
        axis_map.pop("Odometer")

        axis_map_list = list(axis_map.keys())

        # create widget
        y_axis_widget = Select(
            title="Quality Control",
            options=axis_map_list,
            value=axis_map_list[0],
            width=260,
        )

        # data set
        data_dict["x"] = data_dict["PLOTDATE"]
        data_dict["y"] = data_dict[axis_map_list[0]]
        source_visible = ColumnDataSource(data_dict)

        # bokeh Hover
        TOOLTIPS = """
        <table>
          <tr>
            <td><span style="color: #2874a6;">Night</span></td>
            <td>@Night</td>
         </tr>
          <tr>
            <td><span style="color: #2874a6;">Odometer</span></td>
            <td>@Odometer</td>
         </tr>
         <tr>
            <td><span style="color: #2874a6;">QC Value</span></td>
            <td>@y</td>
          </tr>
        </table>
        """

        # plot
        p = figure(
            width=1200,
            height=700,
            tools=TOOLS,
            toolbar_location="left",
            x_axis_label="Night",
            x_axis_type="datetime",
            tooltips=TOOLTIPS,
            title=title,
        )
        p.title.text_font_size = "12pt"
        p.xaxis.axis_label_text_font_size = "12pt"
        p.yaxis.visible = False

    elif "Order" in data_dict:

        # y variable list
        axis_map = data_dict.copy()
        axis_map.pop("Order")

        axis_map_list = list(axis_map.keys())

        # create widget
        y_axis_widget = Select(
            title="Quality Control",
            options=axis_map_list,
            value=axis_map_list[0],
            width=260,
        )

        # data set
        data_dict["x"] = data_dict["Order"]
        data_dict["y"] = data_dict[axis_map_list[0]]
        source_visible = ColumnDataSource(data_dict)

        # bokeh Hover
        TOOLTIPS = """
        <table>
          <tr>
            <td><span style="color: #2874a6;">Order</span></td>
            <td>@x</td>
         </tr>
         <tr>
            <td><span style="color: #2874a6;">QC Value</span></td>
            <td>@y</td>
          </tr>
        </table>
        """

        # plot
        p = figure(
            width=1200,
            height=700,
            tools=TOOLS,
            toolbar_location="left",
            x_axis_label="Order",
            tooltips=TOOLTIPS,
            title=title,
        )
        p.title.text_font_size = "12pt"
        p.xaxis.axis_label_text_font_size = "12pt"
        p.yaxis.visible = False

    elif "Night" in data_dict:

        # night to datetime
        for i in range(len(data_dict["Night"])):
            if "_persi" in data_dict["Night"][i]:
                data_dict["Night"][i] = data_dict["Night"][i][:10]
        data_dict["PLOTDATE"] = pd.to_datetime(data_dict["Night"])

        # y variable list
        axis_map = data_dict.copy()
        axis_map.pop("Night")
        axis_map.pop("PLOTDATE")

        axis_map_list = list(axis_map.keys())

        # create widget
        y_axis_widget = Select(
            title="Quality Control",
            options=axis_map_list,
            value=axis_map_list[0],
            width=260,
        )

        # data set
        data_dict["x"] = data_dict["PLOTDATE"]
        data_dict["y"] = data_dict[axis_map_list[0]]
        source_visible = ColumnDataSource(data_dict)

        # bokeh Hover
        TOOLTIPS = """
        <table>
          <tr>
            <td><span style="color: #2874a6;">Night</span></td>
            <td>@Night</td>
         </tr>
         <tr>
            <td><span style="color: #2874a6;">QC Value</span></td>
            <td>@y</td>
          </tr>
        </table>
        """

        # plot
        p = figure(
            width=1200,
            height=700,
            tools=TOOLS,
            toolbar_location="left",
            x_axis_label="Night",
            x_axis_type="datetime",
            tooltips=TOOLTIPS,
            title=title,
        )
        p.title.text_font_size = "12pt"
        p.xaxis.axis_label_text_font_size = "12pt"
        p.yaxis.visible = False

    else:
        KeyError("Expected Odometer, Order or Night key for x axis")

    p.circle("x", "y", source=source_visible, line_width=2)

    y_axis = LinearAxis(
        axis_label=y_axis_widget.value, axis_label_text_font_size="12pt"
    )
    p.add_layout(y_axis, "left")

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
    callback_y_axis = CustomJS(
        args=dict(source_visible=source_visible, y_axis=y_axis, p=p),
        code=js_code,
    )

    y_axis_widget.js_on_change("value", callback_y_axis)

    # html doc
    parent_link = Div(
        text=f'<a href="../{test_file}">Go back</a>',
        height=25,
    )
    grid_layout = layout([[parent_link], [p, y_axis_widget]])

    output_file(save_path, title=subtest)
    save(grid_layout)

    # keep only subtest dir and file to put in parent html
    html_path = "/".join(save_path.parts[-2:])

    return html_path


def delta_mjd_plot(test_html_path, subtest, cdb_df, title):
    """
    Write an html interactive plot that show the time between the calibration
    file and the output file for a given recipe.

    cdb_df: Pandas MultiIndex DataFrame
    """

    test_dir = test_html_path.parent
    test_file = test_html_path.name
    save_path = Path(test_dir, subtest, subtest + ".html")
    save_dir = save_path.parent

    save_dir.mkdir(exist_ok=True)

    # list unique column names
    col_names = np.unique(cdb_df.columns.get_level_values(0))[::-1]
    # data dict to bokeh
    source = ColumnDataSource(cdb_df.reset_index(level="FILENAME"))
    # remove added underscore
    source.data["FILENAME"] = source.data.pop("FILENAME_")
    # night to datetime
    source.data["PLOTDATE"] = pd.to_datetime(source.data["OBS_DIR"])

    htool = HoverTool(
        tooltips=[
            ("OBS_DIR", "@OBS_DIR"),
            ("FILENAME", "@FILENAME"),
            ("Calib file", "@calf"),
        ]
    )

    # bokeh tools
    TOOLS = [
        "crosshair",
        htool,
        "pan",
        "box_zoom",
        "undo",
        "redo",
        "reset",
        "save",
    ]

    # create widget
    y_axis_widget = Select(
        title="CDBTYPE", options=list(col_names), value=col_names[0], width=260
    )

    # data set (x and y variables)
    source.data["x"] = source.data["PLOTDATE"]
    source.data["y"] = source.data[col_names[0] + "_" + "DELTA_MJD"]
    source.data["calf"] = source.data[col_names[0] + "_" + "CALIB_FILE"]

    # plot
    p = figure(
        width=1200,
        height=700,
        tools=TOOLS,
        toolbar_location="right",
        x_axis_label="OBS_DIR",
        x_axis_type="datetime",
        title=title,
    )
    p.title.text_font_size = "12pt"
    p.xaxis.axis_label_text_font_size = "12pt"
    p.yaxis.visible = False

    p.circle("x", "y", source=source, line_width=2)

    y_axis = LinearAxis(
        axis_label="DELTA_MJD", axis_label_text_font_size="12pt"
    )
    p.add_layout(y_axis, "left")

    # javascript callback
    js_code = """
        var selected_y_axis = cb_obj.value
        var data_visible = source.data
        data_visible.y = data_visible[selected_y_axis + '_' + 'DELTA_MJD']
        data_visible.calf = data_visible[selected_y_axis + '_' + 'CALIB_FILE']
        source.change.emit()
        p.reset.emit()
        """

    callback_y_axis = CustomJS(args=dict(source=source, p=p), code=js_code)

    y_axis_widget.js_on_change("value", callback_y_axis)

    # html doc
    parent_link = Div(
        text=f'<a href="../{test_file}">Go back</a>',
        height=25,
    )
    grid_layout = layout([[parent_link], [p, y_axis_widget]])

    output_file(save_path, title=subtest)
    save(grid_layout)

    # keep only subtest dir and file to put in parent html
    html_path = "/".join(save_path.parts[-2:])

    return html_path
