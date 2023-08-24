#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on {DATE}

@author: cook
"""
import json
from typing import Dict, List, Union

import numpy as np
from tqdm import tqdm
from astropy.time import Time

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = ''


# =============================================================================
# Define functions
# =============================================================================
def python_str_to_html_str(value):
    value = str(value)
    value = value.replace('\n', '<br>')
    value = value.replace('\t', '&nbsp;&nbsp;&nbsp;&nbsp;')
    return value


def filtered_html_table(outlist: Dict[int, Dict[str, Union[str, List[str]]]],
                        col_names: List[str],
                        col_types: List[str]):
    """
    Generate a html page with a table of data that can be filtered by column
    values.

    :param outlist: dictionary of dictionaries containing the data to display
                    each primary key is a number (must contain at least 1 row)
    :param col_names:
    :param col_types:
    :return:
    """
    # -------------------------------------------------------------------------
    print('Cleaning table data')
    # clean up outlist
    for idnumber in tqdm(outlist):
        for c_it, column_name in enumerate(col_names):
            if col_types[c_it] == 'list':
                for r_it, row in enumerate(outlist[idnumber][column_name]):
                    value = python_str_to_html_str(row)
                    outlist[idnumber][column_name][r_it] = value
            else:
                value = python_str_to_html_str(outlist[idnumber][column_name])
                outlist[idnumber][column_name] = value
    # -------------------------------------------------------------------------
    print('Generating html page')
    # get the column headers text in html format
    column_headers = "\n".join([f'<th>{col}</th>' for col in col_names])

    # some require special formatting
    filter_data_cols = []
    render_data_cols = []

    for c_it, column_name in enumerate(col_names):
        if col_types[c_it] == 'list':
            filter_data_cols.append(f'row.{column_name}.join("")'
                                    f'.toLowerCase().includes(filterValue)')
            render_data_cols.append(f'<td>${{row.{column_name}'
                                    f'.join("<br>")}}</td>')
        else:
            filter_data_cols.append(f'\nrow.{column_name}'
                                    f'.toLowerCase().includes(filterValue)')
            render_data_cols.append(f'<td>${{row.{column_name}}}</td>')
    # push into string format
    filter_col_str = ' ||'.join(filter_data_cols)
    render_col_str = '\n'.join(render_data_cols)


    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            table {{
                font-family: Arial, sans-serif;
                border-collapse: collapse;
                width: 100%;
            }}
            th, td {{
                border: 1px solid #dddddd;
                text-align: left;
                padding: 8px;
            }}
            th {{
                background-color: #f2f2f2;
            }}
            #loading {{
                text-align: center;
            }}
        </style>
    </head>
    <body>
        <h1>Job Status Report</h1>
        
        <div>
            <label for="filterSelect">Filter by column:</label>
            <select id="filterSelect"></select>
            <input type="text" id="filterInput" 
                   placeholder="Enter filter value">
            <button id="applyFilter">Apply Filter</button>
        </div>
        
        <br><br>
        
        <table id="jobTable">
            <tr>
                """ + column_headers + """
            </tr>
        </table>
        <div id="loading">Loading...</div>
        <script>
            const outlist = JSON.parse(`""" + json.dumps(outlist) + """`);
            
            const filterSelect = document.getElementById("filterSelect");
            const filterInput = document.getElementById("filterInput");
            const applyFilter = document.getElementById("applyFilter");
            const table = document.getElementById("jobTable");
            const loadingDiv = document.getElementById("loading");
            
            const columns = Object.keys(outlist[1]);
            columns.forEach(column => {
                const option = document.createElement("option");
                option.text = column;
                filterSelect.add(option);
            });
            
            applyFilter.addEventListener("click", loadFilteredRows);
            
            function loadFilteredRows() {
                const filterValue = filterInput.value.toLowerCase();
                const filteredData = Object.values(outlist).filter(row => 
                    """ + filter_col_str + """
                );
                
                renderRows(filteredData);
            }
            
            function renderRows(data) {
                table.innerHTML = `
                    <tr>
                        """ + column_headers + """
                    </tr>
                `;
                
                data.forEach(row => {
                    const newRow = table.insertRow();
                    newRow.innerHTML = `
                        """ + render_col_str + """
                    `;
                });
                
                loadingDiv.style.display = "none";
            }
            
            loadFilteredRows();
        </script>
    </body>
    </html>
    """

    return html_content


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    # Define the start and end times
    start_time = Time("2023-08-01")
    end_time = Time("2023-08-10")

    outlist = dict()
    column_names = ['Timestamp', 'RECIPE', 'PID', 'HAS_ERROR', 'ERROR']
    column_types = ['string', 'string', 'string', 'string', 'list']

    for i in range(10000):
        if i % 2 == 0:
            has_error = True
            error = ['error 1: This is the text for error 1',
                     'error 2: This is the text for error 2']
        else:
            has_error = False
            error = []

        time = start_time + (end_time - start_time) * np.random.rand()

        values = dict(RECIPE="recipe name",
                      PID=f"pid-{i}", HAS_ERROR=str(has_error),
                      ERROR=error,
                      Timestamp=time.iso)
        # outlist dictionary must be in same order as column names
        outlist[i] = dict()
        for col_name in column_names:
            outlist[i][col_name] = values[col_name]

    html_page = filtered_html_table(outlist, col_names=column_names,
                                    col_types=column_types)

    with open("job_report.html", "w") as f:
        f.write(html_page)

    print("HTML page with filtering generated successfully.")


# =============================================================================
# End of code
# =============================================================================
