#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on {DATE}

@author: cook
"""
import json
from astropy.time import Time
import numpy as np

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = ''

# =============================================================================
# Define functions
# =============================================================================
def generate_lazy_loading_html(outlist):
    html_content = """
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            table {
                font-family: Arial, sans-serif;
                border-collapse: collapse;
                width: 100%;
            }
            th, td {
                border: 1px solid #dddddd;
                text-align: left;
                padding: 8px;
            }
            th {
                background-color: #f2f2f2;
            }
            #loading {
                text-align: center;
            }
        </style>
    </head>
    <body>
        <h1>Job Status Report</h1>
        
        <div>
            <label for="filterSelect">Filter by column:</label>
            <select id="filterSelect"></select>
            <input type="text" id="filterInput" placeholder="Enter filter value">
            <button id="applyFilter">Apply Filter</button>
        </div>
        
        <table id="jobTable">
            <tr>
                <th>Timestamp</th>
                <th>Recipe Name</th>
                <th>PID</th>
                <th>Error</th>
                <th>Error Log</th>

            </tr>
        </table>
        
        <div id="loading">Loading...</div>
        
        <script>
            let currentPage = 0;
            const rowsPerPage = 100;
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

            applyFilter.addEventListener("click", filterTable);

            function filterTable() {
                const selectedColumn = filterSelect.value;
                const filterValue = filterInput.value.toLowerCase();
                const rows = table.getElementsByTagName("tr");
                

                
                for (let i = 1; i < rows.length; i++) {
                    const row = rows[i];
                    const cell = row.getElementsByTagName("td")[columns.indexOf(selectedColumn)];
                    
                    console.log('i=', i, cell.textContent || cell.innerText, filterValue);
                    
                    if (cell) {    
                        const cellValue = cell.textContent || cell.innerText;
                        row.style.display = cellValue.toLowerCase().includes(filterValue) ? "" : "none";
                    }                  
                }
            }
            
            function loadMoreRows() {
                const data = Object.values(outlist).slice(currentPage * rowsPerPage, (currentPage + 1) * rowsPerPage);

                data.forEach(row => {
                    const newRow = table.insertRow();
                    newRow.innerHTML = `
                    
                        <td>${row.Timestamp}</td>
                        <td>${row.RECIPE}</td>
                        <td>${row.PID}</td>
                        <td>${row.ERROR.length > 0 ? "Yes" : "No"}</td>
                        <td>${row.ERROR.join("<br>")}</td>


                    `;
                });

                currentPage++;
                loadingDiv.style.display = "none";
            }

            loadMoreRows();

            window.addEventListener("scroll", () => {
                const { scrollTop, clientHeight, scrollHeight } = document.documentElement;
                if (scrollTop + clientHeight >= scrollHeight - 50) {
                    loadMoreRows();
                }
            });
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

    for i in range(10000):
        if i % 2 == 0:
            error = [['error 1'], ['error 2']]
        else:
            error = []

        time = start_time + (end_time - start_time) * np.random.rand()

        outlist[i] = dict(Timestamp=time.iso,
                          RECIPE="recipe name",
                          PID=f"pid-{i}", ERROR=error)

    html_page = generate_lazy_loading_html(outlist)

    with open("lazy_loading_job_report.html", "w") as f:
        f.write(html_page)

    print("Lazy loading HTML page with filtering generated successfully.")

# =============================================================================
# End of code
# =============================================================================
