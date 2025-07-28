#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 2025-07-28

@author: cook
"""
import os
from typing import List, Dict, Tuple
from astropy.table import Table
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
# the PID of the APERO processing run
# Get the PID from here: /cosmos99/spirou/apero-data/spirou_offline/msg/tool/other
APERO_PID = 'PID-00017534957092926630-VLWQ'
# the working directory where the log and report files are located
WORKING_DIR = '/cosmos99/spirou/apero-data/spirou_offline/'
# the paths to the log and report files
PATH_TO_LOG = os.path.join(WORKING_DIR, 'msg/tool/other')
LOG_FILE = 'APEROL-{APERO_PID}_apero_processing.log'
# the path to the report file
PATH_TO_REPORT = os.path.join(WORKING_DIR, 'msg/report/processing')
REPORT_FILE = '{APERO_PID}_apero_processing_ids.txt'


# =============================================================================
# Define functions
# =============================================================================
def get_log_file_errors(lines: List[str]
                        )-> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Function to extract errors from the log file.
    """
    errors = dict()
    error_messages = dict()
    # print progress
    print("Extracting errors from log file...")

    # loop through the lines of the log file
    # look for "W[40-503-00019]: Error found for ID='{id}'"
    # extract the ID and take the next line as the recipe 
    # (split by "apero" and ".py")
    for i, line in tqdm(enumerate(lines)):
        if "W[40-503-00019]: Error found for ID='" in line:
            # Extract the ID
            id_start = line.index("ID='") + 4
            id_end = line.index("'", id_start)
            apero_id = line[id_start:id_end]

            # Get the next line for the recipe
            if i + 1 < len(lines):
                recipe_line = lines[i + 1].strip()
                if 'apero' in recipe_line and '.py' in recipe_line:
                    recipe = recipe_line.split('apero')[1].split('.py')[0].strip()
                    errors[apero_id] = f'apero{recipe}.py'
            # Get the error message (Stored on the next lines between two 
            #  sets of ***** lines)
            error_message = []
            for j in range(i + 2, len(lines)):
                if '*' * 10 in lines[j]:
                    if error_message:
                        # If we hit another set of ***** lines, stop collecting
                        break
                elif len(lines[j].strip()) == 0:
                    # If we hit an empty line, continue to the next line
                    continue
                else:
                    # only keep characters after the PROC message
                    if '|PROC|' in lines[j]:
                        error_message.append(lines[j].split('|PROC|')[-1].strip())
                    else:
                        # Otherwise, just append the line   
                        error_message.append(lines[j].strip())
            if error_message:
                error_messages[apero_id] = ' '.join(error_message).strip()
    # Now we have the errors and their messages, we can return them
    return errors, error_messages


def count_by_recipe(mydict: Dict[str, str]) -> Dict[str, int]:
    """
    Function to count the number of errors by recipe.
    """
    recipe_count = {}
    for recipe in mydict.values():
        if recipe in recipe_count:
            recipe_count[recipe] += 1
        else:
            recipe_count[recipe] = 1
    return recipe_count


def get_report_file_runs(lines: List[str]) -> Dict[str, str]:
    """
    From the file find all lines that start with "id" 
    extract the ID (the number following preceeding an "=" and then extract 
    the recipe name (follows the "=" and starts with "apero" and ends with ".py"))
    """
    runs = dict()
    # print progress
    print("Extracting runs from report file...")
    # loop through the lines of the report file
    for line in tqdm(lines):
        if line.startswith("id"):
            # Extract the ID
            id_start = line.index("id") + 3
            id_end = line.index(" ", id_start)
            apero_id = line[id_start:id_end]

            # Extract the recipe name
            recipe_start = line.index("apero")
            recipe_end = line.index(".py", recipe_start) + 3
            recipe = line[recipe_start:recipe_end].strip()

            runs[apero_id] = recipe

    return runs


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # load the log file)
    log_file_path = os.path.join(PATH_TO_LOG, LOG_FILE.format(APERO_PID=APERO_PID))
    # deal with log file not existing
    if not os.path.exists(log_file_path):
        print(f"Log file does not exist: {log_file_path}")
        raise ValueError(f"Log file does not exist: {log_file_path}")

    # Read the log file
    with open(log_file_path, 'r') as log_file:
        log_lines = log_file.readlines()

    # get the errors from the log file
    errors, error_msgs = get_log_file_errors(log_lines)

    # count how many errors of each recipe there were
    error_counts = count_by_recipe(errors)

    # ----------------------------------------------------------------------
    # load the report file
    report_file_path = os.path.join(PATH_TO_REPORT, REPORT_FILE.format(APERO_PID=APERO_PID))
    # deal with report file not existing
    if not os.path.exists(report_file_path):
        print(f"Report file does not exist: {report_file_path}")
        raise ValueError(f"Report file does not exist: {report_file_path}")
    
    # Read the report file
    with open(report_file_path, 'r') as report_file:
        report_lines = report_file.readlines()
    
    # get the runs from the report file
    runs = get_report_file_runs(report_lines)

    # count how many runs there were
    run_counts = count_by_recipe(runs)


    # -------------------------------------------------------------------------
    # get the unique recipe names from the errors and runs
    unique_recipes = set(error_counts.keys()).union(set(run_counts.keys()))

    # print how many errors and runs there were for each recipe
    print("Errors by recipe:")
    for recipe in unique_recipes:
        error_count = error_counts.get(recipe, 0)
        run_count = run_counts.get(recipe, 'inf')
        print(f"{recipe}: {error_count} errors, {run_count} runs")

    # -------------------------------------------------------------------------
    # save the error messages to a fits Bin Table file 
    # columns = ID, recipe, error message
    error_file_path = os.path.join(PATH_TO_LOG, f'error_table_{APERO_PID}.fits')
    
    error_table = Table(names=('ID', 'recipe', 'error_message'),
                        dtype=('S20', 'S50', 'S200'))
    for apero_id, recipe in errors.items():
        error_message = error_msgs.get(apero_id, '')
        error_table.add_row((apero_id.encode('utf-8'), recipe.encode('utf-8'), 
                             error_message.encode('utf-8')))
    error_table.write(error_file_path, format='fits', overwrite=True)







# =============================================================================
# End of code
# =============================================================================
