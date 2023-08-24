#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
We create a directory of PID files for each night

Each PID file contains a yaml file with the following information:

RECIPE: str
SHORTNAME: str
RUNSTRING: str
GROUP: str
GROUP DATE: str
PID: str
TIME: str
ERRORS: list
WARNINGS: list

Created on 2023-08-24

@author: cook
"""
import glob
import os
from typing import Any, Dict

import yaml
from astropy.time import Time
from tqdm import tqdm


# =============================================================================
# Define variables
# =============================================================================


# =============================================================================
# Define functions
# =============================================================================
def apero_group_to_date(apero_group):
    try:
        # get the unix string
        raw_unix_str = apero_group.split('-')[2]
        # convert to float
        raw_unix_float = float(raw_unix_str)
        # convert to a time
        unix_time = Time(raw_unix_float / 1e7, format='unix')
        # get obs-dir in our standard format
        obs_dir = unix_time.strftime('%Y-%m-%d')
        # return obs-dir
        return obs_dir
    except Exception as _:
        raise ValueError(f'apero group = {apero_group} invalid')


def pid_to_time(pid):
    try:
        # get the unix string
        raw_unix_str = pid.split('-')[1]
        # convert to float
        raw_unix_float = float(raw_unix_str)
        # convert to a time
        unix_time = Time(raw_unix_float / 1e7, format='unix')
        # get obs-dir in our standard format
        time = unix_time.iso
        # return obs-dir
        return time
    except Exception as _:
        raise ValueError(f'pid = {pid} invalid')


def dict_to_yaml(yaml_dict: Dict[str, Any], yaml_file: str):
    """
    Save a dictionary to a yaml file
    """
    # open yaml file
    with open(yaml_file, 'w') as wfile:
        # dump yaml dict to yaml file
        yaml.dump(yaml_dict, wfile)


def yaml_to_dict(yaml_file: str):
    """
    Load a dictionary from a yaml file
    """
    # open yaml file
    with open(yaml_file, 'r') as rfile:
        # load yaml file
        yaml_dict = yaml.load(rfile, Loader=yaml.FullLoader)
    # return yaml dict
    return yaml_dict


def from_outlist(save_path: str, outlist: dict):
    for item in tqdm(outlist):

        # create dictionary for saving to yaml
        yaml_dict = dict()

        # get recipe name
        yaml_dict['RECIPE'] = item['RECIPE']
        # ---------------------------------------------------------------------
        # get shortname
        shortname = ''
        if 'SHORTNAME' in item:
            shortname = item['SHORTNAME']
        elif 'ARGS' in item:
            if 'shortname' in item['ARGS']:
                shortname = item['ARGS']['shortname']
        yaml_dict['SHORTNAME'] = shortname
        # ---------------------------------------------------------------------
        # get runstring
        yaml_dict['RUNSTRING'] = item['RUNSTRING']
        # ---------------------------------------------------------------------
        # get pid
        yaml_dict['PID'] = item['PID']
        # ---------------------------------------------------------------------
        # get time
        yaml_dict['TIME'] = pid_to_time(item['PID'])
        # ---------------------------------------------------------------------
        # get group
        apero_group = ''
        if 'ARGS' in item:
            if 'DRS_GROUP' in item['ARGS']:
                apero_group = item['ARGS']['DRS_GROUP']
        yaml_dict['GROUP'] = apero_group
        # ---------------------------------------------------------------------
        # get group date
        try:
            group_date = apero_group_to_date(apero_group)
            yaml_dict['GROUP_DATE'] = group_date
        except ValueError:
            # get from the PID
            group_date = Time(yaml_dict['TIME']).strftime('%Y-%m-%d')
            yaml_dict['GROUP_DATE'] = group_date
        # ---------------------------------------------------------------------
        # get errors
        errors = []
        if 'ERRORS' in item:
            errors = item['ERRORS']
        yaml_dict['ERRORS'] = errors
        # ---------------------------------------------------------------------
        # get warnings
        warnings = []
        if 'WARNINGS' in item:
            warnings = item['WARNINGS']
        yaml_dict['WARNINGS'] = warnings
        # ---------------------------------------------------------------------
        # save to yaml dict
        to_yaml_dict(save_path, item['PID'], group_date, yaml_dict)


def to_yaml_dict(save_path, save_name, group_date, yaml_dict):
    # ---------------------------------------------------------------------
    # construct path to yaml file
    yaml_file = '{0}.yaml'.format(save_name)
    yaml_path = os.path.join(save_path, group_date)
    # make directory if it does not exist
    if not os.path.exists(yaml_path):
        os.makedirs(yaml_path)
    # save to yaml file
    dict_to_yaml(yaml_dict, os.path.join(yaml_path, yaml_file))


def from_processing_log_file(save_path: str, log_filename: str):
    # load file
    with open(log_filename, 'r') as rfile:
        lines = rfile.readlines()

    # find lines containing "Error found for ID="
    error_lines_start = []
    for row, line in enumerate(lines):
        if 'Error found for ID=' in line:
            error_lines_start.append(row)

    # find the end of the error line
    error_lines_end = error_lines_start[1:] + [len(lines) - 1]

    # now we to make a yaml for each of these
    for it, row in tqdm(enumerate(error_lines_start)):
        # get the recipe
        runstring = lines[row + 1].split('@!|PROC|')[-1]
        # the recipe is the first part of the runstring
        recipe = runstring.split(' ')[0]
        # the short name may be in the runstring after --shortname
        shortname = ''
        if '--shortname=' in runstring:
            shortname = runstring.split('--shortname=')[-1].split(' ')[0]
        elif '--shortname ' in runstring:
            shortname = runstring.split('--shortname ')[-1].split(' ')[0]
        # get the apero log pid
        proc_id = os.path.basename(log_filename).split('.log')[0]
        # the pid cannot be found we have to create one from the filename
        #   and the ID
        idnumber = int(lines[row].split('ID=')[-1].split('\'')[1])
        pid = f'{proc_id}-{idnumber:07d}'
        # get the group from the log_filename
        group = proc_id
        group_date = apero_group_to_date(group)
        # get the time from the group
        time = group_date
        # get the errors
        errors = []
        # errors start from row + 4 and end when we reach another |PROC|
        for row2 in range(row + 4, error_lines_end[it]):
            # get line
            line = lines[row2]
            # if we reach another |PROC| then we are done
            if '|PROC|' in line:
                break
            # append to errors
            errors.append(line)
        # we cannot get the warnings
        warnings = []
        # ---------------------------------------------------------------------
        # create dictionary for saving to yaml
        yaml_dict = dict()
        # get recipe name
        yaml_dict['RECIPE'] = recipe
        # ---------------------------------------------------------------------
        # get shortname
        yaml_dict['SHORTNAME'] = shortname
        # ---------------------------------------------------------------------
        # get runstring
        yaml_dict['RUNSTRING'] = runstring
        # ---------------------------------------------------------------------
        # get pid
        yaml_dict['PID'] = pid
        # ---------------------------------------------------------------------
        # get time
        yaml_dict['TIME'] = time
        # ---------------------------------------------------------------------
        # get group
        yaml_dict['GROUP'] = group
        # ---------------------------------------------------------------------
        # get group date
        yaml_dict['GROUP_DATE'] = group_date
        # ---------------------------------------------------------------------
        # get errors
        yaml_dict['ERRORS'] = errors
        # ---------------------------------------------------------------------
        # get warnings
        yaml_dict['WARNINGS'] = warnings
        # ---------------------------------------------------------------------
        # save to yaml dict
        to_yaml_dict(save_path, pid, group_date, yaml_dict)





# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    _logpath = ('/scratch2/nirps/misc/apero_processing_errors/')

    _save_path = '/scratch2/nirps/misc/apero_processing_errors/yamls/'

    yaml_files = glob.glob(os.path.join(_logpath, '*.log'))

    for _yaml_file in yaml_files:
        from_processing_log_file(_save_path, _yaml_file)

# =============================================================================
# End of code
# =============================================================================
