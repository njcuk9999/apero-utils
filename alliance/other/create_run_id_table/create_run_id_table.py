#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-12-01 at 12:20

@author: cook
"""
from astropy.table import Table, vstack
import yaml
import pandas as pd
from sqlalchemy import create_engine

# =============================================================================
# Define variables
# =============================================================================
# setup for NIRPS and SPIROU
PARAMS = dict()
# ----------------------------------------------------------------------------
# add nirps
PARAMS['NIRPS'] = dict()
# Define the CSV file containing NIRPS observation metadata
PARAMS['NIRPS']['YAML_FILE'] = 'nirps_permissions.yaml'
# define the query to run
PARAMS['NIRPS']['QUERY'] = """
SELECT KW_RUN_ID, KW_PI_NAME AS KW_PI_NAME 
FROM findex_nirps_he_online_db 
GROUP BY KW_RUN_ID  

UNION  

SELECT KW_RUN_ID, KW_PI_NAME AS KW_PI_NAME 
FROM findex_nirps_ha_online_db 
GROUP BY KW_RUN_ID;
"""
# ---- MySQL connection ----
PARAMS['NIRPS']['HOST'] = "rali.astro.umontreal.ca"
PARAMS['NIRPS']['USER'] = "nirps"
PARAMS['NIRPS']['PASSWORD'] = "Covid19!"
PARAMS['NIRPS']['DB'] = "nirps"
# define the permission groups to add
PARAMS['NIRPS']['GROUPS'] = ['ADMIN.UDEM_ASTRO', 'ADMIN.ALLIANCE']
# ----------------------------------------------------------------------------
# add spirou
PARAMS['SPIROU'] = dict()
# Define the CSV file containing NIRPS observation metadata
PARAMS['SPIROU']['YAML_FILE'] = 'spirou_permissions.yaml'
# define the query to run
PARAMS['SPIROU']['QUERY'] = """
        SELECT KW_RUN_ID, KW_PI_NAME AS KW_PI_NAME
        FROM findex_spirou_offline_db
        GROUP BY KW_RUN_ID
        """
# ---- MySQL connection ----
PARAMS['SPIROU']['HOST'] = "cosmos.astro.umontreal.ca"
PARAMS['SPIROU']['USER'] = "spirou"
PARAMS['SPIROU']['PASSWORD'] = "Covid19!"
PARAMS['SPIROU']['DB'] = "spirou"
# define the permission groups to add
PARAMS['SPIROU']['GROUPS'] = ['ADMIN.UDEM_ASTRO', 'ADMIN.ALLIANCE']


# =============================================================================
# define functions
# =============================================================================
def create_run_id_yaml(params, instrument: str):
    # get parameters from params for instrument
    yaml_filename = params[instrument]['YAML_FILE']
    query = params[instrument]['QUERY']
    host = params[instrument]['HOST']
    db = params[instrument]['DB']
    groups = params[instrument]['GROUPS']
    user = params[instrument]['USER']
    password = params[instrument]['PASSWORD']
    # ----------------------------------------------------------------------
    # SQLAlchemy connection string
    engine = create_engine(f"mysql+pymysql://{user}:{password}@{host}/{db}")
    # Execute the query and load the results into a Pandas DataFrame
    df = pd.read_sql(query, engine)
    # remove null rows
    df = df.dropna(subset=['KW_RUN_ID'])
    # Read the CSV file into an Astropy Table
    tbl0 = Table.from_pandas(df)
    # yaml dictionary
    yaml_dict = dict()
    yaml_dict['RUN_ID'] = dict()
    # loop around rows and push into dictionary for yaml writing
    for row in tbl0:
        run_id = str(row['KW_RUN_ID'])
        pi_name = str(row['KW_PI_NAME'])
        # fill out the details of the sub yaml dict
        sub_yaml_dict = dict()
        sub_yaml_dict['PI'] = str(pi_name)
        sub_yaml_dict['GROUPS'] = list(groups)
        sub_yaml_dict['USERS'] = []
        # add back to main yaml dict
        yaml_dict['RUN_ID'][run_id] = sub_yaml_dict
    # write out yaml file
    with open(yaml_filename, 'w') as yaml_file:
        yaml.dump(yaml_dict, yaml_file, default_flow_style=False)


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    for instrument in ['nirps', 'spirou']:
        create_run_id_yaml(PARAMS, instrument.upper())

# =============================================================================
# End of code
# =============================================================================