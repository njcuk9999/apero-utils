# Pierrot, 2024-03-12

import pandas as pd
import numpy as np 
import requests
from io import StringIO
from astropy.table import Table

def fetch_astrometric_database(main_url, pending_url):
    # Fetch main table
    main_request = requests.get(main_url)
    main_dataframe = pd.read_csv(StringIO(main_request.text))

    # Fetch pending table
    pending_request = requests.get(pending_url)
    pending_dataframe = pd.read_csv(StringIO(pending_request.text))

    # Merge main and pending tables, removing duplicates
    astrom_dataframe = pd.concat([main_dataframe, pending_dataframe]).drop_duplicates()

    # Convert to Astropy table
    astrom_table = Table.from_pandas(astrom_dataframe)

    return astrom_table

def check_objects_in_database(csv_file, database):
    df = pd.read_csv(csv_file)
    objects_not_in_database = [] # List of objects not in the database

    # Check if objects are in the database
    for obj_name in df['Object Name']:
        found = False
        for record in database:
            if obj_name.lower() == record['OBJNAME'].lower() or obj_name.lower() in [alias.lower() for alias in record['ALIASES']]:
                found = True
                break
        if not found:
            objects_not_in_database.append(obj_name)

    return objects_not_in_database

MAIN_URL = ('https://docs.google.com/spreadsheets/d/'
            '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/'
            'export?format=csv&gid=0')
PENDING_URL = ('https://docs.google.com/spreadsheets/d/'
               '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/'
               'export?format=csv&gid=623506317')


astrometric_database = fetch_astrometric_database(MAIN_URL, PENDING_URL)
csv_file = "your_file.csv"  # Provide the path to your CSV file

objects_not_in_database = check_objects_in_database(csv_file, astrometric_database)
print("Objects not in the astrometric database:")
for obj in objects_not_in_database:
    print(obj)
