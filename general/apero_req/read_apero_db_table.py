#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-05-01 at 08:48

@author: cook
"""
from typing import Dict
from astropy.table import Table, Column
from sqlalchemy import create_engine, text
import pandas as pd
import numpy as np


# =============================================================================
# Define variables
# =============================================================================
DB_TABLE_NAMES = ['astrom', 'calib', 'findex', 'lang', 'log', 'reject', 'tellu']

DBFILES = dict()
DBFILES['spirou'] = 'spirou_db_table.txt'
DBFILES['nirps'] = 'nirps_db_table.txt'

DBLOGIN = dict()
DBLOGIN['nirps'] = dict()
DBLOGIN['nirps']['type'] = 'mysql+pymysql'
DBLOGIN['nirps']['user'] = 'nirps'
DBLOGIN['nirps']['host'] = 'rali.astro.umontreal.ca'
DBLOGIN['nirps']['password'] = 'Covid19!'
DBLOGIN['nirps']['dbname'] = 'nirps'
DBLOGIN['spirou'] = dict()
DBLOGIN['spirou']['type'] = 'mysql+pymysql'
DBLOGIN['spirou']['user'] = 'spirou'
DBLOGIN['spirou']['host'] = 'cosmos.astro.umontreal.ca'
DBLOGIN['spirou']['password'] = 'Covid19!'
DBLOGIN['spirou']['dbname'] = 'spirou'


# =============================================================================
# Define functions
# =============================================================================
def run_query(dblogin: Dict[str, str], query: str) -> pd.DataFrame:
    # create the url for the database engine
    url = '{type}://{user}:{password}@{host}/{dbname}'.format(**dblogin)
    # create the database engine
    engine = create_engine(url)
    # read the sql with pandas
    with engine.connect() as conn:
        result = pd.read_sql(query, conn)
    # return result
    return result


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":

    for instrument in DBFILES.keys():

        dbfilename = DBFILES[instrument]

        # load db file
        with open(dbfilename) as dbfile:
            lines = dbfile.readlines()

        # Keep only lines with '|'
        clean_lines = [line for line in lines if line.strip().startswith('|')]

        # Join and read using astropy
        clean_text = ''.join(clean_lines)
        table = Table.read(clean_text, format='ascii.fixed_width', delimiter='|')

        # add column table type
        table_type = []
        profile_names = []
        findex_dbs = dict()

        for row in range(len(table)):

            found = False
            table_name = table['table_name'][row]

            for db_table_name in DB_TABLE_NAMES:
                prefix = f'{db_table_name}_'

                if table_name.startswith(prefix):
                    table_type.append(db_table_name)

                    profile_name = table_name.split(prefix)[-1]

                    while profile_name.endswith('_db') or profile_name.endswith('_DB'):
                        profile_name = profile_name.strip('_db').strip('._DB')

                    profile_names.append(profile_name)

                    found = True
                    break

            if not found:
                table_type.append('NULL')
                profile_names.append('NULL')

            if 'findex' in table_name:
                findex_dbs[profile_names[row]] = table_name

        # push into table
        del table['table_name']
        del table['table_schema']
        table.add_column(Column(profile_names, name='profile'), index=0)
        table.add_column(Column(table_type, name='type'), index=1)

        # get the number of raw files for each findex
        findex_counts = dict()

        for findex in findex_dbs:

            findex_table = findex_dbs[findex]

            condition = f'SELECT COUNT(*) as num_raw FROM {findex_table} WHERE BLOCK_KIND="raw"'

            count_result = run_query(DBLOGIN[instrument], condition)

            findex_counts[findex] = int(count_result['num_raw'][0])

        # get the sort order from findex_counts
        # TODO: this doesn't work but I don't know why
        fcount_keys = np.array(list(findex_counts.keys()))
        fcount_values = np.array(list(findex_counts.values()))

        # sort dictionary
        isort = np.argsort(fcount_values)[::-1]
        fcount_keys = fcount_keys[isort]
        fcount_sort = list(range(len(fcount_keys)))
        sort_order = dict(zip(list(fcount_keys), list(fcount_sort)))

        # turn the counts into a column
        raw_counts = []
        sort_col = []
        for row in range(len(table)):
            profile_name = table['profile'][row]
            if profile_name == 'NULL':
                raw_counts.append('NULL')
                sort_col.append(-1)
            else:
                raw_counts.append(findex_counts[profile_name])
                sort_col.append(sort_order[profile_name])

        # add raw counts to table
        table.add_column(Column(raw_counts, name='n_raw'), index=2)
        table.add_column(Column(sort_col, name='sort_col'), index=3)

        # remove null profile names
        table = table[table['profile'] != 'NULL']

        # convert n_raw to int
        table['n_raw'] = np.array(table['n_raw']).astype(int)

        # remove profiles
        remove_profiles = []
        for row in range(len(table)):
            profile_name = table['profile'][row]
            if table['size_mb'][row] == 0:
                remove_profiles.append(profile_name)
                continue
            if table['n_raw'][row] == 0:
                remove_profiles.append(profile_name)
                continue

        # keep only unique profiles
        remove_profiles = np.unique(remove_profiles)

        for remove_profile in remove_profiles:
            table = table[table['profile'] != remove_profile]


        # sort by sort order then profile then type
        table.sort(['sort_col', 'profile', 'type'])
        # remove the sort_col column
        del table['sort_col']
        # round the size_mb column to the nearest Mb
        table['size_mb'] = np.round(table['size_mb'], 2)

        # Escape underscores
        table['profile'] = [s.replace('_', r'\_') for s in table['profile']]

        # write to latex file
        table.write(dbfilename.replace('.txt', '.tex'), format='latex',
                    overwrite=True)

# =============================================================================
# End of code
# =============================================================================
