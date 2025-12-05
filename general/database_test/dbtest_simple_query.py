#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Very basic code to emulate the way apero connects to the MySQL databases
If this works APERO should work.

Created on 2025-06-12 at 10:21

@author: cook
"""
# =============================================================================
# Define variables
# =============================================================================
HOST = 'rali.astro.umontreal.ca'
USER = 'nirps'
PASSWD = 'Covid19!'
DBNAME = 'nirps'
TABLENAME = 'findex_nirps_he_online_db'
# test this one first
MODE = 'mysql.connect'
# test this one second
MODE = 'sqlalchemy'
# -----------------------------------------------------------------------------
# define the command
QUERY = f'SELECT COUNT(*) FROM {TABLENAME}'


# =============================================================================
# Define functions
# =============================================================================
def execute_mysql_connect(_query: str):
    # import mysql
    import mysql.connector as mysql
    # connect and create a MySQL connection object
    conn = mysql.connect(host=HOST, user=USER, passwd=PASSWD,
                         database=DBNAME,
                         connection_timeout=3600)
    # get an entry point from the connection
    cursor = conn.cursor()
    # cmd
    command = str(_query)
    # run the cursor
    cursor.execute(command)
    # get the result
    result = cursor.fetchall()
    rows = list(result)
    cursor.close()
    conn.close()
    return rows


def execute_sqlalchemy(_query: str):
    # import sqlalchemy
    import sqlalchemy
    # create a database engine for sqlalchemy
    dpath = 'mysql+mysqlconnector://{0}:{1}@{2}/{3}'
    dargs = [USER, PASSWD, HOST, DBNAME]
    # connect and create a sqlalchemy engine (connection)
    engine = sqlalchemy.create_engine(dpath.format(*dargs), pool_pre_ping=True)
    # connect to database
    db = engine.connect()
    # cmd
    command = sqlalchemy.text(_query)
    # run the cursor
    result = db.execute(command)

    rows = list(result.fetchall())
    db.close()
    return rows

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    if MODE == 'mysql.connect':
        # execute query
        output = execute_mysql_connect(QUERY)
    elif MODE == 'sqlalchemy':
        # execute query
        output = execute_sqlalchemy(QUERY)
    else:
        raise ValueError(f'Unsupported mode: {MODE}')
    # print result
    print(f'Found {output} rows in {TABLENAME}')


# =============================================================================
# End of code
# =============================================================================
