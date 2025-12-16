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
TABLENAME = 'students'
# test this one first
MODE = 'sqlite3'
# test this one second
# MODE = 'mysql.connect'
# test this one third
# MODE = 'sqlalchemy'
# path for sqlite3 database
SQLITE_DB_PATH = '/home/cook/index.db'
# -----------------------------------------------------------------------------
# define the command
QUERY_CREATE = f"""
CREATE TABLE IF NOT EXISTS {TABLENAME} (
    id INTEGER,
    name TEXT
);
"""

QUERY_INSERT = f"""
INSERT INTO {TABLENAME} (id, name)
VALUES (1, 'Alice');
"""

QUERY_SELECT = f"""
SELECT * FROM {TABLENAME};
"""


QUERIES = [QUERY_CREATE, QUERY_INSERT, QUERY_SELECT]
FETCH = [False, False, True]

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


def execute_sqlite_connect(_query: str, db_path: str, fetch=True):
    import sqlite3
    # connect to SQLite database file
    conn = sqlite3.connect(db_path)
    # get a cursor
    cursor = conn.cursor()
    # execute query
    cursor.execute(str(_query))

    if fetch:
        # fetch results
        rows = cursor.fetchall()
        return rows
    else:
        conn.commit()
    # cleanup
    cursor.close()
    conn.close()



# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------

    for mode in ['sqlite3', 'mysql.connect', 'sqlalchemy']:

        if MODE == 'sqlite3':
            # execute query
            try:
                for it, query in enumerate(QUERIES):
                    output = execute_sqlite_connect(query, SQLITE_DB_PATH,
                                                    fetch=FETCH[it])
            except Exception as e:
                print(e)
        elif MODE == 'mysql.connect':
            # execute query
            try:
                for query in QUERIES:
                    output = execute_mysql_connect(query)
            except Exception as e:
                print(e)
        elif MODE == 'sqlalchemy':
            # execute query
            try:
                for query in QUERIES:
                    output = execute_sqlalchemy(query)
            except Exception as e:
                print(e)
        else:
            raise ValueError(f'Unsupported mode: {MODE}')
        # print result
        print(f'{mode}: Found {len(output)} rows in {TABLENAME}')


# =============================================================================
# End of code
# =============================================================================
