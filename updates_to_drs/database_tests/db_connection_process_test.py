#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-03-30

@author: cook
"""
from astropy.time import Time
import getpass
import mysql.connector as mysql
import os
import time

# =============================================================================
# Define variables
# =============================================================================
# define workspace
WORKSPACE = '/spirou/cook/db_test'
# file for processes
PROCESS_FILE = os.path.join(WORKSPACE, 'process.txt')
# file for connections
CONNECT_FILE = os.path.join(WORKSPACE, 'connect.txt')
# get time now
now = Time.now()

# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # clear + make these files
    if os.path.exists(PROCESS_FILE):
        os.remove(PROCESS_FILE)
    os.system('echo "human_time, unix_time, count" >> {0}'.format(PROCESS_FILE))
    if os.path.exists(CONNECT_FILE):
        os.remove(CONNECT_FILE)
    os.system('echo "count" >> {0}'.format(CONNECT_FILE))
    # get user password
    password = getpass.getpass('Database password: ')
    # get a active database connection
    conn = mysql.connect(host='rali', user='spirou', passwd=password,
                         database='information_schema')
    # start iteration counter
    it = 0
    # try to catch Ctrl+C
    try:
        # loop until Control+C
        while 1:
            # get number of processes
            with conn.cursor() as cursor:
                command = 'SELECT COUNT(*) FROM PROCESSLIST WHERE DB="spirou";'
                cursor.execute(command)
                count = cursor.fetchall()[0][0]
            # get time now
            now = Time.now()
            # print statement
            print('{0}: Iteration {1}: Count = {2}'.format(now.iso, it + 1,
                                                           count))
            # push count to file
            os.system('echo "{0}, {1}, {2}" >> {3}'.format(now.fits, now.unix,
                                                           count, PROCESS_FILE))
            # get number of connections and push to file
            os.system('netstat -nt | wc -l >> {0}'.format(CONNECT_FILE))
            # count how many iterations we've done
            it += 1
            # sleep two seconds
            time.sleep(1)
    except KeyboardInterrupt:
        print('\nEnded')
        conn.close()
    finally:
        conn.close()

# =============================================================================
# End of code
# =============================================================================
