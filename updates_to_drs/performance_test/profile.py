#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-05-17

@author: cook
"""
import psutil
import subprocess

# =============================================================================
# Define variables
# =============================================================================
APERO_COMMAND = 'apero_processing'

NETSTAT_CMD = 'netstat -nt | wc -l'

# =============================================================================
# Define functions
# =============================================================================
def get_system_performance():

    processes = []

    for proc in psutil.process_iter():

        if proc.status() != 'running':
            continue

        for cmd in proc.cmdline():
            if APERO_COMMAND in cmd:
                processes.append(proc.cmdline())
                break

        if proc.username() == 'spirou' and proc.name() =='python':
            print('{0} CPU={1:.2f} MEM={2:.2f}'.format(proc.name(), proc.cpu_percent(), proc.memory_info()[0]/2**30))


def get_database_performance():

    # -------------------------------------------------------------------------
    # net stat command (number of TCP connections)
    # -------------------------------------------------------------------------
    netstatcmd = subprocess.Popen(NETSTAT_CMD, shell=True,
                                  stdout=subprocess.PIPE)
    nconnect = netstatcmd.stdout.read().decode().rstrip()
    try:
        nconnect = int(nconnect)
    except:
        nconnect = 0

    return nconnect




# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # print hello world
    print('Hello World')

# =============================================================================
# End of code
# =============================================================================
