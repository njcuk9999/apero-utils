#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-09-10 at 12:00

@author: cook
"""
from collections import Counter
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy import units as uu
from astropy.table import Table
from astropy.time import Time

# =============================================================================
# Define variables
# =============================================================================
# The URL to google (must have the "sheet_id" and "gid" parts)
GOOGLE_URL = 'https://docs.google.com/spreadsheets/d/{sheet_id}/export?format=csv&gid={gid}'
# allocation sheet info
ALLOCATION_ID = '1s116aabnMH0zJ5YbXrGBWIVP17lmYXl6tYzLteoSEAc'
ALLOCATION_SHEET = '0'
# checks sheet info
CHECKS_ID = '1zvU_XFA1ZOJE111qZKiav7v6ptYveWpDkMHjdhiN06M'
# nirps he sheet info
RAW_SHEET = dict(NIRPS_HA='1017212188', NIRPS_HE='1741354113')
RED_SHEET = dict(NIRPS_HA='279548254', NIRPS_HE='541545005')
COMM_SHEET = dict(NIRPS_HA='346030560', NIRPS_HE='68390787')
# -----------------------------------------------------------------------------
INSTRUMENTS = ['NIRPS_HE', 'NIRPS_HA']
SHORT = ['he', 'ha']
# Start from this date
ALLOCATION_START = '2023-06-01'
# Start plot from this date
COMMENT_START = '2024-01-01'
# Name people who have left
EX_MEMBERS = ['Fred', 'Yuri', 'Rose', 'Olivia']
# Ignore these columns when checking for False in check tables
IGNORE_COLS = ['obsdir', 'date', 'BLANK']


# =============================================================================
# Define functions
# =============================================================================
def read_google_sheet_csv(sheet_id: str, gid: str) -> Table:
    """
    This function reads a Google sheet and returns the content as an
    astropy table

    :param sheet_id:
    :param gid:
    :return: the astropy table
    """
    # Construct the URL
    google_url = GOOGLE_URL.format(sheet_id=sheet_id, gid=gid)
    # print that we are getting table from url
    print(f'Getting table from: {google_url}')
    # Send a GET request to the URL
    response = requests.get(google_url)
    # read the csv file
    table = Table.read(response.text, format='ascii.csv')
    # return the astropy table
    return table


def allocation_histogram(allocation_table: Table):

    # filter table (online prodedure only)
    tasks = allocation_table['Task']
    mask = tasks == 'online procedure'
    # remove who masked rows
    mask &= ~allocation_table['Who'].mask
    # filter the table
    allocation_table = allocation_table[mask]

    # only count days from start date until today
    dates = Time(allocation_table['Date Start'])

    datemask = dates >= Time(ALLOCATION_START)
    datemask &= dates < Time.now()

    # get list of members
    who = allocation_table['Who'][datemask]
    members, exmembers = [], []
    for person in set(who):
        if person in EX_MEMBERS:
            exmembers.append(person)
        else:
            members.append(person)
    # Generate a unique color for each name
    colors = ['red', 'blue', 'orange', 'purple', 'green', 'yellow', 'brown',
              'pink', 'gray', 'cyan', 'magenta', 'olive', 'lime', 'teal'] * 10
    # setup figure
    plt.figure(figsize=(20, 20))
    frames = [plt.subplot2grid((1, 4), (0, 0), colspan=3),
              plt.subplot2grid((1, 4), (0, 3), colspan=1)]
    # store max number of times someone was allocated
    max_entry = 0
    # define groups
    group_names = ['Members', 'Previous Members']
    # loop around members and ex members
    for it, group in enumerate([members, exmembers]):

        groupmask = np.in1d(allocation_table['Who'], group)

        # get the name counts
        name_counts = Counter(allocation_table['Who'][datemask & groupmask])

        # Get the maximum entry
        max_entry = max(max_entry, max(name_counts.values()))

        # Plot the histogram
        frames[it].bar(name_counts.keys(), name_counts.values(),
                color=colors[:len(name_counts)])
        # Add labels and title
        frames[it].set(ylabel='Number of babysitting times',
                       title=group_names[it])

    # set y limits for all frames
    for frame in frames:
        frame.set_ylim(0, max_entry + 1)
    # put the y-axis label and ticks on the right side
    frames[1].yaxis.tick_right()
    frames[1].yaxis.set_label_position('right')

    # construct the title
    targs = [ALLOCATION_START, Time.now().iso.split()[0]]
    title = 'Allocation Histogram\n From {0} to {1}'.format(*targs)
    plt.suptitle(title)
    # save png
    plt.savefig('results/allocation_histogram.png')
    # Show the plot
    plt.show()
    # return the allocation table
    return allocation_table


def trace_falses(raw_table: Table, red_table: Table, comm_table: Table
                 ) -> Tuple[List[str], List[List[str]], List[List[str]]]:

    # deinfe which tables to test
    tables = [raw_table, red_table]
    table_types = ['raw', 'red']
    # storage for dates and tests
    dates = []
    raw_tests = []
    red_tests = []

    # loop around raw and red table
    for table_it, table in enumerate(tables):
        # get a list of check columns
        check_cols = []
        for col in table.columns:
            if col not in IGNORE_COLS:
                check_cols.append(col)
        # find any date that has a False that is not in the comments
        for row, date in enumerate(table['obsdir']):
            # skip dates before start
            if Time(date) < Time(ALLOCATION_START):
                continue
            # skip dates from yesterday onwards
            if Time(date) > Time.now() - 1 * uu.day:
                continue
            # keep track of the check cols (assume they are all true to start)
            passed, fails = True, []
            # loop around check cols
            for col in check_cols:
                if table[col][row] == 'FALSE':
                    passed = False
                    fails.append(col)
            # if some tests have not passed and not been mentioned in
            #  the comments table then add to the list
            if not passed and date not in comm_table['obsdir']:
                # need to deal with date already in list
                if date in dates:
                    index = dates.index(date)
                    if table_types[table_it] == 'raw':
                        raw_tests[index] = list(set(raw_tests[index] + fails))
                    else:
                        red_tests[index] = list(set(red_tests[index] + fails))
                # otherwise we add a new entry
                else:
                    dates.append(date)
                    if table_types[table_it] == 'raw':
                        raw_tests.append(list(set(fails)))
                        red_tests.append([])
                    else:
                        raw_tests.append([])
                        red_tests.append(list(set(fails)))
    # return dates and tests which fail and have no comment to match
    return dates, raw_tests, red_tests

def who_was_allocated(dates: List[str], allocation_table: Table) -> List[str]:

    # get the start and end dates in the allocation table
    start_date = Time(allocation_table['Date Start'])
    end_date = Time(allocation_table['Date End'])
    # storage for who was allocated
    allocated = []
    # loop around dates
    for date in dates:
        # find who was allocated to this date
        mask = (start_date <= Time(date)) & (end_date >= Time(date))
        allocated.append(str(allocation_table['Who'][mask][0]))
    # return the allocated dictionary
    return allocated


def comment_histogram(stat_table: Table):

    # only count days from start date until today
    dates = Time(stat_table['obsdir'])

    datemask = dates >= Time(COMMENT_START)
    datemask &= dates < Time.now()

    # get list of members
    who = stat_table['who'][datemask]
    members, exmembers = [], []
    for person in set(who):
        if person in EX_MEMBERS:
            exmembers.append(person)
        else:
            members.append(person)
    # Generate a unique color for each name
    colors = ['red', 'blue', 'orange', 'purple', 'green', 'yellow', 'brown',
              'pink', 'gray', 'cyan', 'magenta', 'olive', 'lime', 'teal'] * 10
    # setup figure
    plt.figure(figsize=(20, 20))
    frames = [plt.subplot2grid((1, 4), (0, 0), colspan=3),
              plt.subplot2grid((1, 4), (0, 3), colspan=1)]
    # store max number of times someone was allocated
    max_entry = 0
    # define groups
    group_names = ['Members', 'Previous Members']
    # loop around members and ex members
    for it, group in enumerate([members, exmembers]):

        groupmask = np.in1d(stat_table['who'], group)

        # get the name counts
        name_counts = Counter(stat_table['who'][datemask & groupmask])

        # Get the maximum entry
        max_entry = max(max_entry, max(name_counts.values()))

        # Plot the histogram
        frames[it].bar(name_counts.keys(), name_counts.values(),
                color=colors[:len(name_counts)])
        # Add labels and title
        frames[it].set(ylabel='Number of days comments not added on failure',
                       title=group_names[it])

    # set y limits for all frames
    for frame in frames:
        frame.set_ylim(0, max_entry + 1)
    # put the y-axis label and ticks on the right side
    frames[1].yaxis.tick_right()
    frames[1].yaxis.set_label_position('right')

    # construct the title
    targs = [COMMENT_START, Time.now().iso.split()[0]]
    title = 'Comment Histogram\n From {0} to {1}'.format(*targs)
    plt.suptitle(title)
    # save png
    plt.savefig('results/comment_histogram.png')
    # Show the plot
    plt.show()



def main():

    # read the allocation sheet
    allocation_table = read_google_sheet_csv(ALLOCATION_ID, ALLOCATION_SHEET)
    # -------------------------------------------------------------------------
    # run the allocation histogram
    allocation_table = allocation_histogram(allocation_table)
    # -------------------------------------------------------------------------
    # global storage
    stbl = dict()
    stbl['date'] = []
    stbl['who'] = []
    stbl['raw_ha'] = []
    stbl['red_ha'] = []
    stbl['raw_he'] = []
    stbl['red_he'] = []
    # -------------------------------------------------------------------------
    # loop around instruments
    for iit, instrument in enumerate(INSTRUMENTS):
        # read the raw sheet
        raw_table = read_google_sheet_csv(CHECKS_ID, RAW_SHEET[instrument])
        # read the red sheet
        red_table = read_google_sheet_csv(CHECKS_ID, RED_SHEET[instrument])
        # read the comm sheet
        comm_table = read_google_sheet_csv(CHECKS_ID, COMM_SHEET[instrument])
        # trace falses in the raw and red table (should have a matching comment)
        dates, raw_tests, red_tests = trace_falses(raw_table, red_table,
                                                   comm_table)
        # find out who was allocated to these dates
        allocated = who_was_allocated(dates, allocation_table)
        # get short key
        rawkey1 = f'raw_{SHORT[iit]}'
        redkey1 = f'red_{SHORT[iit]}'
        # flip the key
        jjt = (iit % 2) - 1
        rawkey2 = f'raw_{SHORT[jjt]}'
        redkey2 = f'red_{SHORT[jjt]}'
        # loop around rows
        for row in range(len(dates)):
            # add to stats table
            if dates[row] not in stbl['date']:
                stbl['date'].append(dates[row])
                stbl['who'].append(allocated[row])
                stbl[rawkey1].append(raw_tests[row])
                stbl[redkey1].append(red_tests[row])
                stbl[rawkey2].append([])
                stbl[redkey2].append([])
                continue
            # find row in stats table
            index = list(stbl['date']).index(dates[row])
            # update tests
            stbl[rawkey1][index] = list(set(stbl[rawkey1][index] + raw_tests[row]))
            stbl[redkey1][index] = list(set(stbl[redkey1][index] + red_tests[row]))
    # -------------------------------------------------------------------------
    # push into table
    stat_table = Table()
    stat_table['obsdir'] = stbl['date']
    stat_table['who'] = stbl['who']
    stat_table['raw_ha'] = [', '.join(x) for x in stbl['raw_ha']]
    stat_table['raw_he'] = [', '.join(x) for x in stbl['raw_he']]
    stat_table['red_ha'] = [', '.join(x) for x in stbl['red_ha']]
    stat_table['red_he'] = [', '.join(x) for x in stbl['red_he']]
    # filter table (since COMMENT_DATE)
    mask = Time(stat_table['obsdir']) >= Time(COMMENT_START)
    stat_table = stat_table[mask]
    # sort stat table by date
    stat_table = stat_table[np.argsort(stat_table['obsdir'])[::-1]]
    # write to file
    stat_table.write('results/test_table.csv', format='ascii.csv', overwrite=True)
    # -------------------------------------------------------------------------
    # plot the comment histogram
    comment_histogram(stat_table)
    # -------------------------------------------------------------------------
    # one last table: group by who and give days that need reporting
    who_dict = dict()
    who_dict['who'] = []
    who_dict['obsdirs'] = []

    for who in list(set(stat_table['who'])):
        mask = stat_table['who'] == who
        obsdirs = stat_table['obsdir'][mask]
        # sort obsdirs
        obsdirs = sorted(obsdirs)
        # push into dictionary
        who_dict['who'].append(who)
        who_dict['obsdirs'].append(', '.join(obsdirs))
    # push into table
    who_table = Table(who_dict)
    # write to file
    who_table.write('results/obsdir_table.csv', format='ascii.csv', overwrite=True)
    # -------------------------------------------------------------------------
    # return to code
    return stat_table, who_table












# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # run main code
    main()

# =============================================================================
# End of code
# =============================================================================
