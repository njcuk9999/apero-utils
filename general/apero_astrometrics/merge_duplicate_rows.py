#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-02-29 at 08:47

@author: cook
"""
import os
from typing import Any, Dict, List, Tuple

import gspread_pandas as gspd
import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.time import Time

from apero.base import base
from apero.core import constants
from apero.core.core import drs_database
from apero.tools.module.setup import drs_installation
from apero.tools.module.database import drs_astrometrics

# =============================================================================
# Define variables
# =============================================================================
BAD_OBJECTS_PATH = ('/scratch2/spirou/drs-data/spirou_minidata2_07286_jupiter/'
                    'other/astrometrics/')
BAD_OBJECTS_FILE = 'bad_objects_2024-02-27T23:27:01.017.yaml'
OUTPUT_DATABASE_FILE = 'astrometric_database.csv'


# define sheet name
GSHEET_NAME = 'main_list'

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def get_google_sheet() -> pd.DataFrame:
    # load apero parameters
    params = constants.load()
    # add gspread directory and auth files
    drs_astrometrics.gsp_setup()
    # define the sheet id and sheet name (pending)
    sheet_id = params['OBJ_LIST_GOOGLE_SHEET_URL']
    # load google sheet instance
    google_sheet = gspd.spread.Spread(sheet_id)
    # convert google sheet to pandas dataframe
    dataframe = google_sheet.sheet_to_df(index=0, sheet=GSHEET_NAME)
    # return dataframe
    return dataframe


def get_associated_objects(problems: List[str]) -> List[str]:
    # list of associated objects
    associated_objects = []
    # loop around problems
    for problem in problems:
        # look for "OBJNAME:" in problems
        if ' OBJNAME:' in problem:
            # get the object name
            objname = problem.split(' OBJNAME:')[1].strip(']').strip('[').strip(' ')
            # add to associated objects
            associated_objects.append(objname)
        # look for "ALIAS:" in problems
        elif ' ALIAS:' in problem:
            # get the object name
            objname = problem.split(' ALIAS:')[1].strip(']').strip('[').strip(' ')
            # add to associated objects
            associated_objects.append(objname)
        # look for objname between "too close" and "sep"
        elif 'too close' in problem:
            # get the object name
            objname = problem.split('too close')[1].split('sep')[0].strip(' ').strip(' ')
            # add to associated objects
            associated_objects.append(objname)
    # only keep unique objnames
    associated_objects = list(set(associated_objects))
    # return these object names
    return associated_objects


def merge_duplicate_rows(problem_table: Table) -> Tuple[bool, Dict[str, Any]]:
    question_cols = ['OBJNAME', 'RA_DEG', 'DEC_DEG', 'EPOCH', 'PMRA', 'PMDE',
                     'PLX', 'RV', 'TEFF', 'SP_TYPE']

    # make transpose of table
    question_table = dict()
    # first column is the column name
    question_table['COLUMN'] = question_cols
    # columns are the objnames and the rows are the columns
    for row, objname in enumerate(problem_table['OBJNAME']):
        question_table[objname] = []

        # loop around columns
        for column in question_cols:
            # append to question table
            question_table[objname].append(problem_table[column][row])

    # make a table from the question table
    question_table = Table(question_table)

    # store results
    results_table = dict()
    for column in problem_table.colnames:
        results_table[column] = None

    # ----------------------------------------------------------------------
    # print the question table
    question_table.pprint(max_lines=-1, max_width=-1)
    # ----------------------------------------------------------------------
    # ask questions
    # ----------------------------------------------------------------------
    # question1. Are these the same object? (yes/no)
    # -------------------------------------------------------------------------
    question1 = 'Are these the same object? (yes/no)'
    answer1 = drs_installation.ask(question1, dtype='YN')

    if not answer1:
        return False, dict()

    # -------------------------------------------------------------------------
    # question2. If yes which objname do we use
    # -------------------------------------------------------------------------
    question2 = 'Which objname do we use?'
    options2, optiondesc2 = [], []
    for it in range(len(problem_table)):
        options2.append(it + 1)
        optiondesc2.append('{0}. {1}'.format(it + 1, problem_table['OBJNAME'][it]))
    answer2 = drs_installation.ask(question2, dtype=int,
                                   options=options2,
                                   optiondesc=optiondesc2)

    results_table['OBJNAME'] = problem_table['OBJNAME'][answer2 - 1]
    results_table['ORIGINAL_NAME'] = problem_table['ORIGINAL_NAME'][answer2 - 1]
    # combine all alises + objname + original names from all rows
    aliases = []
    for row in range(len(problem_table)):
        aliases.append(problem_table['OBJNAME'][row])
        aliases.append(problem_table['ORIGINAL_NAME'][row])
        aliases += problem_table['ALIASES'][row].split('|')
    # make unique
    aliases = set(aliases)
    # remove blank entries
    for bad_char in ['', ' ']:
        if bad_char in aliases:
            aliases.remove(bad_char)
    # make a list
    aliases = list(aliases)
    # push into aliases
    results_table['ALIASES'] = '|'.join(aliases)
    # -------------------------------------------------------------------------
    # ask the other questions
    # -------------------------------------------------------------------------
    # set up other questions
    question_names = ['RA/Dec', 'EPOCH', 'PMRA/PMDE', 'PLX', 'RV',
                      'TEFF', 'SP_TYPE']
    question_cols = [('RA_DEG', 'DEC_DEG'), 'EPOCH',
                     ('PMRA', 'PMDE'), 'PLX', 'RV', 'TEFF', 'SP_TYPE']
    question_scols = [('RA_SOURCE', 'DEC_SOURCE'), None,
                      ('PMRA_SOURCE', 'PMDE_SOURCE'), 'PLX_SOURCE', 'RV_SOURCE',
                      'TEFF_SOURCE', 'SP_SOURCE']
    # loop around questions
    for question_it in range(len(question_names)):

        # get the columns
        question_name = question_names[question_it]
        question_col = question_cols[question_it]
        question_scol = question_scols[question_it]

        same = True
        # deal with col not being a tuple (force to single element list)
        if not isinstance(question_col, tuple):
            question_col = [question_col]
            question_scol = [question_scol]

        # check if all are the same
        for c_it, col in enumerate(question_col):
            same &= len(set(np.array(problem_table[col]))) == 1
            if question_scol[c_it] is not None:
                scol = question_scol[c_it]
                same &= len(set(np.array(problem_table[scol]))) == 1
        # if same then we can use any column
        if same:
            for c_it, col in enumerate(question_col):
                results_table[col] = problem_table[col][0]
                if question_scol[c_it] is not None:
                    scol = question_scol[c_it]
                    results_table[scol] = problem_table[scol][0]
        else:
            # construct the question
            question = 'Which {0} do we use?'.format(question_name)
            # loop around rows to construct options
            options, optiondesc = [], []
            for it in range(len(problem_table)):
                options.append(it + 1)
                # deal with value
                value = ''
                # deal with having multiple columns
                for c_it, col in enumerate(question_col):
                    if question_scol[c_it] is None:
                        vargs = [col,
                                 problem_table[col][it]]
                        value += '{0}={1}'.format(*vargs)
                    else:
                        vargs = [col, problem_table[col][it],
                                 question_scol[c_it],
                                 problem_table[question_scol[c_it]][it]]

                        value += '{0}={1} ({2}={3}) '.format(*vargs)
                # populate the option description
                optionargs = [it + 1, value, problem_table['OBJNAME'][it]]
                optiondesc.append('{0}. {1} [{2}]'.format(*optionargs))
            # ask the question
            answer = drs_installation.ask(question, dtype=int,
                                          options=options,
                                          optiondesc=optiondesc)
            # push into results table
            for c_it, col in enumerate(question_col):
                results_table[col] = problem_table[col][answer - 1]
                # deal with source column
                if question_scol[c_it] is not None:
                    scol = question_scol[c_it]
                    results_table[scol] = problem_table[scol][answer - 1]
    # -------------------------------------------------------------------------
    # deal with other columns
    # -------------------------------------------------------------------------
    results_table['VSINI'] = ''
    results_table['VSINI_SOURCE'] = ''
    results_table['VSINI_ERR'] = ''
    results_table['u_'] = ''
    results_table['g_'] = ''
    results_table['r_'] = ''
    results_table['i_'] = ''
    results_table['z_'] = ''
    results_table['J'] = ''
    results_table['H'] = ''
    results_table['K'] = ''

    nargs = ['|'.join(problem_table['OBJNAME']), Time.now().iso]
    note = 'Merged from multiple entries [{0}] on {1}'
    results_table['NOTES'] = note.format(*nargs)

    results_table['USED'] = 1

    results_table['KEYWORDS'] = '|'.join(set(np.array(problem_table['KEYWORDS'])))


    # -------------------------------------------------------------------------
    return True, results_table


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # construct the yaml file path
    bad_yaml_file = os.path.join(BAD_OBJECTS_PATH, BAD_OBJECTS_FILE)
    # load the yaml file
    bad_objects = base.load_yaml(bad_yaml_file)

    # store done targets
    done_targets = []
    # get the apero parameters
    params = constants.load()
    # get google sheet as a dataframe
    # get full database
    astrometric_dataframe = get_google_sheet()

    # ---------------------------------------------------------------------
    # # fix for bad K column
    # # TODO: remove this when fixed in google sheet
    # source1 = astrometric_dataframe['SP_SOURCE']
    # source2 = astrometric_dataframe['K']
    # source_comb = []
    # for s1, s2 in zip(source1, source2):
    #     if s1 == s2:
    #         source_comb.append(s1)
    #     elif s1 != '':
    #         source_comb.append(s1)
    #     else:
    #         source_comb.append(s2)
    # astrometric_dataframe['SP_SOURCE'] = source_comb
    # astrometric_dataframe['K'] = [''] * len(source_comb)
    # ----------------------------------------------------------------------
    # loop around bad objects
    for b_it, bad_object in enumerate(bad_objects):

        # print header
        print('\n\n\n\n' + '=' * 50)
        pargs = [bad_object, b_it + 1, len(bad_objects)]
        print('Processing {0} [{1}/{2}]'.format(*pargs))
        print('=' * 50)
        print('\n\n')

        # skip if in done targets
        if bad_object in done_targets:
            print('\t Skipping {0} (already done)'.format(bad_object))
            continue

        # get the problem list for this object
        problem_list = bad_objects[bad_object]

        # try to find other objects associated with this object
        associated_objects = get_associated_objects(problem_list)

        if len(associated_objects) == 0:
            print('No associated objects found for {0}'.format(bad_object))
            continue

        # object list
        obj_list = [bad_object] + associated_objects

        # get mask for the objects list
        obj_mask = astrometric_dataframe['OBJNAME'].isin(obj_list)
        # convert to astropy table
        obj_table = Table.from_pandas(astrometric_dataframe[obj_mask])

        # attempt to merge entries into single row
        #   ask the user to deal with conflicts
        merged, row_dict = merge_duplicate_rows(obj_table)

        # add bad objects to done targets
        done_targets += obj_list

        # do not continue if not merged
        if not merged:
            continue

        # remove old rows from the astromeric dataframe containing these object
        # names
        astrometric_dataframe = astrometric_dataframe[~obj_mask]

        # convert dict into dataframe
        row_dataframe = pd.DataFrame([row_dict])

        # add new row to the dataframe using row_dict
        astrometric_dataframe = pd.concat([astrometric_dataframe,
                                           row_dataframe],
                                          ignore_index=True)
        # reset index
        astrometric_dataframe.reset_index()

    # -------------------------------------------------------------------------
    # construct the output path
    output_database_file = os.path.join(BAD_OBJECTS_PATH, OUTPUT_DATABASE_FILE)
    # print progress
    print('\n\n\n\n' + '=' * 50)
    print('Saving to {0}'.format(output_database_file))
    print('=' * 50)
    # finally save the dataframe to a csv (for opening then save manually to
    # google sheet)
    astrometric_dataframe.to_csv(output_database_file, index=False)


# =============================================================================
# End of code
# =============================================================================
