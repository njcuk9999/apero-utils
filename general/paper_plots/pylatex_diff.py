#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Takes an original latex document and an updated one merges them and
takes the difference

requires latexpand and latexdiff  (and latex)

Created on 2022-10-07 at 15:05

@author: cook
"""
import os

# =============================================================================
# Define variables
# =============================================================================
# path to original
PATH1 = '/scratch2/spirou/misc/APERO_paper_review/original'
# path to update
PATH2 = '/scratch2/spirou/misc/APERO_paper_review/update'
# path to diff
PATH3 = '/scratch2/spirou/misc/APERO_paper_review/diff'
# name of main tex file
MAIN_TEX_FILE = 'main.tex'
# name of expanded tex file (created here)
ALL_TEX_FILE = 'all.tex'
# name of the diff tex file
DIFF_TEX_FILE = 'diff.tex'
# path to latex expand
LATEX_EXPAND = '/scratch/bin/latex_home/latexpand/latexpand'
# path to latex diff
LATEX_DIFF = '/scratch/bin/latex_home/latexdiff/latexdiff'

#==============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # get path to expanded file
    all_path_1 = os.path.join(PATH1, ALL_TEX_FILE)
    all_path_2 = os.path.join(PATH2, ALL_TEX_FILE)

    # get path to diff file
    diff_path = os.path.join(PATH3, DIFF_TEX_FILE)

    # get cwd
    current_dir = os.getcwd()

    # use latex expand to merge all tex files in original
    print('=' * 50 + '\nMERGING ORIGINAL\n' + '=' * 50)
    os.chdir(PATH1)

    cmd = f'{LATEX_EXPAND} {MAIN_TEX_FILE} -o {ALL_TEX_FILE}'
    print('\t', cmd)
    os.system(cmd)

    # use latex expand to merge all tex files in update
    print('=' * 50 + '\nMERGING UPDATE\n' + '=' * 50)
    os.chdir(PATH2)
    cmd = f'{LATEX_EXPAND} {MAIN_TEX_FILE} -o {ALL_TEX_FILE}'
    print('\t', cmd)
    os.system(cmd)

    # take the diff of the two files
    print('=' * 50 + '\nCREATE DIFF\n' + '=' * 50)
    cmd = f'{LATEX_DIFF} {all_path_1} {all_path_2} > {diff_path}'
    print('\t', cmd)
    os.system(cmd)

    # change back to cwd
    os.chdir(current_dir)

# =============================================================================
# End of code
# =============================================================================
