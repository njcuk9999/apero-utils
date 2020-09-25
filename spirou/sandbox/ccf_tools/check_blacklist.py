import os
import numpy as np


def check_blacklist(file_list):
    # pass a file list and remove files that have their odometer
    # included in the blacklisted list available on GitHub

    # adress where the blacklist is kept
    tbl_url = 'https://raw.githubusercontent.com/njcuk9999/apero-utils/master/spirou/sandbox/ccf_tools/blacklist_odo.tbl'

    # if a previous version is on disk, then delete
    fname = tbl_url.split('/')[-1]
    if os.path.isfile(fname):
        os.system('rm ' + fname)

    # get the file from github
    os.system('wget '+tbl_url)

    # append lines read in the file
    lines = []
    with open(fname, 'r') as reader:
        lines.append(reader.read())

    # remove file
    os.system('rm '+fname)

    # join lines
    lines = ''.join(lines)

    # reformat
    lines = np.array(lines.split('\n'))
    keep = np.ones_like(lines,dtype = bool)

    # remove lines that are comments or too short
    for i in range(len(lines)):
        if len(lines[i])<7: # too short to be valid
            keep[i] = False
            continue

        if '#' == lines[i][0]: # is a comment
            keep[i] = False
            continue

        if ',' not in lines[i]: # does not have a comment
            keep[i] = False
            continue

        tmp = lines[i].split(',')[0]
        tmp = ''.join(tmp.split(' '))

        if len(tmp) !=7: # not the length of an odometer
            keep[i] = False
            continue

        lines[i] = tmp

    # keep only lines that math the proper format
    bad_odo = lines[keep]

    # find which files should be kept (they are not in the bad odo list)
    good_file = np.ones_like(file_list,dtype = bool)

    # check if odometer is in file name
    for i in range(len(file_list)):
        for j in range(len(bad_odo)):
            if bad_odo[j] in file_list[i]:
                good_file[i] = False

    # remove bad files from list

    if False in good_file:
        print('We remove the following files : ')
        bad_files = file_list[~good_file]
        for bad_file in bad_files:
            print('\t{0}'.format(bad_file))
        print()
    else:
        print('No file listed in the blacklist')

    return file_list[good_file]