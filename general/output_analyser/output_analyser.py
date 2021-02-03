#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-02-02

@author: cook
"""
from astropy.io import fits
import numpy as np
import os
from pathlib import Path
from tqdm import tqdm
from typing import List, Union

from apero.core import constants

# =============================================================================
# Define variables
# =============================================================================
# set name
__NAME__ = 'output_analyser.py'
__DESCRIPTION__ = 'Produce a human readable summary of all outputs'

# directories to test
test_dirs = dict()
test_dirs['tmp'] = 'DRS_DATA_WORKING'
test_dirs['red'] = 'DRS_DATA_REDUC'
test_dirs['calib'] = 'DRS_CALIB_DB'
test_dirs['tellu'] = 'DRS_TELLU_DB'

# define file type to test
FILE_EXT = '.fits'

# only need to test one file of each type, define here how to find similar
# files (anything after this key that is identical is the same file)
COMMON_KEY = '_pp'
# output filename
output_file = 'output_analysis.txt'


# =============================================================================
# Define functions
# =============================================================================
class FitsExtension:
    def __init__(self, number: int, summary: tuple):
        self.number = number
        self.name = summary[0]
        self.version = summary[1]
        self.type = summary[2]
        self.cards = summary[3]
        self.dimensions = summary[4]
        self.format = summary[5]
        self.keys = []
        self.values = dict()
        self.comments = dict()

    def __str__(self) -> str:
        string_output = '\nFitsExt'
        outputs = self.display()
        for output in outputs:
            string_output += '\n' + output
        return string_output

    def __repr__(self) -> str:
        return self.__str__()

    def add_keys(self, header: Union[fits.Header, None] = None):
        # deal with no header
        if header is None:
            return
        # else add keys to header
        else:
            # loop around keys
            for key in header:
                # only save string representation of value
                self.values[key] = str(header[key])
                # add comment (copy as string)
                self.comments[key] = str(header.comments[key])
            # add keys list
            self.keys = list(header.keys())

    def display(self) -> List[str]:
        """
        Display a human readable version of extension

        :return:
        """
        output = []
        output += [' - EXTENSION {0}'.format(self.number)]
        output += ['         TYPE:   {0}'.format(self.type)]
        output += ['         NCARDS: {0}'.format(self.cards)]
        output += ['         DIM:    {0}'.format(self.dimensions)]
        return output


class FitsInfo:
    def __init__(self, filename):
        self.filename = Path(filename)
        self.basename = self.filename.name
        self.ext = self.filename.suffix
        self.extensions = []
        # update parameters using hdu._summary
        self.get_info()

    def __str__(self) -> str:
        return 'FitsFile[{0}]'.format(self.basename)

    def __repr__(self) -> str:
        return self.__str__()

    def get_info(self):
        # open hdu
        with fits.open(self.filename) as hdulist:
            # loop around each extension
            for idx, hdu in enumerate(hdulist):
                # take from fits.info()
                summary = hdu._summary()
                # make extension instance
                extension = FitsExtension(idx, summary)
                # add header keys
                if hdu.header is not None:
                    extension.add_keys(hdu.header)
                # add to extensions container
                self.extensions.append(extension)

    def pprint(self):
        string_output = self.__str__()
        outputs = self.display()
        for output in outputs:
            string_output += '\n' + output
        return string_output

    def display(self) -> List[str]:
        output = []
        # loop around extensions
        for extension in self.extensions:
            ext_outputs = extension.display()
            for ext_output in ext_outputs:
                output += ['\t' + ext_output]
        return output


def find_unique_files(test_path: str) -> List[str]:
    # store all files
    allfiles, common_texts = [], []
    # get all files in this path matching file extension
    for root, dirs, files in os.walk(test_path):
        for filename in files:
            # cond1: filename ends with extension
            cond1 = filename.endswith(FILE_EXT)
            # cond2: filename contains common key
            cond2 = COMMON_KEY in filename
            # add to all files if conditions are satisfied
            if cond1 and cond2:
                # append file
                allfiles.append(os.path.join(root, filename))
                # get comment text for all files
                common_suffix = filename.split(COMMON_KEY)[-1]
                # add to array
                cargs = [COMMON_KEY, common_suffix]
                common_texts.append('{0}{1}'.format(*cargs))
    # only keep one of each type of file
    unique_texts = np.unique(common_texts)
    unique_files = []
    # find first file with each unique text
    for unique_text in unique_texts:
        # get position
        pos = np.where(np.array(common_texts) == unique_text)[0][0]
        # add to unique files
        unique_files.append(allfiles[pos])
    # return unique_files
    return unique_files


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # get parameters
    params = constants.load('SPIROU')
    # storage of strings
    outstrings = ''
    # loop around directories
    for test_dir in test_dirs:
        # log progress
        print('\nProcessing {0} directory\n'.format(test_dir))
        # get path using params
        test_path = params[test_dirs[test_dir]]
        # add header to outstrings
        outstrings += '\n\n' + '=' * 100
        outstrings += '\n\t {0}: {1}'.format(test_dir, test_path)
        outstrings += '\n' + '=' * 100
        # find unique examples of each "COMMON KEY" file type
        ufiles = find_unique_files(str(test_path))
        # loop around each file and make fits instance
        for ufile in tqdm(ufiles):
            # make fits info instance
            fitsfile = FitsInfo(ufile)
            # store display
            outstrings += '\n\n' + fitsfile.pprint()
    # log progress
    print('\nSaving analysis to {0}'.format(output_file))
    # save outstrings to file
    with open(output_file, 'w') as ofile:
        ofile.write(outstrings)


# =============================================================================
# End of code
# =============================================================================
