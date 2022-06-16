#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-06-16

@author: cook
"""
from astropy.io import fits
import glob
import numpy as np
import os
from tqdm import tqdm


# =============================================================================
# Define variables
# =============================================================================
# set the instrument name
INSTRUMENT = 'NIRPS'
# set the instrument mode (None if no mode)
MODE = 'HA'
# whether to add additional printout
DEBUG = False
# whether to update raw files when mistake was found
UPDATE = True
# extra debug not to save file (as test)
TEST = True
# define for use in if statements below
SKIP_KEYS = dict()

if INSTRUMENT == 'NIRPS' and MODE == 'HA':
    REF_DIR = '/nirps_raw/nirps/apero-data/drs-misc/ref_files/NIRPS_HA'
    TEST_DIR = '/nirps_raw/nirps/raw-data/nirps_ha/'
    HEADER_KEYS = ['HIERARCH ESO DPR TYPE', 'HIERARCH ESO DPR CATG']
    MODE_HDR_KEY = 'HIERARCH ESO INS MODE'
    SKIP_KEYS['HIERARCH ESO DPR TYPE'] = ['DARK', 'OBJECT,SKY', 'OBJECT,FP',
                                          'OBJECT,DARK', 'FLUX,STD,SKY',
                                          'EFF,SKY,SKY']
    SKIP_KEYS['HIERARCH ESO DPR CATG'] = ['TEST', 'SCIENCE']
    REMOVE_SUFFIX = '_HA'
    # translate between dprtype and header
    TRANSLATE = dict()
    TRANSLATE['DARK_FLAT'] = {'HIERARCH ESO DPR TYPE': ['ORDERDEF,DARK,LAMP', 'FLAT,DARK,LAMP']}
    TRANSLATE['DARK_FP'] =  {'HIERARCH ESO DPR TYPE': ['CONTAM,DARK,FP']}
    TRANSLATE['FLAT_DARK'] =  {'HIERARCH ESO DPR TYPE': ['ORDERDEF,LAMP,DARK', 'FLAT,LAMP,DARK']}
    TRANSLATE['FLAT_FLAT'] =  {'HIERARCH ESO DPR TYPE': ['FLAT,LAMP,LAMP']}
    TRANSLATE['FP_DARK'] = {'HIERARCH ESO DPR TYPE': ['CONTAM,FP,DARK']}
    TRANSLATE['FP_FP'] =  {'HIERARCH ESO DPR TYPE': ['WAVE,FP,FP']}
    TRANSLATE['FP_HC'] =  {'HIERARCH ESO DPR TYPE': ['WAVE,FP,UN1']}
    TRANSLATE['HC_FP'] =  {'HIERARCH ESO DPR TYPE': ['WAVE,UN1,FP']}
    TRANSLATE['HC_HC'] =  {'HIERARCH ESO DPR TYPE': ['WAVE,UN1,UN1']}

    FORCE_KEYS = {'HIERARCH ESO DPR CATG': 'TEST',
                  'MTLUPDATE': 'True'}


elif INSTRUMENT == 'NIRPS' and MODE == 'HE':
    REF_DIR = '/nirps_raw/nirps/apero-data/drs-misc/ref_files/NIRPS_HE'
    TEST_DIR = '/nirps_raw/nirps/raw-data/nirps_he/'
    HEADER_KEYS = ['HIERARCH ESO DPR TYPE', 'HIERARCH ESO DPR CATG']
    MODE_HDR_KEY = 'HIERARCH ESO INS MODE'
    SKIP_KEYS['HIERARCH ESO DPR TYPE'] = ['DARK', 'OBJECT,SKY', 'OBJECT,FP',
                                          'OBJECT,DARK', 'FLUX,STD,SKY',
                                          'EFF,SKY,SKY']
    SKIP_KEYS['HIERARCH ESO DPR CATG'] = ['TEST', 'SCIENCE']
    REMOVE_SUFFIX = '_HE'
    # translate between dprtype and header
    TRANSLATE = dict()
    TRANSLATE['DARK_FLAT'] = {'HIERARCH ESO DPR TYPE': ['ORDERDEF,DARK,LAMP', 'FLAT,DARK,LAMP']}
    TRANSLATE['DARK_FP'] =  {'HIERARCH ESO DPR TYPE': ['CONTAM,DARK,FP']}
    TRANSLATE['FLAT_DARK'] =  {'HIERARCH ESO DPR TYPE': ['ORDERDEF,LAMP,DARK', 'FLAT,LAMP,DARK']}
    TRANSLATE['FLAT_FLAT'] =  {'HIERARCH ESO DPR TYPE': ['FLAT,LAMP,LAMP']}
    TRANSLATE['FP_DARK'] = {'HIERARCH ESO DPR TYPE': ['CONTAM,FP,DARK']}
    TRANSLATE['FP_FP'] =  {'HIERARCH ESO DPR TYPE': ['WAVE,FP,FP']}
    TRANSLATE['FP_HC'] =  {'HIERARCH ESO DPR TYPE': ['WAVE,FP,UN1']}
    TRANSLATE['HC_FP'] =  {'HIERARCH ESO DPR TYPE': ['WAVE,UN1,FP']}
    TRANSLATE['HC_HC'] =  {'HIERARCH ESO DPR TYPE': ['WAVE,UN1,UN1']}

    FORCE_KEYS = {'HIERARCH ESO DPR CATG': 'TEST',
                  'MTL_EDIT': 'True'}

else:
    raise ValueError(f'Unsupported instrument={INSTRUMENT} mode={MODE}')

# set the sub-directory to test
OBS_DIR = '2022-06-11'


# =============================================================================
# Define functions
# =============================================================================
class FileMatch:
    def __init__(self, keys, no_match=False):
        self.keys = keys
        self.no_match = no_match
        self.header_keys = dict()
        self.match_keys = dict()
        self.match_num = -1
        self.full_path = ''
        self.fraction = np.nan
        self.correct = False

    def __str__(self) -> str:
        """

        :return:
        """
        string = ''
        for key in self.keys:
            string += f'HDR[{key}]={self.header_keys[key]} '
        for key in self.keys:
            string += f'MATCH[{key}]={self.match_keys[key]} '
        string += f'FRAC={100 * self.fraction:.2f}%'
        return string

    def __repr__(self) -> str:
        return self.__str__()

    def info(self):
        string = ''
        for key in self.keys:
            string += f'\n\tHDR[{key}]={self.header_keys[key]} '
        for key in self.keys:
            string += f'\n\tMATCH[{key}]={self.match_keys[key]} '
        string += f'\n\tFRAC={100 * self.fraction:.2f}%'
        return string

    def populate(self, hdr):

        for key in self.keys:
            self.header_keys[key] = hdr[key]

            if self.no_match:
                self.match_keys[key] = None
            else:
                self.match_keys[key] = self.keys[key]

    def is_correct(self):
        for key in self.keys:
            if self.header_keys[key] in self.match_keys[key]:
                self.correct = True
            else:
                self.correct = False


def im2mask(image ,percentile=95):
    # get 95th percentile mask of a 512x512 sampling of the image
    # this gets 1 in 8x8 pixel and speeds computation. It also allows
    # to keep everything tidy in a dictionary without abusing the memory
    y, x = np.indices([512,512]) * 8 + 4
    image2 = image[y, x]
    return image2 > np.nanpercentile(image2, percentile)


def make_ref_cude():
    ref_files = glob.glob('{}/ref_*.fits'.format(REF_DIR))
    ref_files.sort()

    ref_cube = np.zeros([512, 512, len(ref_files)])
    types = np.array(ref_files)
    for ifile in tqdm(range(len(ref_files)), leave=False):
        ref_image = fits.getdata(ref_files[ifile])
        ref_cube[:, :, ifile] = im2mask(ref_image)
        basename = os.path.basename(ref_files[ifile])
        kind = basename.split('ref_')[1].split('.')[0]

        if REMOVE_SUFFIX is not None:
            if kind.endswith(REMOVE_SUFFIX):
                kind = kind[:-len(REMOVE_SUFFIX)]

        types[ifile] = kind

    return types, ref_cube


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    test_dir = os.path.join(TEST_DIR, OBS_DIR)
    files = glob.glob(test_dir + '/*.fits')

    print('Making reference cube from reference files')
    types, ref_cube = make_ref_cude()

    # storage of good files
    good_files = dict()
    # storage of bad files
    bad_files = dict()


    print(f'Analysing files for dir: {test_dir}')
    for filename in tqdm(files):
        # get basename (for dict key)
        basename = os.path.basename(filename)
        # get image and header
        image = fits.getdata(filename)
        hdr = fits.getheader(filename)
        # ---------------------------------------------------------------------
        # deal with wrong mode
        if MODE_HDR_KEY is not None:
            if hdr[MODE_HDR_KEY] != MODE:
                continue
        # ---------------------------------------------------------------------
        # deal with files we should not check
        no_check = False
        # loop around skip keys
        for skey in SKIP_KEYS:
            # loop around values of keys to check
            for value in SKIP_KEYS[skey]:
                # if value if
                if hdr[skey] == value:
                    no_check = True
                    break
            if no_check:
                break
        # skip file
        if no_check:
            continue

        # print header type
        if DEBUG:
            # print file
            print('\tFile: {0}'.format(filename))
            for key in HEADER_KEYS:
                print(f'\t\tHeader {key}: {hdr[key]}')
        # ---------------------------------------------------------------------
        # get brightest pixels from file of interest
        mask_from_file = im2mask(image)
        norm = np.nansum(mask_from_file)

        # set up the fraction
        frac_similar = np.zeros(len(types), dtype=float)

        # loop around reference files
        for ifile in range(len(types)):
            # calculate how similar the image is
            tmp = mask_from_file * ref_cube[:, :, ifile]
            frac_similar[ifile] = np.nansum(tmp) / norm
        # sort matches by highest similarity
        matches = np.argsort(-frac_similar)
        # ---------------------------------------------------------------------

        # print match
        if DEBUG:
            pargs = [types[matches[0]], 100 * frac_similar[matches[0]],
                     types[matches[1]], 100 * frac_similar[matches[1]]]
            pmsg = ('\t\t 1st match = {0} [{1:.2f}%]'
                    '\n\t\t 2nd match = {2}, [{3:.2f}%]')
            print(pmsg.format(*pargs))

        # ---------------------------------------------------------------------
        # compare best match to header
        # ---------------------------------------------------------------------
        match1 = types[matches[0]]

        if match1 in TRANSLATE:
            # set up file/hdr match 1
            file_match1 = FileMatch(TRANSLATE[match1])
            file_match1.populate(hdr)
            # test match 1
            file_match1.is_correct()
            file_match1.full_path = filename
            file_match1.match_num = 1
            file_match1.fraction = frac_similar[matches[0]]

            if file_match1.correct:
                good_files[basename] = file_match1
            else:
                bad_files[basename] = file_match1
            continue
        # ---------------------------------------------------------------------
        # if we have got to here file is bad and no match
        # ---------------------------------------------------------------------
        # else we have no match
        else:
            no_match = FileMatch(keys=HEADER_KEYS, no_match=True)
            no_match.populate(hdr)
            no_match.full_path = filename
            bad_files[basename] = no_match

    # -------------------------------------------------------------------------
    # Update raw files
    # -------------------------------------------------------------------------
    if UPDATE:

        for bad_file in bad_files:
            # get file instance
            file_match = bad_files[bad_file]
            # Get absolute filepath
            abs_path = os.path.join(test_dir, bad_file)
            # ask user to verify
            msg = f'\n\nFile: {file_match.full_path}'
            msg += file_match.info()
            msg += '\n\n\tAccept match? [Y]es or [N]o:\t'

            userinput = input(msg)

            if 'Y' in userinput.upper():

                # open fits file
                with fits.open(abs_path) as hdulist:
                    # loop around match keys
                    for key in file_match.match_keys:
                        # get match key list
                        match_keys = file_match.match_keys[key]

                        # if we have more than one entry we must choose
                        if len(match_keys) > 1:
                            msg = '\n\tChoose key:'
                            for m_it in range(len(match_keys)):
                                msg += f'\n\t\t{m_it + 1}.{match_keys[m_it]}'

                            match_key = match_keys[m_it - 1]
                        # else we know the match key
                        else:
                            match_key = match_keys[0]

                        hdulist[0].header[key] = match_key

                    # deal with force keys
                    for fkey in FORCE_KEYS:
                        hdulist[0].header[fkey] = FORCE_KEYS[fkey]

                    print('\n\t Overwriting file')
                    if not TEST:
                        hdulist.writeto(file_match.full_path, overwrite=True)



# =============================================================================
# End of code
# =============================================================================
