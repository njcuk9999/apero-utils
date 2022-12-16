#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-10-27 at 10:50

@author: cook
"""
import os
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import LinearStretch
from astropy.visualization import imshow_norm, ZScaleInterval
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
NAME1 = 'spirou@rali'
NAME2 = 'cook@jupiter'
NAME3 = 'cook@nb19'
# Define which reduction is the reference reduction
REF_NAME = str(NAME1)
# -----------------------------------------------------------------------------
# just add another entry here
#  i.e. paths[NAME3] = path/to/reduced/dir
paths = dict()
paths[NAME1] = '/scratch2/spirou/drs-data/minidata2_2022-12-15_rali/red'
paths[NAME2] = '/scratch2/spirou/drs-data/minidata2_07XXX/reduced'
paths[NAME3] = '/scratch2/spirou/drs-data/minidata2_2022-12-15_nb19/red'
# -----------------------------------------------------------------------------
E2DS_SUFFICES = ['_AB.fits', '_A.fits', '_B.fits', '_C.fits']
NORDERS = 49


# =============================================================================
# Define classes
# =============================================================================
class Comparison:
    def __init__(self, name1, name2, filename1, filename2):
        self.name1 = name1
        self.name2 = name2
        self.filename1 = filename1
        self.filename2 = filename2
        self.basename = os.path.basename(filename1)
        self.loaded = False
        self.kind1 = 'image'
        self.kind2 = 'image'
        # two types of possible files pp or e2ds
        self.image_kind1 = id_image_kind(self.filename1)
        self.image_kind2 = id_image_kind(self.filename2)

        self.image1 = None
        self.table1 = None
        self.image2 = None
        self.table2 = None
        self.hdr1 = None
        self.hdr2 = None
        self.exact = False
        self.stats = dict()

    def __str__(self) -> str:
        return f'Comparison[{self.basename}]'

    def __repr__(self) -> str:
        return self.__str__()

    def load(self):
        # if we have already loaded the data don't reload
        if self.loaded:
            return
        # load data 1
        lout1 = load_data(self.filename1)
        self.kind1, self.image1, self.table1, self.hdr1 = lout1
        # load data 2
        lout2 = load_data(self.filename2)
        self.kind2, self.image2, self.table2, self.hdr2 = lout2
        # check files are same type
        if self.kind1 != self.kind2:
            emsg = (f'File1 and File2 are not the same format\n'
                    f'\tFile1[{self.kind1}]: {self.filename1}\n'
                    f'\tFile2[{self.kind2}]: {self.filename2}\n')
            raise ValueError(emsg)
        # mark data as loaded
        self.loaded = True

    def unload(self):
        # if we have already unloaded the data don't unload again
        if not self.loaded:
            return
        self.loaded = False
        self.image1 = None
        self.table1 = None
        self.image2 = None
        self.table2 = None
        self.hdr1 = None
        self.hdr2 = None

    def comparison(self):
        self.load()
        # get exact match
        if self.kind1 == 'image':
            # if its an e2ds we use full image
            if self.image_kind1 == 'e2ds' and self.image1.shape[0] == NORDERS:
                orders = np.arange(self.image1.shape[0])
                order_str = '[{0}]'
                odesc = 'E2DS'
            # else we do all y-axis values at once
            else:
                orders = [np.arange(self.image1.shape[0])]
                order_str = ''
                odesc = 'PP'
            # loop around "orders" (or full image)
            for order_num in tqdm(orders, leave=False, desc=odesc):
                ostr = order_str.format(order_num)

                self.stats[f'STD1{ostr}'] = np.nanstd(self.image1)
                self.stats[f'STD2{ostr}'] = np.nanstd(self.image2)
                self.stats[f'MED1{ostr}'] = np.nanmedian(self.image1)
                self.stats[f'MED2{ostr}'] = np.nanmedian(self.image2)

                diff = self.image1 - self.image2
                if np.nansum(diff) == 0:
                    self.exact = True
                    self.stats[f'DIFF_STD{ostr}'] = 0
                    self.stats[f'DIFF_MEAN{ostr}'] = 0
                else:
                    self.exact = False
                    self.stats[f'DIFF_STD{ostr}'] = np.nanstd(diff)
                    self.stats[f'DIFF_MEAN{ostr}'] = np.nanmedian(diff)

        else:
            self.exact = True
            odesc = 'TABLE'
            for col in tqdm(self.table1.colnames, leave=False, desc=odesc):
                diff = self.table1[col] - self.table2[col]
                if np.nansum(diff) != 0:
                    self.exact = False
                    return
        # unload data
        self.unload()

    def compare_image(self):
        self.load()

        if self.kind1 != 'image':
            return

        ratio = self.image1 / self.image2
        diff = self.image1 - self.image2

        fig, frames = plt.subplots(ncols=2, nrows=1)

        _, _ = imshow_norm(ratio, frames[0], origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           interpolation='None', rasterized=True)
        _, _ = imshow_norm(diff, frames[1], origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           interpolation='None', rasterized=True)

        frames[0].set(title='Ratio')
        frames[1].set(title='Diff')

        plt.show()
        plt.close()
        self.unload()


# =============================================================================
# Define functions
# =============================================================================
def get_files(path: str) -> List[str]:
    fits_files = []
    for root, dir, files in os.walk(path):
        for filename in files:
            if filename.endswith('.fits'):
                fits_files.append(os.path.join(root, filename))
    return fits_files


def load_data(filename):
    with fits.open(filename) as ffile:
        if len(ffile) > 1:
            dataext = 1
        else:
            dataext = 0
        kind = 'image' if ffile[dataext].is_image else 'table'
        if kind == 'image':
            image = np.array(ffile[dataext].data)
            table = None
        else:
            image = None
            table = Table(ffile[dataext].data)
        hdr = fits.Header(ffile[0].header)

    return kind, image, table, hdr


def id_image_kind(filename: str):
    for suffix in E2DS_SUFFICES:
        if filename.endswith(suffix):
            return 'e2ds'

    return 'pp'


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":

    # storage
    exact_files = dict()
    non_exact_files = dict()
    failed_files = dict()
    # -------------------------------------------------------------------------
    for name in paths:
        # skip ref path
        if name == REF_NAME:
            continue
        # set names
        name1 = REF_NAME
        name2 = name
        # get paths
        path1 = paths[name1]
        path2 = paths[name2]
        # get all fits files
        files1 = get_files(path1)
        files2 = get_files(path2)
        # ---------------------------------------------------------------------
        # storage for linked files
        linked_files: List[Comparison] = []
        # get link files
        print('Linking all files')
        for filename1 in tqdm(files1):
            for filename2 in files2:
                if os.path.basename(filename1) == os.path.basename(filename2):
                    linked_files.append(Comparison(name1, name2, filename1,
                                                   filename2))

        # ---------------------------------------------------------------------
        # find exact matches
        exact_matches = []
        print('Finding exact matches')
        for it in tqdm(range(len(linked_files))):
            try:
                linked_file = linked_files[it]
                linked_file.comparison()
                if linked_file.exact:
                    exact_matches.append(it)
            except Exception as _:
                if name in failed_files:
                    failed_files[name].append(linked_files[it])
                else:
                    failed_files[name] = [linked_files[it]]
        msg = '\t Found {0} exact matches. {1} files did not match exactly'
        margs = [len(exact_matches), len(linked_files) - len(exact_matches)]
        print(msg.format(*margs))
        # ---------------------------------------------------------------------
        # only keep those that don't match
        exact_files[name] = np.array(linked_files)[np.array(exact_matches)]
        non_exact_files[name] = np.array(linked_files)[~np.in1d(np.arange(len(linked_files)), exact_files)]









# =============================================================================
# End of code
# =============================================================================
