"""
Code to tell us the maximum flux in a given set of files
and related this to the ND

auther: cook
date: 2022-11-23
"""
import glob
import os
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

# ==================================================================
# VARIABLES
# ==================================================================
# Define the time of the first file to use (should be slightly earlier)
#   Set to None for no limit
#   In fits format YYYY-MM-DDThh:mm:ss.ss
TIME_START = '2022-11-24T01:10:00.00'
# Define the time of the last file to use (should be slightly later)
#   Set to None for no limit
#   In fits format YYYY-MM-DDThh:mm:ss.ss
TIME_END = None
# Define the directory to find the files
PATH = '/data/NIRPS/INS_ROOT/SYSTEM/DETDATA/'
# Define the file string to select files (including wildcards)
FILESTR = 'NIRPS_*CONTAM*.fits'
# The percentiles to take (between 0 and 100)
PERCENTILES = [99, 95, 90, 50]
# colors for each percentile
PCOLORS = ['r', 'g', 'purple', 'b']
# ------------------------------------------------------------------
# date header key (fits format)
KW_DATE = 'DATE'
# mode header key
KW_MODE = 'HIERARCH ESO INS MODE'
# dprtype header key
KW_DPRTYPE = 'HIERARCH ESO DPR TYPE'
# exposure time header key
KW_EXPTIME = 'EXPTIME'
# ND filter 1 header key
KW_ND_FILT1 = 'HIERARCH ESO INS FILT1 ID'
# ND filter 2 header key
KW_ND_FILT2 = 'HIERARCH ESO INS FILT2 ID'
# ND samples header key
KW_ND_SAMPLES = 'HIERARCH ESO DET NDSAMPLES'
# ------------------------------------------------------------------
# Valid fiber A flux types
FIBER_A_DPRTYPE = ['CONTAM,TUN,DARK', 'CONTAM,FP,DARK', 'CONTAM,UN1,DARK']
# Valid fiber B flux types
FIBER_B_DPRTYPE = ['CONTAM,DARK,TUN', 'CONTAM,DARK,FP', 'CONTAM,DARK,UN1']


# ==================================================================
# CLASSES
# ==================================================================
class NIRPSException(Exception):
    pass


class Exposure:
    def __init__(self, filename):
        # the filename
        self.filename = filename
        self.basename = os.path.basename(self.filename)
        # the header
        self.header = None
        # the data
        self.data = None
        # the mode
        self.mode = None
        # the fiber
        self.fiber = None
        # the exposure time
        self.exptime = None
        # the ND filter value (for this fiber)
        self.nd = None
        # the number of samples
        self.nsamples = None
        # the number of saturated pixels
        self.n_saturated = 0
        # the fraction of saturated pixels
        self.sat_frac = 0
        # record the percentile values in ADU/s
        self.pvalues = np.zeros(len(PERCENTILES))
        # record the percentiles values in ADU
        self.pflux = np.zeros(len(PERCENTILES))

    def __str__(self):
        return 'EXP[{}]'.format(self.basename)

    def __repr__(self):
        return self.__str__()

    def from_header(self):
        # load file header
        header = fits.getheader(self.filename)
        try:
            # reject files outside time limits (if set)
            if time_start is not None:
                if Time(header[KW_DATE], format='fits') < time_start:
                    emsg = '\tSkipping: File taken too early: ' + header[KW_DATE]
                    print(emsg)
                    raise NIRPSException(emsg)
            if time_end is not None:
                if Time(header[KW_DATE], format='fits') > time_end:
                    emsg = '\tSkipping: File taken too late: ' + header[KW_DATE]
                    print(emsg)
                    raise NIRPSException(emsg)
            # get mode
            self.mode = header[KW_MODE].strip()
            # get dprtype
            self.dprtype = header[KW_DPRTYPE].strip()
            # get exposure time
            self.exptime = header[KW_EXPTIME]
            # get number of samples
            self.nsamples = header[KW_ND_SAMPLES]
            # get the dprtype
            if self.dprtype in FIBER_A_DPRTYPE:
                self.fiber = 'A'
                self.nd = self._nd_from_raw(header[KW_ND_FILT1])
            elif self.dprtype in FIBER_B_DPRTYPE:
                self.fiber = 'B'
                self.nd = self._nd_from_raw(header[KW_ND_FILT2])
            else:
                emsg = ('\tSkipping: DPRTYPE not in fiber A or B '
                        'list: {}'.format(self.dprtype))
                print(emsg)
                raise NIRPSException(emsg)
        except Exception as e:
            raise NIRPSException(e)

    def _nd_from_raw(self, rawvalue):
        try:
            return float(rawvalue[2:])
        except Exception as e:
            raise NIRPSException(e)

    def from_data(self):
        # make sure we always close the file
        with fits.open(self.filename) as hdulist:
            # image data
            image = np.array(hdulist[1].data)
            # number of reads data
            nreads = np.array(hdulist[3].data)
            # create a mask where nreads = nsamples
            mask = nreads == self.nsamples
            # record the number of saturated pixels
            self.n_saturated = np.sum(~mask)
            # fraction of saturated pixels
            self.sat_frac = self.n_saturated / np.product(image.shape)
            # find the percentiles in ADU
            self.pvalues = np.nanpercentile(image, PERCENTILES)
            # convert to ADU/s
            self.pflux = self.pvalues * self.exptime

    def add_to_dict(self, sdict):
        sdict = add_append(sdict, 'FILE', self.basename)
        sdict = add_append(sdict, 'MODE', self.mode)
        sdict = add_append(sdict, 'FIBER', self.fiber)
        sdict = add_append(sdict, 'DPRTYPE', self.dprtype)
        sdict = add_append(sdict, 'EXPTIME', self.exptime)
        sdict = add_append(sdict, 'ND', self.nd)
        sdict = add_append(sdict, 'SAT_FRAC', self.sat_frac)
        for percentile in PERCENTILES:
            sdict = add_append(sdict, 'P[{}]'.format(percentile), self.sat_frac)

        return sdict


# ==================================================================
# FUNCTIONS
# ==================================================================
def add_append(sdict, key, value):
    if key in sdict:
        sdict[key].append(value)
    else:
        sdict[key] = [value]
    return sdict


def plot_grid(storage):
    # close plots
    plt.close()
    # loop around each dprtype
    for key in storage:

        # get nd value
        nd = np.array(list(map(lambda x: getattr(x, 'nd'),
                               storage[key])))
        # get flux values
        flux = np.array(list(map(lambda x: getattr(x, 'pflux'),
                                 storage[key])))
        # get the fraction of saturated pixels
        fsat = np.array(list(map(lambda x: getattr(x, 'sat_frac'),
                                 storage[key])))
        # get the exposure time
        exptime = np.array(list(map(lambda x: getattr(x, 'exptime'),
                                    storage[key])))
        # set up plots
        fig, frames = plt.subplots(ncols=2, nrows=2, figsize=(20, 20))
        # sort by nd
        sortmask = np.argsort(nd)
        nd = nd[sortmask]
        flux = flux[sortmask]
        fsat = fsat[sortmask]
        exptime = exptime[sortmask]
        # -------------------------------------------------------------
        # loop around percentiles
        for it, percentile in enumerate(PERCENTILES):
            pflux = flux[:, it]
            # plot percentile values against nd
            frames[0][0].scatter(nd, pflux,
                                 label='Percentile[{}]'.format(percentile),
                                 facecolor='None', edgecolor=PCOLORS[it],
                                 marker='o')
            # plot percentile values against exposure time              
            frames[1][0].scatter(exptime, pflux,
                                 label='Percentile[{}]'.format(percentile),
                                 facecolor='None', edgecolor=PCOLORS[it],
                                 marker='o')

        # plot fraction of saturated pixels
        frames[0][1].plot(nd, fsat, marker='o', ls='None')
        frames[1][1].plot(exptime, fsat, marker='o', ls='None')
        # -------------------------------------------------------------        
        # legend and labels
        frames[0][0].legend(loc=0)
        frames[1][0].legend(loc=0)
        frames[0][0].set(xlabel='ND', ylabel='Flux / ADU')
        frames[1][0].set(xlabel='EXPTIME [s]', ylabel='Flux / ADU')
        frames[0][1].set(xlabel='ND', ylabel='Fraction of saturated pixels')
        frames[1][1].set(xlabel='EXPTIME',
                         ylabel='Fraction of saturated pixels')
        plt.subplots_adjust(right=0.99, left=0.1, wspace=0.4)
        plt.suptitle(key)
        plt.savefig('Flux-tweak-{0}.jpg'.format(key))
        plt.show()
        plt.close()


# ==================================================================
# START OF CODE
# ==================================================================
if __name__ == '__main__':

    # get files
    files = glob.glob(os.path.join(PATH, FILESTR))
    # ------------------------------------------------------------------
    # set up astropy time for start and end
    if TIME_START is not None:
        time_start = Time(TIME_START, format='fits')
    else:
        time_start = None
    if TIME_END is not None:
        time_end = Time(TIME_END, format='fits')
    else:
        time_end = None
    # ------------------------------------------------------------------
    # storage key = DPRTYPE
    storage = OrderedDict()
    # ------------------------------------------------------------------
    # loop around files
    for filename in files:
        # print progress
        print('Analysising file {}'.format(filename))
        # --------------------------------------------------------------
        # create exposure class
        exposure = Exposure(filename)
        # populate exposure using header
        try:
            exposure.from_header()
        except NIRPSException:
            continue
        # populate the number of saturated and percentile valeus
        exposure.from_data()
        # --------------------------------------------------------------
        # get name for storage
        name = '{0}-{1} {2}'.format(exposure.mode,
                                    exposure.fiber, exposure.dprtype)
        # save to storage
        if name not in storage:
            storage[name] = [exposure]
        else:
            storage[name].append(exposure)
    # ------------------------------------------------------------------
    # plot
    print('Plotting grids')
    plot_grid(storage)
    # ------------------------------------------------------------------
    # get data and put it in dictionary
    storage_dict = OrderedDict()
    for key in storage:
        for exposure in storage[key]:
            storage_dict = exposure.add_to_dict(storage_dict)
    # ------------------------------------------------------------------
    # save to table
    print('Writing results table')
    table = Table(storage_dict)
    table.write('Flux-tweak.txt', format='ascii.fixed_width', overwrite=True)
    table.write('Flux-tweak.fits', format='fits', overwrite=True)

# ==================================================================
# END OF CODE
# ==================================================================
