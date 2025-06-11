#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-06-10 at 12:31

@author: cook
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import constants
from astropy import units as u
from tqdm import tqdm
import warnings

# =============================================================================
# Define variables
# =============================================================================
BASENAME = 'NIRPS.2023-03-03T09:53:21.430'

ORDER_NUM = 60

speed_of_light_ms = constants.c.value

# Set global font size
mpl.rcParams.update({'font.size': 24})
# Set background color in plots
PLOT_BACKGROUND_COLOR = '#FEFDE1'
PLOT_BACKGROUND_COLOR2 = '#E1F0FE'


# =============================================================================
# Define classes
# =============================================================================
class File:
    def __init__(self):
        self.path = ''
        self.obsdir = None
        self.basename = ''
        self.prefix = ''
        self.suffix = ''
        self.format = None

    def get_abspath(self):
        basename = self.prefix + self.basename + self.suffix
        if self.obsdir is None:
            abspath = os.path.join(self.path,basename)
        else:
            abspath = os.path.join(self.path, self.obsdir, basename)
        if not os.path.exists(abspath):
            raise ValueError(f'File {abspath} does not exist.')
        return abspath

    def load_data(self, **kwargs):
        """Load the file data."""
        abspath = self.get_abspath()

        if self.format == 'image':
            return fits.getdata(abspath, **kwargs)
        elif self.format == 'table':
            return Table.read(abspath, **kwargs)

    def load_header(self, **kwargs):
        abspath = self.get_abspath()
        return fits.getheader(abspath, **kwargs)


# =============================================================================
# Define functions
# =============================================================================
def apero2eso_ord(order_num: int, fmt='apero') -> int:
    """Convert APERO order number to ESO order number."""

    # load csv
    order_def = Table.read('order_def_apero_eso.csv')

    if fmt == 'apero':
        order_dict = dict(zip(order_def['APERO_PYTHON'],
                              order_def['ESO_PYTHON']))
    elif fmt == 'eso':
        order_dict = dict(zip(order_def['ESO_PYTHON'],
                              order_def['APERO_PYTHON']))
    else:
        raise ValueError(f'Format {fmt} not recognized. Use "apero" or "eso".')

    if order_num not in order_dict:
        raise ValueError(f'Order number {order_num} not found in order definition.')

    if np.isnan(order_dict[order_num]):
        raise ValueError(f'Order number {order_num} is not defined in order definition.')

    return int(order_dict[order_num])


def get_echelle_ord(order_num: int, fmt='apero') -> int:
    # load csv
    order_def = Table.read('order_def_apero_eso.csv')

    if fmt == 'apero':
        order_dict = dict(zip(order_def['APERO_PYTHON'],
                              order_def['ECHELLE_ORDER']))
    elif fmt == 'eso':
        order_dict = dict(zip(order_def['ESO_PYTHON'],
                              order_def['ECHELLE_ORDER']))
    else:
        raise ValueError(f'Format {fmt} not recognized. Use "apero" or "eso".')

    if order_num not in order_dict:
        raise ValueError(f'Order number {order_num} not found in order definition.')

    if np.isnan(order_dict[order_num]):
        raise ValueError(f'Order number {order_num} is not defined in order definition.')

    return int(order_dict[order_num])


def doppler_shift(wavegrid: np.ndarray, velocity: float) -> np.ndarray:
    """
    Apply a doppler shift

    :param wavegrid: wave grid to shift
    :param velocity: float, velocity expressed in m/s

    :return: np.ndarray, the updated wave grid
    """
    # relativistic calculation (1 - v/c)
    part1 = 1 - (velocity / speed_of_light_ms)
    # relativistic calculation (1 + v/c)
    part2 = 1 + (velocity / speed_of_light_ms)
    # return updated wave grid
    return wavegrid * np.sqrt(part1 / part2)



def plot_spec_zoom(wavemap1, spectrum1, wavemap2, spectrum2,
                   zoom, color1, color2, label1, label2):
    plt.close()
    fig, frames = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))

    if color1 == 'blank':
        color1a, color1b = PLOT_BACKGROUND_COLOR, PLOT_BACKGROUND_COLOR2
    else:
        color1a, color1b = color1, color1
    if color2 == 'blank':
        color2a, color2b = PLOT_BACKGROUND_COLOR, PLOT_BACKGROUND_COLOR2
    else:
        color2a, color2b = color2, color2

    frames[0].plot(wavemap1, spectrum1, color=color1a, label=label1)
    frames[0].plot(wavemap2, spectrum2, color=color2a, label=label2)

    frames[1].plot(wavemap1, spectrum1, color=color1b, label=label1)
    frames[1].plot(wavemap2, spectrum2, color=color2b, label=label2)
    frames[1].set_facecolor('yellow')

    frames[1].set_xlim(*zoom)
    frames[0].set_yticklabels([])
    frames[1].set_yticklabels([])
    frames[0].set(xlabel='Wavelength (nm)', ylabel='flux')
    frames[1].set(xlabel='Wavelength (nm)', ylabel='flux')

    frames[0].set_facecolor(PLOT_BACKGROUND_COLOR)
    frames[1].set_facecolor(PLOT_BACKGROUND_COLOR2)

    frames[0].legend(loc='lower center')

    return frames


# =============================================================================
# Define files
# =============================================================================
apero = File()
apero.path = '/cosmos99/nirps/apero-data/nirps_he_online/objects/GL699'
apero.obsdir = None
apero.basename = BASENAME
apero.suffix = 'e.fits'
apero.format = 'image'

eso = File()
eso.path = '/data/cook/nirps_comp/geneva/GL699/'
eso.obsdir = '2023-03-02'
eso.basename = BASENAME
eso.prefix = 'r.'
eso.suffix = '_S2D_BLAZE_A.fits'
eso.format = 'image'


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # load data
    apero_sp = apero.load_data(extname='FluxA')
    apero_wave = apero.load_data(extname='WaveA')
    apero_hdr = apero.load_header(ext=1)

    eso_sp = eso.load_data(extname='SCIDATA')
    eso_wave = eso.load_data(extname='WAVEDATA_VAC_BARY') / 10.0
    eso_hdr = eso.load_header(ext=0)

    apero_ordnum = int(ORDER_NUM)
    eso_ordnum = apero2eso_ord(ORDER_NUM)
    ech_ordnum = get_echelle_ord(ORDER_NUM)
    # get berv
    apero_berv = apero_hdr['BERV']
    eso_berv = eso_hdr['HIERARCH ESO QC BERV']
    eso_wave_bc = doppler_shift(eso_wave, eso_berv * 1000)

    # -------------------------------------------------------------------------
    frames_blank = plot_spec_zoom(apero_wave[apero_ordnum],
                                  apero_sp[apero_ordnum],
                                  eso_wave[apero_ordnum],
                                  eso_sp[apero_ordnum],
                                  zoom=(1624, 1628),
                                  color1='blank', color2='blank',
                                  label1=f'APERO [{apero_ordnum}]',
                                  label2=f'ESO [{apero_ordnum}]')

    frames_blank[0].text(0.5, 0.5, 'Full order', ha='center', va='center',
                         transform=frames_blank[0].transAxes)
    frames_blank[1].text(0.5, 0.5, 'Zoom in', ha='center', va='center',
                         transform=frames_blank[1].transAxes)
    plt.suptitle('Blank')
    plt.show()

    # -------------------------------------------------------------------------
    frames_sameord = plot_spec_zoom(apero_wave[apero_ordnum],
                                    apero_sp[apero_ordnum],
                                    eso_wave[apero_ordnum],
                                    eso_sp[apero_ordnum],
                                    zoom=(1624, 1628),
                                    color1='orange', color2='blue',
                                    label1=f'APERO [{apero_ordnum}]',
                                    label2=f'ESO [{apero_ordnum}]')
    plt.suptitle('Extracted order "60"')
    plt.show()

    # -------------------------------------------------------------------------
    frames_ech = plot_spec_zoom(apero_wave[apero_ordnum],
                                apero_sp[apero_ordnum],
                                eso_wave[eso_ordnum],
                                eso_sp[eso_ordnum],
                                zoom=(1624, 1628),
                                color1='orange', color2='blue',
                                label1=f'APERO [{apero_ordnum}]',
                                label2=f'ESO [{eso_ordnum}]')
    plt.suptitle(f'Echelle order {ech_ordnum}')
    plt.show()
    # -------------------------------------------------------------------------
    frames_bshift = plot_spec_zoom(apero_wave[apero_ordnum],
                                   apero_sp[apero_ordnum],
                                   eso_wave_bc[eso_ordnum],
                                   eso_sp[eso_ordnum],
                                   zoom=(1625, 1626),
                                   color1='orange', color2='blue',
                                   label1=f'APERO [{apero_ordnum}]',
                                   label2=f'ESO [{eso_ordnum}]')
    plt.suptitle(f'Echelle order {ech_ordnum}, ESO BERV SHIFT={eso_berv:.2f} km/s')
    plt.show()

    # -------------------------------------------------------------------------



# =============================================================================
# End of code
# =============================================================================
