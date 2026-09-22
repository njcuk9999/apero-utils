#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2026-09-22

@author: cook
"""
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# =============================================================================
# Define variables
# =============================================================================
# Input wave
WAVEFILE1 = '/scratch2/spirou/misc/wave_spline/input/REF_WAVE_2400416c_C.fits'
# Output weave
WAVEFILE2 = '/scratch2/spirou/misc/wave_spline/output/REF_WAVE_2400416c_C.fits'

WAVESIZE1 = 4088
WAVESIZE2 = 8176

MODE = 'MAKE'
# MODE = 'TEST'

# =============================================================================
# Define functions
# =============================================================================
def make_wave_spline():
    """
    Make a wave spline file from a wave file

    :return: None
    """
    # load the wavefile1 data
    with fits.open(WAVEFILE1) as wavefile1:
        wave1 = wavefile1[0].data
        # get the shape of wave1 should be (number of orders, WAVESIZE1)
        shape1 = wave1.shape

        # check that the shape is correct
        if shape1[1] != WAVESIZE1:
            emsg = ('The input wave file {0} has shape {1} but expected shape '
                    '{2}')
            eargs = [WAVEFILE1, shape1, (WAVESIZE1,)]
            raise ValueError(emsg.format(*eargs))

        # create the wave2 array
        wave2 = np.zeros([shape1[0], WAVESIZE2])

        # set every other column of wave2 to be wave1 and the other
        # columns to be wave1 + gradient/2
        wave2[:, 0::2] = wave1
        wave2[:, 1::2] = wave1 + np.gradient(wave1, axis=1) / 2.0

        # Need to save wave2 exactly back into wave1 file but into WAVEFILE2
        # preserve all headers and other extensions

        # create a new HDUList
        hdulist = fits.HDUList()
        # copy the primary HDU
        primary_hdu = fits.PrimaryHDU(header=wavefile1[0].header,
                                      data=wave2)
        hdulist.append(primary_hdu)
        # write to WAVEFILE2
        hdulist.writeto(WAVEFILE2, overwrite=True)


def test_wave_spline():
    """
    Test the wave spline file

    :return: None
    """
    # load the wavefile1 data
    with fits.open(WAVEFILE1) as wavefile1:
        wave1 = wavefile1[0].data

    # load the wavefile2 data
    with fits.open(WAVEFILE2) as wavefile2:
        wave2 = wavefile2[0].data

    # make a fake data set of size (number of orders, WAVESIZE1)
    # sine function
    data1 = np.arange(0, np.prod(wave1.shape), 1).reshape(wave1.shape)
    # make a fake data set of using the same make_wave_spline logic
    # sine function (only diff with data1 is that it is double the size)
    data2 = np.arange(0, np.prod(wave1.shape), 0.5).reshape(wave2.shape)

    # plot wave1 and wave2
    plot_spec(wave1, data1, label='Wave1', colours=['orange', 'purple'],
              marker='x')
    plot_spec(wave2, data2, label='Wave2', colours=['red', 'green'],
              marker='+')
    plt.legend()
    plt.show()


def plot_spec(wave, spectrum, label, colours, marker):
    for order_num in range(spectrum.shape[0]):
        colour = colours[0] if order_num % 2 == 0 else colours[1]
        plt.plot(wave[order_num], spectrum[order_num], color=colour,
                 label=label if order_num in [0, 1] else None,
                 marker=marker)



# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    if MODE == 'MAKE':
        make_wave_spline()
    elif MODE == 'TEST':
        test_wave_spline()
    else:
        raise ValueError(f'Unknown MODE: {MODE}')




# =============================================================================
# End of code
# =============================================================================

