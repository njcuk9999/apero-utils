#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-01-24 at 17:09

@author: cook
"""
import os

import matplotlib.pyplot as plt
from astropy.io import fits

# =============================================================================
# Define variables
# =============================================================================
# mode = 'HA'
mode = 'HE'

if mode == 'HA':
    path = '/nirps_raw/nirps/apero-data/drs-data/nirps_ha_202211/red/2022-12-02'
    filename = 'NIRPS_2022-12-03T03_32_41_149'
    wavefile = '/nirps_raw/nirps/apero-data/drs-data/nirps_he_202211/calib/NIRPS_2022-11-25T12_00_56_233_pp_e2dsff_A_wavesol_ref_A.fits'
else:
    path = '/nirps_raw/nirps/apero-data/drs-data/nirps_he_202211/red/2022-12-05'
    filename = 'NIRPS_2022-12-06T08_19_48_184'
    wavefile = '/nirps_raw/nirps/apero-data/drs-data/nirps_ha_202211/calib/NIRPS_2022-11-25T12_24_04_050_pp_e2dsff_A_wavesol_ref_A.fits'

pclean_ext = '_pp_tellu_pclean_A.fits'
e2ds_ext = '_pp_e2dsff_A.fits'
tcorr_ext = '_pp_e2dsff_tcorr_A.fits'
recon_ext = '_pp_e2dsff_recon_A.fits'


# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def function1():
    return 0


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":

    hdr = fits.getheader(os.path.join(path, filename + e2ds_ext))
    e2ds = fits.getdata(os.path.join(path, filename + e2ds_ext))
    tcorr = fits.getdata(os.path.join(path, filename + tcorr_ext))
    recon = fits.getdata(os.path.join(path, filename + recon_ext))
    skycorr_sci = fits.getdata(os.path.join(path, filename + pclean_ext),
                               extname='SKYCORR_SCI')
    finite_res = fits.getdata(os.path.join(path, filename + pclean_ext),
                              extname='FINITE_RES')
    wavemap = fits.getdata(wavefile)

    objname = hdr['DRSOBJN']

    plt.close()
    fig = plt.figure()
    frame1 = plt.subplot2grid((5, 1), (0, 0), rowspan=2)
    frame2 = plt.subplot2grid((5, 1), (2, 0), rowspan=1, sharex=frame1)
    frame3 = plt.subplot2grid((5, 1), (3, 0), rowspan=1, sharex=frame1)
    frame4 = plt.subplot2grid((5, 1), (4, 0), rowspan=1, sharex=frame1)

    for order_num in range(wavemap.shape[0]):

        if order_num == 0:
            label1a, label1b = 'e2dsff', 'tcorr'
            label2, label3, label4 = 'skycorr', 'recon', 'finite res'
        else:
            label1a, label1b = None, None
            label2, label3, label4 = None, None, None

        frame1.plot(wavemap[order_num], e2ds[order_num], label=label1a,
                    color='k', alpha=0.8, lw=0.75)
        frame1.plot(wavemap[order_num], tcorr[order_num], label=label1b,
                    color='r', alpha=0.8, lw=0.75)

        frame2.plot(wavemap[order_num], skycorr_sci[order_num], label=label2,
                    color='orange', lw=0.75)

        frame3.plot(wavemap[order_num], recon[order_num], label=label3,
                    color='b', lw=0.75)

        frame4.plot(wavemap[order_num], finite_res[order_num], label=label4,
                    color='purple', lw=0.75)

    frame1.legend(loc=0)
    frame2.legend(loc=0)
    frame3.legend(loc=0)
    frame4.legend(loc=0)

    plt.suptitle(f'{objname}     [{filename}]')
    plt.show()
    plt.close()

# =============================================================================
# End of code
# =============================================================================
