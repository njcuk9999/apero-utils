#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Check (and correct) the FP peak counting from order to order in the FP line
tables of a wave reference

The FP peak number at a fixed reference pixel varies smoothly from one order
to the next. An order where the peaks were miscounted stands out as an offset
by an integer number of peaks from a smooth polynomial across orders, and its
peak numbers are shifted back by that integer.

Created on 2026-10-07 at 10:53

@author: artigau
"""
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import glob
from typing import Tuple
from aperocore import math as mp


# =============================================================================
# Define functions
# =============================================================================
def table_adjust_fp_peak_numbers(tbl: Table, nsigcut: float = 3.0) -> Table:
    """
    Check the FP peak numbers of an FP line table from order to order and
    correct the orders where the peaks were miscounted

    The FP peak number at a fixed reference pixel (near the order center)
    varies smoothly from one order to the next. We fit a robust polynomial to
    this peak number as a function of order, and an order that is offset from
    the fit by an integer number of peaks has its peak numbers shifted back by
    that integer. Orders are corrected one at a time, starting with the worst
    one, until a correction no longer reduces the dispersion of the residuals

    :param tbl: astropy.table.Table, the FP line table, must have the columns
                'ORDER', 'PIXEL_MEAS' and 'PEAK_NUMBER'
    :param nsigcut: float, the threshold sigma above which an order is
                    considered an outlier in the robust polynomial fit
    :return: astropy.table.Table, the FP line table with the corrected
             'PEAK_NUMBER' values (the input table is also updated in place)
    """
    # reference pixel (near the order center) at which we evaluate the FP peak
    #   number of each order
    pixref = (np.max(tbl['PIXEL_MEAS'])- np.min(tbl['PIXEL_MEAS']))//2

    # ---------------------------------------------------------------------
    # storage for the (fractional) peak number at pixref for each order
    midcount = np.zeros(np.max(tbl['ORDER'])+1)
    # storage for the correction of the peak numbers for each order (an
    #   integer number of peaks)
    err = np.zeros(np.max(tbl['ORDER'])+1,dtype=int)
    # ---------------------------------------------------------------------
    # loop around orders
    for order_num in np.unique(tbl['ORDER']):
        # get the FP lines for this order
        tbl2 = tbl[tbl['ORDER'] == order_num]
        # get the measured pixel position of the lines
        pix = np.array(tbl2['PIXEL_MEAS'].data)
        # get the FP peak number of the lines
        count = np.array(tbl2['PEAK_NUMBER'].data)

        # only keep lines with a valid position and a valid peak number
        valid = ~np.isnan(pix) & ~np.isnan(count)
        pix = pix[valid]
        count = count[valid]

        # find the index of the last line before pixref (lines are sorted
        #   by increasing pixel position so the index in the masked vector
        #   is also the index in pix)
        i_prev = np.argmax(pix[pix < pixref])

        # fit a straight line to the peak number as a function of pixel
        #   position for the 5 lines around pixref
        fit = np.polyfit(pix[i_prev-2:i_prev+3], count[i_prev-2:i_prev+3], deg=1)

        # evaluate the straight line to find the count at pixel pixref
        count_pixref = np.polyval(fit, pixref)

        # store the count at pixref for this order
        midcount[order_num] = count_pixref

    # ---------------------------------------------------------------------
    # iteratively correct the orders that are offset by an integer number
    #   of peaks, one order at a time, starting with the worst one
    # ---------------------------------------------------------------------
    # flag to keep correcting orders
    redo = True
    # iteration counter (number of corrections applied)
    ite=0
    # loop until correcting the worst order no longer improves the fit
    while redo:
        # robust fit of a smooth polynomial to the count at pixref as a
        #   function of order number
        fit, _ = mp.robust_polyfit(np.arange(len(midcount)), (midcount),
                                   degree=11, nsigcut=3.0)

        # get the residuals to the fit (in units of FP peaks)
        res = (midcount) - np.polyval(fit, np.arange(len(midcount)))
        plt.plot(res, marker='o', label=f'Iteration {ite}')
        # find the order with the largest residual
        worst = np.argmax(np.abs(res))

        # copy the counts and shift the worst order by its residual rounded
        #   to the nearest integer number of peaks
        midcount_tmp = np.array(midcount)
        midcount_tmp[worst] -= np.round(res[worst])
        # get the residuals to the same fit once the worst order is shifted
        res2 = (midcount_tmp) - np.polyval(fit, np.arange(len(midcount_tmp)))
        # if the shift reduces the dispersion of the residuals we keep it
        if np.nanstd(res2) < np.nanstd(res):
            # update the counts
            midcount = midcount_tmp
            # add the shift to the correction for this order
            err[worst] -= np.round(res[worst])
            # increment the iteration counter
            ite+=1
        # otherwise there is nothing left to correct and we stop here
        else:
            redo = False

    # ---------------------------------------------------------------------
    # add the legend and show the residuals for all iterations
    plt.legend()
    plt.show()

    # ---------------------------------------------------------------------
    # apply the corrections (if any) to the peak numbers in the table
    if np.any(err != 0):
        # loop around orders
        for iord in np.unique(tbl['ORDER']):
            # only deal with orders that need a correction
            if err[iord] != 0:
                # print the corrected count and the correction applied
                # find the lines for this order
                g = tbl['ORDER'] == iord
                # shift the peak numbers of this order by the correction
                tbl['PEAK_NUMBER'][g] += err[iord]

    return tbl

# =============================================================================
# Start of code
# =============================================================================
# get the FP line tables of the wave reference (one per fiber)
files =glob.glob('2F3798BAE7a_pp_e2dsff_*_waveref_fplines_*.fits')

# loop around the FP line tables
for file in files:
    # load the FP line table
    print(file)
    tbl = Table.read(file)
    table_adjust_fp_peak_numbers(tbl, nsigcut=3.0)
