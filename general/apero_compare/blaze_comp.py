#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-08-15 at 11:45

@author: cook
"""
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os


# =============================================================================
# Define variables
# =============================================================================
BLAZE_FILE = 'calib/B7F5E274D4_pp_blaze_AB.fits'

NORDERS = 49

PLOT_ALL = True

PLOT_DIR = '/scratch3/lbl/misc/blazecomp/'

# -----------------------------------------------------------------------------
# define names
name1 = 'spirou@rali'
name2 = 'cook@jupiter'
name3 = 'cook@nb19'
name4 = 'spirou@maestria'
name5 = 'newworlds'
name6 = 'lam'
name7 = 'cfht'
# reference name
ref_name = name1
# choose which ones to run
names = [name1, name2, name3, name4, name5, name6, name7]
# define data directory
outpaths = dict()
outpaths[name1] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_rali/'
outpaths[name2] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_jupiter/'
outpaths[name3] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_nb19/'
outpaths[name4] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_maestria/'
outpaths[name5] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_newworld/'
outpaths[name6] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_lam/'
outpaths[name7] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_cfht/'
# add a color for each reduction (i.e. b, g, r, k, m, c, orange, purple)
COLORS = dict()
COLORS[name1] = 'b'
COLORS[name2] = 'r'
COLORS[name3] = 'g'
COLORS[name4] = 'orange'
COLORS[name5] = 'purple'
COLORS[name6] = 'k'
COLORS[name7] = 'm'



# =============================================================================
# Define functions
# =============================================================================
def blaze_plot(frames, datadiff, data, refname, name):

    frames[0].plot(data, label=name, color=COLORS[name])

    frames[1].plot(datadiff, label=f'{name} - {refname}',
                   color=COLORS[name])

    return frames


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # find the biggest diff
    diff = dict()
    flux = dict()
    frac = dict()
    maxdiff = dict()
    maxfracdiff = dict()
    # load ref
    blaze_ref = fits.getdata(os.path.join(outpaths[ref_name], BLAZE_FILE))
    # loop around orders and calculate biggest diff
    for order_num in range(NORDERS):
        # set up figure
        if PLOT_ALL:
            fig, frames = plt.subplots(ncols=1, nrows=2, figsize=(20, 10))
        else:
            frames = None
        # compare blaze files
        for name in names:
            # load blaze for file
            blaze_local = fits.getdata(os.path.join(outpaths[name], BLAZE_FILE))

            blaze_diff = blaze_local[order_num] - blaze_ref[order_num]
            # plot
            if PLOT_ALL:
                blaze_plot(frames, blaze_diff, blaze_local[order_num],
                           ref_name, name)

            argmax = np.nanargmax(abs(blaze_diff))

            if order_num in diff:
                diff[order_num].append(blaze_diff[argmax])
                flux[order_num].append(blaze_local[order_num][argmax])
                frac[order_num].append(blaze_diff[argmax] / blaze_local[order_num][argmax])
            else:
                diff[order_num] = [blaze_diff[argmax]]
                flux[order_num] = [blaze_local[order_num][argmax]]
                frac[order_num] = [blaze_diff[argmax] / blaze_local[order_num][argmax]]

        if PLOT_ALL:
            # add legend labels
            frames[0].legend(loc=0)
            frames[1].legend(loc=0)
            frames[0].set(xlabel='pixel', ylabel='blaze')
            frames[1].set(xlabel='pixel', ylabel='blaze diff')
            plt.suptitle(f'Order {order_num}')

            plt.savefig(os.path.join(PLOT_DIR, 'blazecomp_order_{0}.png'.format(order_num)))
            plt.close()

        # store the order with the biggest diff
        maxdiff[order_num] = np.nanmax(np.abs(diff[order_num]))
        maxfracdiff[order_num] = np.nanmax(np.abs(frac[order_num]))

    # find the order with the biggest diff
    maxdifforder = np.argmax(list(maxdiff.values()))
    print('Order with biggest diff = {0}'.format(maxdifforder))

    maxfracdifforder = np.argmax(list(maxfracdiff.values()))
    print('Order with biggest frac diff = {0}'.format(maxfracdifforder))

    for order_num in range(NORDERS):
        print(f'\n\nOrder {order_num}')
        for it, name in enumerate(names):
            print(f'\t{name:16s}   max diff = {diff[order_num][it]:.5f}   '
                  f'flux at max diff = {flux[order_num][it]:.5f}   '
                  f'frac diff = {frac[order_num][it]:.5e}')


# =============================================================================
# End of code
# =============================================================================
