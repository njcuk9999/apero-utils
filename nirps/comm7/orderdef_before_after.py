#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-11-20 at 19:19

@author: cook
"""
import os
import warnings

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import LinearStretch
from astropy.visualization import imshow_norm, ZScaleInterval
from mpl_toolkits.axes_grid1 import make_axes_locatable

# =============================================================================
# Define variables
# =============================================================================
fibers = ['A', 'B']
modes = ['HA', 'HE']
# --------------------------------------------------------------------------------------
labels = ['center', 'left', 'right', 'red', 'blue']
boxsize = 5
function = np.nanmedian
func_name = 'Med'
# --------------------------------------------------------------------------------------
PLOT = False
PLOT_PATH = '/nirps_raw/nirps/misc/comm7'
# --------------------------------------------------------------------------------------
RAW_DIR_HE = '/nirps_raw/nirps/raw-data/nirps_he/'
RAW_DIR_HA = '/nirps_raw/nirps/raw-data/nirps_ha/'
DPRKEY = 'HIERARCH ESO DPR TYPE'
# --------------------------------------------------------------------------------------
# switch files key = (mode, fiber)
CASE = 5
old_files = dict()
new_files = dict()
xpos_all = dict()
ypos_all = dict()

# 2022-11-17 vs 2022-11-19
if CASE == 1:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_15_02_461.fits'
    new_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_32_10_890.fits'
    xpos_all[('HE', 'A')] = [2124, 2259, 2358, 3515, 235]
    ypos_all[('HE', 'A')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_16_42_792.fits'
    new_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_33_51_221.fits'
    xpos_all[('HE', 'B')] = [2145, 2280, 2380, 3535, 260]
    ypos_all[('HE', 'B')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_22_57_043.fits'
    new_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_28_44_649.fits'
    xpos_all[('HA', 'A')] = [2124, 2259, 2358, 3515, 235]
    ypos_all[('HA', 'A')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_24_37_374.fits'
    new_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_30_24_979.fits'
    xpos_all[('HA', 'B')] = [2145, 2280, 2380, 3535, 260]
    ypos_all[('HA', 'B')] = [1952, 4022, 85, 1977, 1900]
# 2022-11-19 vs 2022-11-21
elif CASE == 2:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_32_10_890.fits'
    new_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-21/NIRPS_2022-11-21T16_46_59_021.fits'
    xpos_all[('HE', 'A')] = [2124, 2259, 2358, 3515, 235]
    ypos_all[('HE', 'A')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_33_51_221.fits'
    new_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-21/NIRPS_2022-11-21T16_48_39_353.fits'
    xpos_all[('HE', 'B')] = [2145, 2280, 2380, 3535, 260]
    ypos_all[('HE', 'B')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_28_44_649.fits'
    new_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-21/NIRPS_2022-11-21T16_43_38_361.fits'
    xpos_all[('HA', 'A')] = [2124, 2259, 2358, 3515, 235]
    ypos_all[('HA', 'A')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_30_24_979.fits'
    new_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-21/NIRPS_2022-11-21T16_45_18_691.fits'
    xpos_all[('HA', 'B')] = [2145, 2280, 2380, 3535, 260]
    ypos_all[('HA', 'B')] = [1952, 4022, 85, 1977, 1900]
# 2022-11-19 vs 2022-11-21
elif CASE == 3:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_32_10_890.fits'
    new_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-21/NIRPS_2022-11-22T03_36_09_070.fits'
    xpos_all[('HE', 'A')] = [2124, 2259, 2358, 3515, 235]
    ypos_all[('HE', 'A')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_33_51_221.fits'
    new_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-21/NIRPS_2022-11-22T03_37_49_401.fits'
    xpos_all[('HE', 'B')] = [2145, 2280, 2380, 3535, 260]
    ypos_all[('HE', 'B')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_28_44_649.fits'
    new_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-21/NIRPS_2022-11-22T03_32_42_829.fits'
    xpos_all[('HA', 'A')] = [2124, 2259, 2358, 3515, 235]
    ypos_all[('HA', 'A')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_30_24_979.fits'
    new_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-21/NIRPS_2022-11-22T03_34_23_159.fits'
    xpos_all[('HA', 'B')] = [2145, 2280, 2380, 3535, 260]
    ypos_all[('HA', 'B')] = [1952, 4022, 85, 1977, 1900]
# 2022-11-17 vs 2022-11-21
elif CASE == 4:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_15_02_461.fits'
    new_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-21/NIRPS_2022-11-22T03_36_09_070.fits'
    xpos_all[('HE', 'A')] = [2124, 2259, 2358, 3515, 235]
    ypos_all[('HE', 'A')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_16_42_792.fits'
    new_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-21/NIRPS_2022-11-22T03_37_49_401.fits'
    xpos_all[('HE', 'B')] = [2145, 2280, 2380, 3535, 260]
    ypos_all[('HE', 'B')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_22_57_043.fits'
    new_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-21/NIRPS_2022-11-22T03_32_42_829.fits'
    xpos_all[('HA', 'A')] = [2124, 2259, 2358, 3515, 235]
    ypos_all[('HA', 'A')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_24_37_374.fits'
    new_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-21/NIRPS_2022-11-22T03_34_23_159.fits'
    xpos_all[('HA', 'B')] = [2145, 2280, 2380, 3535, 260]
    ypos_all[('HA', 'B')] = [1952, 4022, 85, 1977, 1900]
# 2022-11-21 vs 2022-11-22[am]
elif CASE == 5:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-21/NIRPS_2022-11-22T03_36_09_070.fits'
    new_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-22/NIRPS_2022-11-22T17_12_31_818.fits'
    xpos_all[('HE', 'A')] = [2124, 2259, 2358, 3515, 235]
    ypos_all[('HE', 'A')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-21/NIRPS_2022-11-22T03_37_49_401.fits'
    new_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-22/NIRPS_2022-11-22T17_14_12_148.fits'
    xpos_all[('HE', 'B')] = [2145, 2280, 2380, 3535, 260]
    ypos_all[('HE', 'B')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-21/NIRPS_2022-11-22T03_32_42_829.fits'
    new_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-22/NIRPS_2022-11-22T17_09_11_157.fits'
    xpos_all[('HA', 'A')] = [2124, 2259, 2358, 3515, 235]
    ypos_all[('HA', 'A')] = [1952, 4022, 85, 1977, 1900]
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-21/NIRPS_2022-11-22T03_34_23_159.fits'
    new_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-22/NIRPS_2022-11-22T17_10_51_488.fits'
    xpos_all[('HA', 'B')] = [2145, 2280, 2380, 3535, 260]
    ypos_all[('HA', 'B')] = [1952, 4022, 85, 1977, 1900]
else:
    raise ValueError(f'Case = {CASE} not supported')


# =============================================================================
# Define functions
# =============================================================================
def plot_image(image, title, side='bottom', pad=0.0):
    if side in ['right', 'left']:
        orientation = 'vertical'
    else:
        orientation = 'horizontal'
    cmap1 = matplotlib.cm.get_cmap('gist_heat')

    fig, frame = plt.subplots(ncols=1, nrows=1)
    im, _ = imshow_norm(image, frame, origin='lower', aspect='auto',
                        interval=ZScaleInterval(), stretch=LinearStretch(),
                        cmap=cmap1, interpolation='None', rasterized=True)
    divider = make_axes_locatable(frame)
    cax = divider.append_axes(side, '5%', pad=pad)
    cbar = fig.colorbar(im, cax=cax, orientation=orientation)
    cbar.ax.tick_params(labelsize=8)

    plt.suptitle(title)


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # storage of text
    summary_text = []

    for mode in modes:
        for fiber in fibers:
            old_file = old_files[(mode, fiber)]
            new_file = new_files[(mode, fiber)]
            xpos = xpos_all[(mode, fiber)]
            ypos = ypos_all[(mode, fiber)]

            # load old and new file
            old_data = fits.getdata(old_file)
            new_data = fits.getdata(new_file)
            old_hdr = fits.getheader(old_file)
            new_hdr = fits.getheader(new_file)

            xindices, yindices = np.indices(old_data.shape)

            print('\n')
            print("*" * 50)
            print(f'* NIRPS_{mode} Fiber {fiber}')
            print('*' * 50)
            print(f'b={boxsize}')
            print(f'old={old_hdr["ORIGFILE"]}\n\tEXPTIME={old_hdr["EXPTIME"]}  '
                  f'DPRTYPE={old_hdr[DPRKEY]}   DATE={old_hdr["DATE"]}')
            print(f'new={new_hdr["ORIGFILE"]}\n\tEXPTIME={new_hdr["EXPTIME"]}  '
                  f'DPRTYPE={new_hdr[DPRKEY]}   DATE={new_hdr["DATE"]}')

            for pos in range(len(xpos)):
                print(f'{labels[pos]}  \tx,y={xpos[pos]:4d},{ypos[pos]:4d}')

            print('\n')

            summary_stat = f'{mode} mode fiber {fiber}: '
            stats = []

            for pos in range(len(xpos)):

                # mask positions
                mask = xindices > xpos[pos] - boxsize
                mask &= xindices < xpos[pos] + boxsize
                mask &= yindices > ypos[pos] - boxsize
                mask &= yindices < ypos[pos] + boxsize

                # sum pixels in old
                old_value = function(old_data[mask])
                new_value = function(new_data[mask])

                # ratio old to new
                ratio = new_value / old_value

                # print result
                print(f'{func_name}[{labels[pos]}]'
                      f'\told={old_value:8.3f} new={new_value:8.3f} ratio={ratio:.5f}')

                # summary stats
                if ratio < 1:
                    rstat = (1 - ratio) * 100
                    stat = f' {rstat:.2f}% less flux at {labels[pos]}'
                else:
                    rstat = (ratio - 1) * 100
                    stat = f' {rstat:.2f}% more flux at {labels[pos]}'
                stats.append(stat)

            # full image 90% percentile
            old_90 = np.nanpercentile(old_data, 90)
            new_90 = np.nanpercentile(new_data, 90)

            ratio_90 = new_90 / old_90

            print(f'Percentile[90]',
                  f'\told={old_90:8.3f} new={new_90:8.3f} ratio={ratio_90:.5f}')

            # summary stats
            if ratio_90 < 1:
                rstat = (1 - ratio_90) * 100
                stat = f' {rstat:.2f}% less flux 90th Percentile'
            else:
                rstat = (ratio_90 - 1) * 100
                stat = f' {rstat:.2f}% more flux 90th Percentile'
            stats.append(stat)

            # add to the stats
            summary_stat += ', '.join(stats)

            summary_text.append(summary_stat)

            # diff image
            with warnings.catch_warnings(record=True) as _:
                ratio_full = (new_data / old_data)

            if PLOT:
                plot_image(ratio_full, title=f'NIRPS_{mode} Fiber {fiber}',
                           side='right')
                plot_path = os.path.join(PLOT_PATH, f'FW_before_after_NIRPS_{mode}_fiber_{fiber}.pdf')
                plt.savefig(plot_path)
                plt.show()

    # print the summary
    print('\n\n Summary: \n\n')
    for stext in summary_text:
        print(stext + '\n')

# =============================================================================
# End of code
# =============================================================================
