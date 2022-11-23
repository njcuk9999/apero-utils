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
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

# =============================================================================
# Define variables
# =============================================================================
fibers = ['A', 'B']
modes = ['HA', 'HE']
# --------------------------------------------------------------------------------------
boxsize = 20
percentile = 95
func_name = 'P95'

# DATATYPE 'RAW' or 'PP'
DATATYPE = 'PP'
# --------------------------------------------------------------------------------------
PLOT = True
PLOT_PATH = '/nirps_raw/nirps/misc/comm7'
# --------------------------------------------------------------------------------------
if DATATYPE == 'RAW':
    DIR_HE = '/nirps_raw/nirps/raw-data/nirps_he/'
    DIR_HA = '/nirps_raw/nirps/raw-data/nirps_ha/'
elif DATATYPE == 'PP':
    DIR_HE = '/nirps_raw/nirps/apero-data/drs-data/nirps_he_202211/tmp/'
    DIR_HA = '/nirps_raw/nirps/apero-data/drs-data/nirps_ha_202211/tmp/'
else:
    raise ValueError(f'DATATYPE = {DATATYPE} not valid')
DPRKEY = 'HIERARCH ESO DPR TYPE'
# --------------------------------------------------------------------------------------
# switch files key = (mode, fiber)
CASE = 8
old_files = dict()
new_files = dict()

edge = 100
size = 500
# iamge shape (y, x)
image_shape = [4096, 4096]

if DATATYPE == 'RAW':
    # center raw image
    center_min = [image_shape[0] // 2 - size, image_shape[1] // 2 - size]
    center_max = [image_shape[0] // 2 + size, image_shape[1] // 2 + size]
    # order left raw image
    left_min = [image_shape[0] - size - edge, edge]
    left_max = [image_shape[0] - edge, image_shape[1] - edge]
    # order right raw image
    right_min = [edge, edge]
    right_max = [size - edge, image_shape[1] - edge]
    # blue det raw image
    blue_min = [edge, edge]
    blue_max = [image_shape[0] - edge, size + edge]
    # red det raw image
    red_min = [edge, image_shape[1] - size - edge]
    red_max = [image_shape[0] - edge, image_shape[1] - edge]
else:
    # center raw image
    center_min = [image_shape[0] // 2 - size, image_shape[1] // 2 - size]
    center_max = [image_shape[0] // 2 + size, image_shape[1] // 2 + size]
    # order left raw image
    left_min = [edge, image_shape[1] - size - edge]
    left_max = [image_shape[0] - edge, image_shape[1] - edge]
    # order right raw image
    right_min = [edge, edge]
    right_max = [image_shape[0] - edge, size + edge]
    # blue det raw image
    blue_min = [image_shape[0] - size - edge, edge]
    blue_max = [image_shape[0] - edge, image_shape[1] - edge]
    # red det raw image
    red_min = [edge, edge]
    red_max = [size + edge, image_shape[1] - edge]

labels = ['center', 'left', 'right', 'red', 'blue']
min_points = [center_min, left_min, right_min, blue_min, red_min]
max_points = [center_max, left_max, right_max, blue_max, red_max]

# CASE 1: 2022-11-17 vs 2022-11-19
if CASE == 1:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_15_02_461.fits'
    new_files[('HE', 'A')] = DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_32_10_890.fits'
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_16_42_792.fits'
    new_files[('HE', 'B')] = DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_33_51_221.fits'
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_22_57_043.fits'
    new_files[('HA', 'A')] = DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_28_44_649.fits'
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_24_37_374.fits'
    new_files[('HA', 'B')] = DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_30_24_979.fits'
# CASE 2: 2022-11-19 vs 2022-11-21
elif CASE == 2:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_32_10_890.fits'
    new_files[('HE', 'A')] = DIR_HE + '2022-11-21/NIRPS_2022-11-21T16_46_59_021.fits'
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_33_51_221.fits'
    new_files[('HE', 'B')] = DIR_HE + '2022-11-21/NIRPS_2022-11-21T16_48_39_353.fits'
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_28_44_649.fits'
    new_files[('HA', 'A')] = DIR_HA + '2022-11-21/NIRPS_2022-11-21T16_43_38_361.fits'
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_30_24_979.fits'
    new_files[('HA', 'B')] = DIR_HA + '2022-11-21/NIRPS_2022-11-21T16_45_18_691.fits'
# CASE 3: 2022-11-19 vs 2022-11-21
elif CASE == 3:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_32_10_890.fits'
    new_files[('HE', 'A')] = DIR_HE + '2022-11-21/NIRPS_2022-11-22T03_36_09_070.fits'
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_33_51_221.fits'
    new_files[('HE', 'B')] = DIR_HE + '2022-11-21/NIRPS_2022-11-22T03_37_49_401.fits'
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_28_44_649.fits'
    new_files[('HA', 'A')] = DIR_HA + '2022-11-21/NIRPS_2022-11-22T03_32_42_829.fits'
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_30_24_979.fits'
    new_files[('HA', 'B')] = DIR_HA + '2022-11-21/NIRPS_2022-11-22T03_34_23_159.fits'
# CASE 4: 2022-11-17 vs 2022-11-21
elif CASE == 4:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_15_02_461.fits'
    new_files[('HE', 'A')] = DIR_HE + '2022-11-21/NIRPS_2022-11-22T03_36_09_070.fits'
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_16_42_792.fits'
    new_files[('HE', 'B')] = DIR_HE + '2022-11-21/NIRPS_2022-11-22T03_37_49_401.fits'
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_22_57_043.fits'
    new_files[('HA', 'A')] = DIR_HA + '2022-11-21/NIRPS_2022-11-22T03_32_42_829.fits'
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_24_37_374.fits'
    new_files[('HA', 'B')] = DIR_HA + '2022-11-21/NIRPS_2022-11-22T03_34_23_159.fits'
# CASE 5: 2022-11-21 vs 2022-11-22[am]
elif CASE == 5:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = DIR_HE + '2022-11-21/NIRPS_2022-11-22T03_36_09_070.fits'
    new_files[('HE', 'A')] = DIR_HE + '2022-11-22/NIRPS_2022-11-22T17_12_31_818.fits'
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = DIR_HE + '2022-11-21/NIRPS_2022-11-22T03_37_49_401.fits'
    new_files[('HE', 'B')] = DIR_HE + '2022-11-22/NIRPS_2022-11-22T17_14_12_148.fits'
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = DIR_HA + '2022-11-21/NIRPS_2022-11-22T03_32_42_829.fits'
    new_files[('HA', 'A')] = DIR_HA + '2022-11-22/NIRPS_2022-11-22T17_09_11_157.fits'
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = DIR_HA + '2022-11-21/NIRPS_2022-11-22T03_34_23_159.fits'
    new_files[('HA', 'B')] = DIR_HA + '2022-11-22/NIRPS_2022-11-22T17_10_51_488.fits'
# CASE 6: 2022-11-17 vs 2022-11-22[am]
elif CASE == 6:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_15_02_461.fits'
    new_files[('HE', 'A')] = DIR_HE + '2022-11-22/NIRPS_2022-11-22T17_12_31_818.fits'
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_16_42_792.fits'
    new_files[('HE', 'B')] = DIR_HE + '2022-11-22/NIRPS_2022-11-22T17_14_12_148.fits'
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_22_57_043.fits'
    new_files[('HA', 'A')] = DIR_HA + '2022-11-22/NIRPS_2022-11-22T17_09_11_157.fits'
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_24_37_374.fits'
    new_files[('HA', 'B')] = DIR_HA + '2022-11-22/NIRPS_2022-11-22T17_10_51_488.fits'
# CASE 7: 2022-11-17 vs 2022-11-22[17:35]
elif CASE == 7:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_15_02_461.fits'
    new_files[('HE', 'A')] = DIR_HE + '2022-11-22/NIRPS_2022-11-22T20_37_39_024.fits'
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_16_42_792.fits'
    new_files[('HE', 'B')] = DIR_HE + '2022-11-22/NIRPS_2022-11-22T20_39_19_354.fits'
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_22_57_043.fits'
    new_files[('HA', 'A')] = DIR_HA + '2022-11-22/NIRPS_2022-11-22T20_34_18_362.fits'
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_24_37_374.fits'
    new_files[('HA', 'B')] = DIR_HA + '2022-11-22/NIRPS_2022-11-22T20_35_58_692.fits'
# CASE 8: 2022-11-17 vs 2022-11-22[18:29]
elif CASE == 8:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_15_02_461.fits'
    new_files[('HE', 'A')] = DIR_HE + '2022-11-22/NIRPS_2022-11-22T21_33_45_669.fits'
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_16_42_792.fits'
    new_files[('HE', 'B')] = DIR_HE + '2022-11-22/NIRPS_2022-11-22T21_35_25_999.fits'
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_22_57_043.fits'
    new_files[('HA', 'A')] = DIR_HA + '2022-11-22/NIRPS_2022-11-22T21_30_25_007.fits'
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_24_37_374.fits'
    new_files[('HA', 'B')] = DIR_HA + '2022-11-22/NIRPS_2022-11-22T21_32_05_338.fits'
else:
    raise ValueError(f'Case = {CASE} not supported')


# =============================================================================
# Define functions
# =============================================================================
def plot_ratio(image1, image2, ratio, title, side='bottom', pad=0.0,
               minpoints=None, maxpoints=None, labels=None):
    if side in ['right', 'left']:
        orientation = 'vertical'
    else:
        orientation = 'horizontal'
    cmap1 = matplotlib.cm.get_cmap('gist_heat')

    images = [image1, ratio, image2]
    titles = ['OLD', 'RATIO', 'NEW']

    fig, frames = plt.subplots(ncols=len(images), nrows=1, figsize=(15, 5),
                               sharex='all', sharey='all')

    for it in range(len(images)):

        frame = frames[it]
        im, _ = imshow_norm(images[it], frame, origin='lower', aspect='auto',
                            interval=ZScaleInterval(), stretch=LinearStretch(),
                            cmap=cmap1, interpolation='None', rasterized=True)

        if minpoints is not None and maxpoints is not None and labels is not None:
            for pos in range(len(labels)):
                frame.add_patch(Rectangle((min_points[pos][1], min_points[pos][0]),
                                          max_points[pos][1] - min_points[pos][1],
                                          max_points[pos][0] - min_points[pos][0],
                                          edgecolor='yellow', facecolor='None',
                                          fill=False))

        divider = make_axes_locatable(frame)
        cax = divider.append_axes(side, '5%', pad=pad)
        cbar = fig.colorbar(im, cax=cax, orientation=orientation)
        cbar.ax.tick_params(labelsize=8)

        frame.set(title=titles[it])

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

            if DATATYPE == 'RAW':
                old_file = old_files[(mode, fiber)]
                new_file = new_files[(mode, fiber)]
            else:
                old_file = old_files[(mode, fiber)].replace('.fits', '_pp.fits')
                new_file = new_files[(mode, fiber)].replace('.fits', '_pp.fits')

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

            for pos in range(len(labels)):
                print(f'{labels[pos]}'
                      f'\txmin,ymin={min_points[pos][1]},{min_points[pos][0]}'
                      f'\txmax,ymax={max_points[pos][1]},{max_points[pos][0]}')

            print('\n')

            summary_stat = f'{mode} mode fiber {fiber}: '
            stats = []

            for pos in range(len(labels)):

                # mask positions
                mask = xindices > min_points[pos][1]
                mask &= xindices < max_points[pos][1]
                mask &= yindices > min_points[pos][0]
                mask &= yindices < max_points[pos][0]

                # sum pixels in old
                old_value = np.nanpercentile(old_data[mask], percentile)
                new_value = np.nanpercentile(new_data[mask], percentile)

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
            old_90 = np.nanpercentile(old_data, percentile)
            new_90 = np.nanpercentile(new_data, percentile)

            ratio_90 = new_90 / old_90

            print(f'FULL image P[{percentile}]',
                  f'\told={old_90:8.3f} new={new_90:8.3f} ratio={ratio_90:.5f}')

            # summary stats
            if ratio_90 < 1:
                rstat = (1 - ratio_90) * 100
                stat = f' {rstat:.2f}% less flux {percentile}th Percentile'
            else:
                rstat = (ratio_90 - 1) * 100
                stat = f' {rstat:.2f}% more flux {percentile}th Percentile'
            stats.append(stat)

            # add to the stats
            summary_stat += ', '.join(stats)

            summary_text.append(summary_stat)

            # diff image
            with warnings.catch_warnings(record=True) as _:
                ratio_full = (new_data / old_data)

            if PLOT:
                plot_ratio(old_data, new_data, ratio_full, title=f'NIRPS_{mode} Fiber {fiber}',
                           side='right', minpoints=min_points, maxpoints=max_points, labels=labels)

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
