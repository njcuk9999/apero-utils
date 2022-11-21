#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-11-20 at 19:19

@author: cook
"""
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
labels = ['center', 'left', 'right', 'red', 'blue']
boxsize = 5
function = np.nanmedian
func_name = 'Med'
PLOT = True

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

            if fiber == 'A' and mode == 'HE':
                old_file = '/nirps_raw/nirps/raw-data/nirps_he/2022-11-17/NIRPS_2022-11-17T22_15_02_461.fits'
                new_file = '/nirps_raw/nirps/raw-data/nirps_he/2022-11-19/NIRPS_2022-11-19T23_32_10_890.fits'
                xpos = [2124, 2259, 2358, 3515, 235]
                ypos = [1952, 4022, 85, 1977, 1900]
            elif fiber == 'B' and mode == 'HE':
                old_file = '/nirps_raw/nirps/raw-data/nirps_he/2022-11-17/NIRPS_2022-11-17T22_16_42_792.fits'
                new_file = '/nirps_raw/nirps/raw-data/nirps_he/2022-11-19/NIRPS_2022-11-19T23_33_51_221.fits'
                xpos = [2145, 2280, 2380, 3535, 260]
                ypos = [1952, 4022, 85, 1977, 1900]
            elif fiber == 'A' and mode == 'HA':
                old_file = '/nirps_raw/nirps/raw-data/nirps_ha/2022-11-17/NIRPS_2022-11-17T15_22_57_043.fits'
                new_file = '/nirps_raw/nirps/raw-data/nirps_ha/2022-11-19/NIRPS_2022-11-19T23_28_44_649.fits'
                xpos = [2124, 2259, 2358, 3515, 235]
                ypos = [1952, 4022, 85, 1977, 1900]
            elif fiber == 'B' and mode == 'HA':
                old_file = '/nirps_raw/nirps/raw-data/nirps_ha/2022-11-17/NIRPS_2022-11-17T15_24_37_374.fits'
                new_file = '/nirps_raw/nirps/raw-data/nirps_ha/2022-11-19/NIRPS_2022-11-19T23_30_24_979.fits'
                xpos = [2145, 2280, 2380, 3535, 260]
                ypos = [1952, 4022, 85, 1977, 1900]
            else:
                raise ValueError(f'mode={mode} fiber={fiber} not supported')

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
            print(f'old={old_hdr["ORIGFILE"]}  EXPTIME={old_hdr["EXPTIME"]}')
            print(f'new={new_hdr["ORIGFILE"]}  EXPTIME={new_hdr["EXPTIME"]}')

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
            diff = new_data - old_data

            if PLOT:
                plot_image(diff, title=f'NIRPS_{mode} Fiber {fiber}')

    # print the summary
    print('\n\n Summary: \n\n')
    for stext in summary_text:
        print(stext + '\n')

# =============================================================================
# End of code
# =============================================================================
