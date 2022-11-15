#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare HC 2D preprocessed images with each other (side-by-side)
with zoom in

Created on 2022-11-15 at 13:58

@author: cook
"""
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ImageNormalize
from astropy.visualization import LinearStretch
from astropy.visualization import ZScaleInterval
from mpl_toolkits.axes_grid1 import make_axes_locatable

# =============================================================================
# Define variables
# =============================================================================
path2 = ('/scratch2/nirps_ha/drs-data/nirps_ha_202211/tmp/2022-11-07/'
         'NIRPS_2022-11-07T15_34_08_871_pp.fits')
path1 = ('/scratch2/nirps_ha/drs-data/nirps_ha_202205/tmp/2022-06-09/'
         'NIRPS_2022-06-10T12_29_47_622_pp.fits')
extent = [2250, 2550, 2300, 2600]
titles = ['2022-06-09', '2022-11-07',
          '2022-06-09  [Zoom]', '2022-11-07  [Zoom]']

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    image1 = fits.getdata(path1)
    image2 = fits.getdata(path2)

    image3 = image1[extent[2]:extent[3], extent[0]:extent[1]]
    image4 = image2[extent[2]:extent[3], extent[0]:extent[1]]

    fig, frames = plt.subplots(ncols=2, nrows=2, figsize=(12, 12))

    frames = frames.ravel()
    images = [image1, image2, image3, image4]
    extents = [None, None, extent, extent]

    for it, frame in enumerate(frames):
        image = images[it]
        extent = extents[it]

        norm = ImageNormalize(image, interval=ZScaleInterval(),
                              stretch=LinearStretch())

        divider = make_axes_locatable(frame)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        im = frame.imshow(image, origin='lower', norm=norm,
                          cmap='inferno', extent=extents[it])
        frame.set(title=titles[it])

        fig.colorbar(im, cax=cax, orientation='vertical')

    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)
    plt.savefig('/scratch2/nirps/nirps_ha_HC_HC_comp.pdf', dpi=300)

# =============================================================================
# End of code
# =============================================================================
