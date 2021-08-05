import matplotlib
matplotlib.use('Qt5Agg')


from astropy.io import fits
from astropy.visualization import imshow_norm, ZScaleInterval, LinearStretch
import glob
import matplotlib.pyplot as plt
import numpy as np
import os

from apero.core import constants
from apero.io import drs_image


NIGHT = '2020-08-31'
PLOT_PATH = '/data/spirou/drs-data/misc/paper_plots'
# PLOT_PATH = '/scratch2/drs-data/misc/paper_plots'

PLOTS = []
# PLOTS.append('SIZE_GRID')
PLOTS.append('BADMAP')



def plot_size_grid(params):

    odocode = '2510288a'
    hashcode = '3444961B5D'
    # get file paths
    raw_file = os.path.join(params['DRS_DATA_RAW'], NIGHT,
                            '{0}.fits'.format(odocode))
    pp_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                           '{0}_pp.fits'.format(hashcode))
    e2ds_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                             '{0}_pp_e2dsff_AB.fits'.format(hashcode))
    # get
    print('Loading raw image')
    raw_image = fits.getdata(raw_file)
    print('Loading pp image')
    pp_image = fits.getdata(pp_file)
    print('Loading E2DS image')
    e2ds_image = fits.getdata(e2ds_file)

    # rotation to match HARPS orientation (expected by DRS)
    # image1 = drs_image.rotate_image(raw_image, params['RAW_TO_PP_ROTATION'])
    # flip image
    image2 = drs_image.flip_image(params, pp_image)
    # get resize size
    sargs = dict(xlow=params['IMAGE_X_LOW'], xhigh=params['IMAGE_X_HIGH'],
                 ylow=params['IMAGE_Y_LOW'], yhigh=params['IMAGE_Y_HIGH'])
    # resize flat
    image3 = drs_image.resize(params, image2, **sargs)

    print('Plotting size_grid plot')
    fig = plt.figure(figsize=(10, 15))
    size = (3, 2)
    frame1 = plt.subplot2grid(size, (0, 0))
    frame2 = plt.subplot2grid(size, (0, 1))
    frame3 = plt.subplot2grid(size, (1, 0))
    frame4 = plt.subplot2grid(size, (1, 1))
    frame5 = plt.subplot2grid(size, (2, 0), colspan=2)

    cmap = matplotlib.cm.get_cmap('inferno').copy()
    cmap.set_bad(color='green')

    # top left raw image
    im, norm = imshow_norm(raw_image, frame1, origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           cmap=cmap, interpolation='None', rasterized=True)
    # top right: pp image (rotated)
    im, norm = imshow_norm(pp_image, frame2, origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           cmap=cmap, interpolation='None', rasterized=True)
    # middle left: flipped image
    im, norm = imshow_norm(image2, frame3, origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           cmap=cmap, interpolation='None', rasterized=True)
    # middle right: flipped + resized image
    im, norm = imshow_norm(image3, frame4, origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           cmap=cmap, interpolation='None', rasterized=True)
    # bottom: e2ds
    im, norm = imshow_norm(e2ds_image, frame5, origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           cmap=cmap, interpolation='None', rasterized=True)

    # add labels
    frame1.set(xlim=(0, 4096), ylim=(0, 4096))
    frame1.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    frame1.set_title('raw (4096x4096)', loc='left',
                     x=0.05, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')


    frame2.set(xlim=(0, 4096), ylim=(0, 4096))
    frame2.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    frame2.set_title('pre-processed (4096x4096)', loc='left',
                     x=0.05, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')

    frame3.set(xlim=(0, 4096), ylim=(0, 4096))
    frame3.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    frame3.set_title('pre-processed flipped (4096x4096)', loc='left',
                     x=0.05, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')

    frame4.set(xlim=(-4, 4092), ylim=(0, 4096))
    frame4.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    frame4.set_title('pre-processed, flipped, resized (3100x4088)', loc='left',
                     x=0.05, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')

    frame5.set(xlim=(0, 4088), ylim=(0, 49))
    frame5.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    frame5.set_title('Extracted (E2DS) 49x4088', loc='left',
                     x=0.025, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')

    plt.subplots_adjust(wspace=0.0075, hspace=0.01,
                        left=0.01, right=0.99, top=0.975, bottom=0.025)

    outfile = os.path.join(PLOT_PATH, 'size_grid.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile)
    print('Showing graph')
    plt.show()
    plt.close()



def plot_badpix_plot(params):

    hashcode = 'DB67D5C4F5'
    xlow = params['IMAGE_X_LOW']
    xhigh = params['IMAGE_X_HIGH']
    ylow = params['IMAGE_Y_LOW']
    yhigh = params['IMAGE_Y_HIGH']

    dark_file = glob.glob(os.path.join(params['DRS_DATA_REDUC'], 'other',
                                       '*d_pp_dark_master.fits'))[0]
    bad_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                            '{0}_pp_badpixel.fits'.format(hashcode))

    dark_image = fits.getdata(dark_file)
    bad_image = fits.getdata(bad_file).astype(float)

    # fill bad image
    bad_image_full = np.zeros_like(dark_image)
    bad_image_full[ylow:yhigh, xlow:xhigh] = bad_image
    dark_image = drs_image.flip_image(params, dark_image)

    cmap1 = matplotlib.cm.get_cmap('Greys_r').copy()
    cmap2 = matplotlib.cm.get_cmap('Greys').copy()

    plt.close()
    fig, frames = plt.subplots(figsize=(20, 10), ncols=2, nrows=1)
    frame0 = frames[0]
    frame1 = frames[1]

    im, norm = imshow_norm(dark_image, frame0, origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           cmap=cmap1, interpolation='None', rasterized=True)
    im, norm = imshow_norm(bad_image_full, frame1, origin='lower', aspect='auto',
                           cmap=cmap2, interpolation='None', rasterized=True)



    frame0.tick_params(axis='both', which='both', bottom=False, top=False,
                      left=False, right=False, labelleft=False,
                      labelbottom=False)
    frame1.tick_params(axis='both', which='both', bottom=False, top=False,
                      left=False, right=False, labelleft=False,
                      labelbottom=False)

    frame0.hlines(y=ylow, xmin=xlow, xmax=xhigh, color='r', lw=2)
    frame0.hlines(y=yhigh, xmin=xlow, xmax=xhigh, color='r', lw=2)
    frame0.vlines(x=xlow, ymin=ylow, ymax=yhigh, color='r', lw=2)
    frame0.vlines(x=xhigh, ymin=ylow, ymax=yhigh, color='r', lw=2)
    frame1.hlines(y=ylow, xmin=xlow, xmax=xhigh, color='r', lw=2)
    frame1.hlines(y=yhigh, xmin=xlow, xmax=xhigh, color='r', lw=2)
    frame1.vlines(x=xlow, ymin=ylow, ymax=yhigh, color='r', lw=2)
    frame1.vlines(x=xhigh, ymin=ylow, ymax=yhigh, color='r', lw=2)

    frame0.set(xlim=[0, 4096], ylim=[0, 4096])
    frame1.set(xlim=[0, 4096], ylim=[0, 4096])

    plt.subplots_adjust(wspace=0, hspace=0, left=0.01, right=0.99,
                        bottom=0.01, top=0.99)

    outfile = os.path.join(PLOT_PATH, 'badmap.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile)
    print('Showing graph')
    plt.show()
    plt.close()










if __name__ == '__main__':

    params = constants.load()

    if 'SIZE_GRID' in PLOTS:
        plot_size_grid(params)
    if 'BADMAP' in PLOTS:
        plot_badpix_plot(params)