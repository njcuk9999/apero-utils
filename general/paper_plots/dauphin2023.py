import os

import matplotlib
from astropy.io import fits
from astropy.visualization import LinearStretch
from astropy.visualization import imshow_norm, ZScaleInterval

# noinspection PyPep8Naming
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from apero.core import constants

# =============================================================================
# Define variables
# =============================================================================
NIGHT = '2020-08-31'
PLOT_PATH = '/spirou/cook/paper_plots/dauphin2023/'
ODOCODE = '2510385o'

# =============================================================================
# PLOT functions
# =============================================================================
def _norm_image(image, frame, cmap):
    im, norm = imshow_norm(image, frame, origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           cmap=cmap, interpolation='None', rasterized=True)
    return im, norm


def _add_colorbar(fig, im, frame, side='bottom', pad=0.0):
    if side in ['right', 'left']:
        orientation = 'vertical'
    else:
        orientation = 'horizontal'

    divider = make_axes_locatable(frame)
    cax = divider.append_axes(side, '5%', pad=pad)
    cbar = fig.colorbar(im, cax=cax, orientation=orientation)
    cbar.ax.tick_params(labelsize=8)

    return cbar


def plot_fiber_layout(params):
    odocode = ODOCODE

    zoom0area = [1700, 2000, 2900, 3200]

    # get file paths
    raw_file = os.path.join(params['DRS_DATA_RAW'], NIGHT,
                            '{0}.fits'.format(odocode))

    # get
    print('Loading raw image')
    raw_image = fits.getdata(raw_file)

    print('Plotting raw_zoom plot')
    plt.close()

    fig, frame1 = plt.subplots(ncols=1, nrows=1, figsize=(12, 12))

    cmap = matplotlib.cm.get_cmap('inferno').copy()
    cmap.set_bad(color='green')

    # -------------------------------------------------------------------------
    # zoom raw
    zraw_image = raw_image[zoom0area[2]:zoom0area[3], zoom0area[0]:zoom0area[1]]
    # left raw image
    _, _ = _norm_image(zraw_image, frame1, cmap)
    # add labels
    frame1.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)

    frame1.set_title('raw', loc='left',
                     x=0.05, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')

    frame1.text(50, 140, 'A', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame1.text(64, 115, 'B', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame1.text(79, 90, 'C', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)

    frame1.text(50 + 54, 140, 'A', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame1.text(64 + 54, 115, 'B', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame1.text(79 + 54, 90, 'C', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)

    frame1.text(50 + 111, 140, 'A', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame1.text(64 + 111, 115, 'B', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame1.text(79 + 111, 90, 'C', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    # -------------------------------------------------------------------------
    plt.subplots_adjust(wspace=0.05, hspace=0.05,
                        left=0.01, right=0.99, top=0.99, bottom=0.01)
    # save file
    outfile = os.path.join(PLOT_PATH, 'raw_zoom.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile, dpi=300)
    print('Showing graph')
    plt.show()
    plt.close()



def plot_raw_features(params):
    # set file to use
    odocode = ODOCODE
    # get file paths
    raw_file = os.path.join(params['DRS_DATA_RAW'], NIGHT,
                            '{0}.fits'.format(odocode))
    # get raw image
    print('Loading raw image')
    raw_image = fits.getdata(raw_file)
    # plot feature grid
    plot_feature_grid(raw_image, outname='raw_full.pdf')

def plot_feature_grid(image, outname):
    # -------------------------------------------------------------------------
    print('Plotting raw_full plot')
    # set up grid
    plt.close()
    fig, frame1 = plt.subplots(ncols=1, nrows=1, figsize=(12, 12))
    # get colour map
    cmap = matplotlib.cm.get_cmap('inferno').copy()
    cmap.set_bad(color='green')

    cmap0 = matplotlib.cm.get_cmap('Greys_r').copy()
    cmap0.set_bad(color='green')

    # top left raw image
    _ = _norm_image(image, frame1, cmap)
    # add labels
    frame1.set(xlim=(0, image.shape[1]), ylim=(0, image.shape[0]))
    frame1.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    frame1.set_title('raw ({0}x{1})'.format(*image.shape), loc='left',
                     x=0.05, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')
    # -------------------------------------------------------------------------
    plt.subplots_adjust(wspace=0.05, hspace=0.05,
                        left=0.01, right=0.99, top=0.99, bottom=0.01)
    # -------------------------------------------------------------------------
    outfile = os.path.join(PLOT_PATH, outname)
    print('Saving to file: ' + outfile)
    plt.savefig(outfile, dpi=300)
    print('Showing graph')
    plt.show()
    plt.close()


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':

    _params = constants.load()

    plot_fiber_layout(_params)

    plot_raw_features(_params)

# =============================================================================
# End of code
# =============================================================================