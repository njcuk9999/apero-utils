from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
import os
import numpy as np
from astropy.visualization import ZScaleInterval, LinearStretch, ImageNormalize
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
FIBERS = ['A', 'B', 'C']
LOC_FIBER_NAME = dict()
LOC_FIBER_NAME['A'] = 'AB'
LOC_FIBER_NAME['B'] = 'AB'
LOC_FIBER_NAME['C'] = 'C'

TMP_PATH = '/data/spip/drs-data/spip_test/tmp/'
CALIB_PATH = '/data/spip/drs-data/spip_test/calib/'

LOC_FILE = '*_loco_{fiber}.fits'

PP_FILE = '2025-07-07/189403a_pp.fits'

# Defines the resized image
IMAGE_X_LOW = 4
IMAGE_X_HIGH = 4092
IMAGE_Y_LOW = 250
IMAGE_Y_HIGH = 3350


CUT_DOWN_SIZE_X = [1700, 2000]
CUT_DOWN_SIZE_Y = [2300, 2600]

# =============================================================================
# Define functions
# =============================================================================
def ab_to_a_b(data, fiber):

    if fiber in ['AB', 'C']:
        return data

    if fiber == 'A':
        return data[1::2, :]
    else:
        return data[:-1:2, :]


def cut_down(image: np.ndarray):
    return image[CUT_DOWN_SIZE_Y[0]:CUT_DOWN_SIZE_Y[1],
           CUT_DOWN_SIZE_X[0]:CUT_DOWN_SIZE_X[1]]

def pp_to_loco(image: np.ndarray):
    # flip image in both directions
    newimage = image[::-1, ::-1]
    # resize image
    xtake = np.arange(IMAGE_X_LOW, IMAGE_X_HIGH)
    ytake = np.arange(IMAGE_Y_LOW, IMAGE_Y_HIGH)
    newimage = np.take(np.take(newimage, xtake, axis=1), ytake, axis=0)
    # return new image
    return newimage


def plot_with_labels(frame, x, y, label, xlim, ylim):
    # Plot using 'A' as marker
    frame.plot(x, y, linestyle='None')  # no line between points
    for xi, yi in zip(x, y):
        if not (xlim[0] <= xi <= xlim[1] and ylim[0] <= yi <= ylim[1]):
            continue
        if xi % 64 == 0:
            frame.text(xi, yi, label, fontsize=14, ha='center', va='center',
                       bbox=dict(facecolor='white', edgecolor='black', pad=1.0))


# =============================================================================
# main code
# =============================================================================
if __name__ == '__main__':
    # storage
    loc_data = dict()
    # -------------------------------------------------------------------------
    # find files for each fiber
    for fiber in FIBERS:
        # get the base name for this loc
        _basename = LOC_FILE.format(fiber=LOC_FIBER_NAME[fiber])
        # get all files that match this fiber
        _files = glob.glob(os.path.join(CALIB_PATH, _basename))
        # assume the first one is good
        _loc_data = fits.getdata(_files[0])
        # get AB to A+B
        loc_data[fiber] = ab_to_a_b(_loc_data, fiber)
    # -------------------------------------------------------------------------
    # load pp file
    pp_data = fits.getdata(os.path.join(TMP_PATH, PP_FILE))
    # convert pp image to loco frame
    pp_image = pp_to_loco(pp_data)
    # cut down the image
    pp_image = cut_down(pp_image)
    # -------------------------------------------------------------------------
    # plotting
    # -------------------------------------------------------------------------
    fig, frame = plt.subplots(ncols=1, nrows=1, figsize=(10, 10))
    # -------------------------------------------------------------------------
    # Apply DS9-style zscale and linear stretch
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(pp_image)
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())
    # Mask NaNs
    pp_image_masked = np.ma.masked_invalid(pp_image)

    # Plot pp_image using imshow with normalization
    # Choose a colormap and set NaN color to green
    cmap = plt.cm.get_cmap('inferno').copy()
    cmap.set_bad(color='green')

    plt.imshow(pp_image_masked, origin='lower', cmap=cmap, norm=norm,
               extent=CUT_DOWN_SIZE_X+CUT_DOWN_SIZE_Y)
    # -------------------------------------------------------------------------
    # plot the orders
    for fiber in FIBERS:
        print(f'Plotting fiber {fiber}')
        loco = loc_data[fiber]

        for order_num in tqdm(range(loco.shape[0])):

            xx = np.arange(loco.shape[1])
            yy = loco[order_num]
            # cut down
            cut_mask = yy > CUT_DOWN_SIZE_Y[0]
            cut_mask &= yy < CUT_DOWN_SIZE_Y[1]
            xx_cut = xx[cut_mask]
            yy_cut = yy[cut_mask]

            plot_with_labels(frame, xx_cut, yy_cut, label=fiber,
                             xlim=CUT_DOWN_SIZE_X, ylim=CUT_DOWN_SIZE_Y)

    frame.set(xlim=CUT_DOWN_SIZE_X, ylim=CUT_DOWN_SIZE_Y)

    plt.show()
    plt.close()


# =============================================================================
# end of  code
# =============================================================================
