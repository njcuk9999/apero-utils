from astropy.io import fits
import matplotlib.pyplot as plt
import argparse
import glob
import os
import numpy as np
from astropy.visualization import ZScaleInterval, LinearStretch, ImageNormalize
from tqdm import tqdm
from astroquery.vizier import Vizier


# =============================================================================
# Define variables
# =============================================================================

PATH = '/data/spip/drs-data/spip_test/red/2025-07-07'

SCI_FIBER = 'AB'

HC_FILE = '*c_pp_e2dsff_{FIBER}.fits'

WAVE_FILE = '*_pp_e2dsff_{FIBER}_wave_night_{FIBER}.fits'

E2DSFF_FILE = 'DEBUG*a_pp_e2dsll_{FIBER}.fits'

E2DSLL_FLAT = 'DEBUG_*_pp_e2dsll_{FIBER}.fits'

# fetch that info from params
NORDERS = 49

# Added wavelength beyond overlap
WAVE_OVERLAP_EXTRA = 10

# Orders to plot for HC overlap (set to None for all)
WAVE_ORDERS = [1, 10, 20, 30, 40, 47]

# HC model from Vizier
# BY default this is the full table reference for J/ApJS/195/24, Table 2
VIZIER_HC_MODEL = 'J/ApJS/195/24/table2'

HC_WAVE_COL = 'lambda'
HC_SPECIES_COL = 'Ion'
HC_SPECIES = ['UI', 'UII', 'Ne']
HC_SPECIES_COLORS = dict()
HC_SPECIES_COLORS['UI'] = 'blue'
HC_SPECIES_COLORS['UII'] = 'cyan'
HC_SPECIES_COLORS['Ne'] = 'green'



# =============================================================================
# Define functions
# =============================================================================
def get_hc_model():
    # No row limit: get all rows
    Vizier.ROW_LIMIT = -1
    # Query the entire table
    tables = Vizier.get_catalogs(VIZIER_HC_MODEL)
    # return this table
    return tables[0]


def get_hc_lines(table, species, wavemin, wavemax):
    # get a mask just for this species
    smask = table[HC_SPECIES_COL] == species
    # only keep those in the wavemap wavelength range
    smask &= table[HC_WAVE_COL] >= wavemin
    smask &= table[HC_WAVE_COL] <= wavemax

    if np.sum(smask) == 0:
        return None
    else:
        return table[smask]

def plot_hc_overlap():

    # find hc file
    hc_glob = os.path.join(PATH, HC_FILE.format(FIBER=SCI_FIBER))
    hc_files = glob.glob(hc_glob)

    # find wave sols
    wave_glob = os.path.join(PATH, WAVE_FILE.format(FIBER=SCI_FIBER))
    wave_files = glob.glob(wave_glob)

    if len(hc_files) == 0:
        print(f'No HC file found for {hc_glob}')
        return
    if len(wave_files) == 0:
        print(f'No wave file {wave_glob}')
        return

    # load hc file
    hc_data = fits.getdata(hc_files[0])
    # load wave file
    wave_data = fits.getdata(wave_files[0])

    # get the HC model
    hc_model_table = get_hc_model()


    for order_num in range(len(wave_data))[:-1]:

        if WAVE_ORDERS is not None:
            if order_num not in WAVE_ORDERS:
                continue

        plt.close()
        fig, frame = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

        if order_num % 2 == 0:
            color1 = 'orange'
            color2 = 'purple'
        else:
            color1 = 'purple'
            color2 = 'orange'


        wave1 = wave_data[order_num]
        wave2 = wave_data[order_num + 1]
        flux1 = hc_data[order_num]
        flux2 = hc_data[order_num + 1]

        cut1 = np.nanmin(wave2[np.isfinite(flux2)]) - WAVE_OVERLAP_EXTRA
        cut2 = np.nanmax(wave1[np.isfinite(flux1)]) + WAVE_OVERLAP_EXTRA


        mask1 = (wave1 > cut1) & (wave1 < cut2)
        mask2 = (wave2 > cut1) & (wave2 < cut2)

        # plot all species
        for species in HC_SPECIES:
            # get all HC lines from the model for this species at these
            #  wavelengths
            hc_table_ord = get_hc_lines(hc_model_table, species, cut1, cut2)
            # deal with no HC lines for this species
            if hc_table_ord is None:
                continue
            # loop around all hc lines
            for row in range(len(hc_table_ord)):
                if row == 0:
                    label = f'{VIZIER_HC_MODEL} [{species}]'
                else:
                    label = None
                # plot the hc lines
                frame.axvline(x=hc_table_ord[HC_WAVE_COL][row],
                              color=HC_SPECIES_COLORS[species],
                              ls='--', label=label, alpha=0.25)

        frame.plot(wave1[mask1], np.sqrt(np.abs(flux1[mask1])),
                   color=color1, label=f'Order {order_num}',
                   linestyle='-', marker='.')
        frame.plot(wave2[mask2], np.sqrt(np.abs(flux2[mask2])),
                   color=color2, label=f'Order {order_num + 1}',
                   linestyle='-', marker='.')

        frame.set(xlabel='Wavelength [nm]', ylabel=r'sqrt(Flux)',
                  title=f'Orders [{order_num}, {order_num+1}]')

        frame.legend(loc=0)


        plt.show()


def plot_e2dsff():

    # find hc file
    e2dsff_glob = os.path.join(PATH, E2DSFF_FILE.format(FIBER=SCI_FIBER))
    e2dsff_files = glob.glob(e2dsff_glob)

    if len(e2dsff_files) == 0:
        print(f'No E2DSFF file found for {e2dsff_glob}')
        return

    e2dsff_data = fits.getdata(e2dsff_files[0])

    plt.close()
    fig, frame = plt.subplots(ncols=1, nrows=1, figsize=(10, 10))
    # -------------------------------------------------------------------------
    # Apply DS9-style zscale and linear stretch
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(e2dsff_data)
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())
    # Mask NaNs
    pp_image_masked = np.ma.masked_invalid(e2dsff_data)

    # Plot pp_image using imshow with normalization
    # Choose a colormap and set NaN color to green
    cmap = plt.cm.get_cmap('inferno').copy()
    cmap.set_bad(color='green')

    plt.imshow(pp_image_masked, aspect='auto', interpolation='nearest',
               origin='lower', cmap=cmap, norm=norm)

    plt.show()
    plt.close()




def plot_e2dsll_flat():
    # find e2dsll flat file
    e2dsll_flat_glob = os.path.join(PATH, E2DSLL_FLAT.format(FIBER=SCI_FIBER))
    e2dsll_flat_files = glob.glob(e2dsll_flat_glob)

    if len(e2dsll_flat_files) == 0:
        print(f'No E2DSLL flat file found for {e2dsll_flat_glob}')
        return

    # load wave file
    image = fits.getdata(e2dsll_flat_files[0])

    # find the middle of the array
    mid = image.shape[1]//2

    # median trace profile of the center of the image
    med =  np.nanmedian(image[:,mid-100:mid+100],axis=1)

    # Reshape the median profile to have NORDERS rows
    med = med.reshape(NORDERS, med.shape[0]//NORDERS)

    # normalize each order to a mean of 1
    for iord in range(NORDERS):
        segment = med[iord]
        segment /= np.nansum(segment)

    # we find the flux at the edges of the trace for each order. The total
    # edge flux should be small and account for <1% of the total flux.
    flux_edge = med[:,0] + med[:,-1]
    flux_left = med[:,0]
    flux_right = med[:,-1]

    fig, frames = plt.subplots(figsize=(12, 6), nrows=1, ncols=2)

    frames[0].imshow(med.T, aspect='auto', cmap='inferno',
               interpolation='nearest', origin='lower',
               extent=[0, NORDERS, 0, med.shape[1]])
    #ax[0].colorbar(label='Normalized Flux')
    frames[0].set_title('Line List Image')
    frames[0].set_ylabel('Pixel')
    frames[0].set_xlabel('Order Number')

    frames[1].plot(flux_edge, label='Total Edge Flux')
    frames[1].plot(flux_left, label='Left Edge Flux',alpha=0.5)
    frames[1].plot(flux_right, label='Right Edge Flux',alpha=0.5)
    frames[1].set_xlabel('Order Number')
    frames[1].set_ylabel('Edge Flux')
    frames[1].set_title('Edge Flux of Each Order in the Line List')
    frames[1].grid()
    frames[1].axhline(y=0.01, color='r', linestyle='--', label='1% Threshold')
    frames[1].legend()
    plt.tight_layout()
    plt.show()


# =============================================================================
# main code
# =============================================================================
if __name__ == '__main__':

    # check
    parser = argparse.ArgumentParser(description='Calibration checks')
    parser.add_argument('--hc', action='store_true', help='Run HC overlap')
    parser.add_argument('--e2dsff', action='store_true', help='Run E2DSFF plot')
    parser.add_argument('--e2dsll', action='store_true', help='Run E2DSLL plot')
    args = parser.parse_args()

    # check wavelength (with HC overlap)
    if args.hc:
        plot_hc_overlap()
    # plot e2dsff
    if args.e2dsff:
        plot_e2dsff()
    # check extraction
    if args.e2dsll:
        plot_e2dsll_flat()

