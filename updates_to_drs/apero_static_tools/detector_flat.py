from astropy.io import fits  
import matplotlib.pyplot as plt  
import numpy as np 
from scipy.ndimage import zoom as scizoom  # For resizing (zooming) images
import os  
from astropy.table import Table  
from tqdm import tqdm

############################ NUMERICAL PARAMETERS ############################

namp = 32  # Number of amplifiers on the detector (used for amplifier folding/corrections)

binsize = 32  # Size of bins for spatial binning (e.g., 32x32 pixels per bin)

frac_flat_bad = 0.2  # Fractional threshold for flat field bad pixels (relative deviation from 1)

dark_threshold = 1  # Threshold (in e-/s) for dark current to consider a pixel "hot"

doplot = True  # If True, plots will be shown for visual inspection

########################## FILES #############################################

# List of dark frame FITS files (used for master dark and dark current fitting)
# Should have a range of exposure times to cover the dark current behavior

dark_files = [
    'NIRPS_PATCHED_2022-11-16T17_32_50_169.fits',
    'NIRPS_PATCHED_2022-11-16T17_46_51_728.fits',
    'NIRPS_PATCHED_2022-11-16T18_00_53_288.fits',
    'NIRPS_PATCHED_2022-11-16T18_14_54_848.fits',
    'NIRPS_PATCHED_2022-11-16T18_28_56_408.fits',
    'NIRPS_PATCHED_2022-11-17T00_46_33_198.fits',
    'NIRPS_PATCHED_2022-11-17T00_47_45_651.fits',
    'NIRPS_PATCHED_2022-11-17T00_48_58_105.fits',
    'NIRPS_PATCHED_2022-11-17T00_50_10_559.fits',
    'NIRPS_PATCHED_2022-11-17T00_51_23_012.fits'
]

# List of LED flat field FITS files (used for master flat)
led_files = [
    'NIRPS_PATCHED_2022-11-16T19_54_34_959.fits',
    'NIRPS_PATCHED_2022-11-16T19_55_30_693.fits',
    'NIRPS_PATCHED_2022-11-16T19_56_26_427.fits',
    'NIRPS_PATCHED_2022-11-16T19_57_22_161.fits',
    'NIRPS_PATCHED_2022-11-16T19_58_17_894.fits',
    'NIRPS_PATCHED_2022-11-16T19_59_58_223.fits',
    'NIRPS_PATCHED_2022-11-16T20_01_38_543.fits',
    'NIRPS_PATCHED_2022-11-16T20_03_18_863.fits'
]

path_to_data = 'data/'

# not used, just for reference
sample_dark = 'sample_dark.fits'  # Output file for the median dark frame

sample_dark_table = 'sample_dark_table.csv'  # Output CSV for hot pixel table

# detector high-passed flat field
flat_outname = 'flat_field.fits'  # Output file for the processed flat field

# median LED frame without filtering. Not used, just for reference
sample_led = 'sample_led.fits'  # Output file for the median LED frame

# amplifier bias model
dark_intercept = 'dark_intercept.fits'  # Output file for dark current intercept map
dark_slope = 'dark_slope.fits'  # Output file for dark current slope map

####### END OF USER INPUT #####################################################

# add path to data files
dark_files = [os.path.join(path_to_data, f) for f in dark_files]
led_files = [os.path.join(path_to_data, f) for f in led_files]


def sigma(im):
    """
    Compute a robust estimate of the standard deviation of an image,
    using the 16th and 84th percentiles (ignoring NaNs).
    """
    n1, p1 = np.nanpercentile(im, [16, 84])
    return (p1 - n1) / 2

def medbin(im, binsize=32):
    """
    Bin an image into non-overlapping bins of size binsize x binsize,
    taking the median value in each bin (ignoring NaNs).
    Returns the binned image.
    """
    ny, nx = im.shape  # Get image dimensions
    n_bins_y = ny // binsize  # Number of bins along y
    n_bins_x = nx // binsize  # Number of bins along x
    binned_image = np.zeros((n_bins_y, n_bins_x), dtype=np.float32)
    for i in range(binned_image.shape[0]):
        for j in range(binned_image.shape[1]):
            # Take the median of each bin, ignoring NaNs
            binned_image[i, j] = np.nanmedian(im[i * binsize:(i + 1) * binsize, j * binsize:(j + 1) * binsize])
    return binned_image

# Convert file lists to numpy arrays for easier indexing
dark_files = np.array(dark_files)
led_files = np.array(led_files)

# Pre-allocate array for exposure times of dark frames
exptimes = np.zeros(len(dark_files))
# Loop through all dark files to extract their exposure times from FITS headers
for i, file in enumerate(dark_files):
    print('Reading dark frame:', file)
    hdr = fits.getheader(file)
    exptimes[i] = hdr['EXPTIME']

# If the median dark frame does not exist, create it from the longest exposure darks
if not os.path.exists(sample_dark):
    print('Computing median dark frame from the longest exposure times...')
    long_darks = dark_files[exptimes == np.max(exptimes)]  # Select only the longest exposure darks
    sz = fits.getdata(long_darks[0]).shape  # Get image shape from first file
    cube_dark = np.zeros((len(long_darks), sz[0], sz[1]), dtype=float)  # Allocate cube for stacking
    for i, file in enumerate(long_darks):
        cube_dark[i] = fits.getdata(file)  # Read each dark frame into the cube
    dark = np.nanmedian(cube_dark, axis=0)  # Take the median across the stack
    fits.writeto(sample_dark, dark, overwrite=True)  # Save the median dark frame

dark = fits.getdata(sample_dark)  # Load the median dark frame

# If the median LED frame does not exist, create it by stacking and normalizing all LED frames
if not os.path.exists(sample_led):
    print('Creating LED frame table...')
    cube = np.zeros((len(led_files), dark.shape[0], dark.shape[1]), dtype=float)  # Allocate cube for LED frames
    for i, file in enumerate(led_files):
        tmp = fits.getdata(file) - dark  # Subtract dark from each LED frame
        tmp /= np.nanmedian(tmp)  # Normalize by median value
        cube[i] = tmp  # Store in cube
    led = np.nanmedian(cube, axis=0)  # Take median across all LED frames
    fits.writeto(sample_led, led, overwrite=True)  # Save the median LED frame

led = fits.getdata(sample_led)  # Load the median LED frame

def amp_fold(im, namp):
    """
    Fold the image data into a 2D array with 'namp' amplifiers.
    Each amplifier is assumed to be a vertical region in the image.
    For odd-numbered amplifiers, the region is flipped horizontally.
    Returns the median of all amplifier regions.
    """
    sz = im.shape
    amp_size = im.shape[0] // namp  # Size of each amplifier region
    cubeamp = np.zeros((namp, sz[0], sz[1] // namp))
    for i in range(namp):
        slice = im[:, i * amp_size:(i + 1) * amp_size]  # Extract amplifier region
        slice -= np.nanmedian(slice)  # Normalize by median
        if i % 2 == 0:
            slice2 = slice  # Even amplifiers: no flip
        else:
            slice2 = slice[:, ::-1]  # Odd amplifiers: flip horizontally
        cubeamp[i, :, :] = slice2  # Store in cube
    return np.nanmedian(cubeamp, axis=0)  # Return median across amplifiers

# If the flat field output does not exist, compute it from the LED frame
if not os.path.exists(flat_outname):
    print('Computing flat field image...')
    n1, med, p1 = np.nanpercentile(led, [16, 50, 84])  # Compute robust statistics
    s1 = (p1 - n1) / 2  # Robust standard deviation
    # Identify bad pixels in the flat field (outliers and negatives)
    bad = (led < n1 - 3 * s1) | (led > p1 + 3 * s1) | (led < 0)
    led[bad] = np.nan  # Set bad pixels to NaN
    led /= np.nanmedian(led)  # Normalize by median

    # Prepare for iterative large-scale flat correction
    binned_image = np.zeros((led.shape[0] // binsize, led.shape[1] // binsize))
    prev = np.zeros_like(binned_image)  # Previous iteration's binned image
    sig_step = np.inf  # Convergence metric
    ite = 0  # Iteration counter

    # Iteratively divide by a smoothed version of the flat to remove large-scale structure
    while (sig_step > 1e-4):
        print(f'Iteration {ite + 1}, previous step: {sig_step:.4e}')
        binned_image = medbin(led)
        recon = scizoom(binned_image, binsize, order=1)  # Upsample binned image
        led /= recon  # Divide by smoothed image
        sig_step = np.nanmedian(np.abs(binned_image - prev))  # Check for convergence
        prev = np.array(binned_image)  # Store for next iteration
        ite += 1

    # After convergence, mask any pixels that deviate too much from 1
    n1, med, p1 = np.nanpercentile(led, [16, 50, 84])
    sig = (p1 - n1) / 2
    bad = np.abs(led - 1) > frac_flat_bad  # Pixels too far from 1 are bad
    led[bad] = np.nan  # Mask bad pixels
    fits.writeto(flat_outname, led, overwrite=True)  # Save the flat field

flat = fits.getdata(flat_outname)  # Load the flat field

# If the dark current intercept and slope maps do not exist, fit them from all dark frames
if not os.path.exists(dark_intercept) or not os.path.exists(dark_slope):
    cube_amps = np.zeros((len(dark_files), dark.shape[0], dark.shape[1] // namp), dtype=float)
    intercept = np.zeros((dark.shape[0], dark.shape[1] // namp), dtype=float)
    slope = np.zeros((dark.shape[0], dark.shape[1] // namp), dtype=float)
    for i, file in enumerate(dark_files):
        print('Reading dark frame:', file)
        im = fits.getdata(file)
        amp = amp_fold(im, namp)  # Fold into amplifier regions
        cube_amps[i] = amp  # Store in cube

    # For each pixel in the amplifier-folded image, fit a line (dark current vs. exposure time)
    for i in tqdm(range(dark.shape[0]), leave=False, desc='Fitting dark frames'):
        for j in tqdm(range(dark.shape[1] // namp), leave=False):
            coeffs = np.polyfit(exptimes, cube_amps[:, i, j], 1)  # Linear fit
            intercept[i, j] = coeffs[1]  # Offset (bias)
            slope[i, j] = coeffs[0]  # Slope (dark current rate)
    fits.writeto(dark_intercept, intercept, overwrite=True)
    fits.writeto(dark_slope, slope, overwrite=True)

dark_intercept = fits.getdata(dark_intercept)  # Load intercept map
dark_slope = fits.getdata(dark_slope)  # Load slope map

dark = fits.getdata(dark_files[0])  # Load a sample dark frame
hdr = fits.getheader(dark_files[0])  # Load its header
recon_amp = dark_intercept + dark_slope * hdr['EXPTIME']  # Reconstruct amplifier dark current map

recon_map = np.zeros_like(dark)  # Prepare full-size reconstructed dark map
for amp in range(namp):
    if amp % 2 == 0:
        flip = 1  # Even amplifiers: no flip
    else:
        flip = -1  # Odd amplifiers: flip
    # Fill in the amplifier region with the reconstructed amp map (flipped if needed)
    recon_map[:, amp * (dark.shape[1] // namp):(amp + 1) * (dark.shape[1] // namp)] = recon_amp[:, ::flip]

if doplot:
    # Plot the amplifier region, reconstructed amp, and their difference for visual inspection
    amp = 0  # Select amplifier 0 for plotting
    tmp = dark[:, amp * (dark.shape[1] // namp):(amp + 1) * (dark.shape[1] // namp)]
    fig, ax = plt.subplots(1, 3, figsize=(10, 8), sharex=True, sharey=True)
    ax[0].imshow(tmp, origin='lower', cmap='gray', vmin=-0.1, vmax=0.1, aspect='auto')
    ax[1].imshow(recon_amp, origin='lower', cmap='gray', vmin=-0.1, vmax=0.1, aspect='auto')
    ax[2].imshow(tmp - recon_amp, origin='lower', cmap='gray', vmin=-0.1, vmax=0.1, aspect='auto')
    plt.tight_layout()
    plt.show()

base_name_sample_dark = os.path.basename(dark_files[0])
# Save the corrected dark frame (dark minus reconstructed map) and the reconstructed map itself
fits.writeto('corr_' + base_name_sample_dark, dark - recon_map, header=hdr, overwrite=True)
fits.writeto('map_' + base_name_sample_dark, recon_map, header=hdr, overwrite=True)

if not os.path.exists(sample_dark_table):
    # Identify hot pixels: those with dark current above threshold and not NaN in the flat
    ypix, xpix = hot_pixels = np.where((dark > dark_threshold) & (np.isnan(flat) == False))

    sig = sigma(dark)  # Compute robust standard deviation of the dark frame
    nsig = dark[hot_pixels] / sig  # Normalize hot pixel values by standard deviation

    tbl = Table([nsig, xpix, ypix], names=('nsig', 'xpix', 'ypix'))  # Create table of hot pixel coordinates
    tbl.write(sample_dark_table, format='csv', overwrite=True)  # Save hot pixel table as CSV