import glob
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

def sigma(x):
    """Calculate the robust standard deviation of a dataset."""
    # Median Absolute Deviation (MAD) scaled to approximate std deviation for normal dist.
    return np.nanmedian(np.abs(x - np.nanmedian(x))) * 1.4826


# --- Configuration ---
# Threshold for flagging bad files (in units of sigma)
bad_nsig = 10.0

# Find all FITS files starting with 'NIRPS' in the current directory
files = glob.glob('NIRPS*.fits')

# List of FITS header keywords to extract
ccf_keys = ['RV_OBJ','CCFMFWHM']

# -- End of Configuration ---

# Dictionary to store extracted data
dict_keys = dict()

# Store file names
dict_keys['FILES'] = files

# Initialize empty lists for each header keyword
for key in ccf_keys:
    dict_keys[key] = []

# Loop over all files and extract header values
for file in files:
    # Read FITS header from extension 1
    hdr = fits.getheader(file,exten = 1)
    for key in ccf_keys:
        # If key exists in header, append its value; else, append NaN
        if key in hdr:
            dict_keys[key].append(hdr[key])
        else:
            dict_keys[key].append(np.nan)

# Convert lists to numpy arrays for easier math
for key in ccf_keys:
    dict_keys[key] = np.array(dict_keys[key])

# Compute median and robust std deviation for each parameter
median_rv = np.nanmedian(dict_keys['RV_OBJ']    )
median_fwhm = np.nanmedian(dict_keys['CCFMFWHM'])
sig_rv = sigma(dict_keys['RV_OBJ'])
sig_fwhm = sigma(dict_keys['CCFMFWHM'])

# Calculate the "distance" in sigma units for each file from the median
nsig = np.sqrt( ((dict_keys['RV_OBJ']-median_rv) / sig_rv)**2 + ((dict_keys['CCFMFWHM']-median_fwhm) / sig_fwhm)**2 )

# Identify files that are farther than bad_nsig from the median (outliers)
bad = nsig > bad_nsig

# Print number of bad files and their names
print(f'Number of bad files: {np.sum(bad)}')
bad_files = np.array(dict_keys['FILES'])[bad]

for bad_file in bad_files:
    print(f'Bad file: {bad_file}')

# Plotting the results
fig, ax = plt.subplots(figsize=(6, 10), nrows = 2, ncols = 1)

for i in range(2):
    # Plot good files (not flagged as bad) in blue
    ax[i].plot(dict_keys['RV_OBJ'][~bad], dict_keys['CCFMFWHM'][~bad], 'o', markersize=5,label = 'Good files')
    # Plot bad files (flagged as outliers) in red
    ax[i].plot(dict_keys['RV_OBJ'][bad], dict_keys['CCFMFWHM'][bad], 'o', markersize=5, color='red', label = 'Bad files')
    # Draw an ellipse representing the bad threshold region (centered at median, axes = 2*bad_nsig*sigma)
    ellipse = Ellipse((median_rv, median_fwhm), width=2*bad_nsig*sig_rv, height=2*bad_nsig*sig_fwhm, edgecolor='green', facecolor='none', linestyle='--',label='Bad threshold')
    ax[i].add_patch(ellipse)
    # Label axes
    ax[i].set_xlabel('RV_OBJ (km/s)')
    ax[i].set_ylabel('CCFMFWHM (km/s)')
    # Set plot title
    ax[i].set_title('CCF MFWHM vs RV_OBJ')
    # Add grid for readability
    ax[i].grid()

# For the second subplot, zoom in around the median ± 1.5*bad_nsig*sigma
ax[1].set_xlim([-1.5*bad_nsig*sig_rv + median_rv, 1.5*bad_nsig*sig_rv + median_rv])
ax[1].set_ylim([-1.5*bad_nsig*sig_fwhm + median_fwhm, 1.5*bad_nsig*sig_fwhm + median_fwhm])
# Add legend to the second subplot
ax[1].legend(loc='upper right')

# Adjust layout to prevent overlap
plt.tight_layout()
# Save the figure to a PNG file
plt.savefig('qc_ccffwhm.png', dpi=300)
# Display the plot on screen
plt.show()