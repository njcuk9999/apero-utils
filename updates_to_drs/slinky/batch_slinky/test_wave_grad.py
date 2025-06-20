"""
This script processes FITS files to extract and plot the gradient of the wavelength data.
The color of each plot is determined by the Modified Julian Date (MJD) of the observation.
A colorbar is added to represent the MJD values.
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob
from tqdm import tqdm
import matplotlib.cm as cm
import matplotlib.colorbar as colorbar
from etienne_tools import read_pickle, save_pickle
import os

# Load FITS files matching the specified pattern
files = np.array(glob.glob('*_pp_e2dsff_A_wave_night_A.fits'))
files = files[0:200]

# Extract MJD values with a progress bar
mjdmid = []
for i in tqdm(range(len(files)), desc="Extracting MJD values"):

    pickle_name = files[i].replace('.fits','.pkl')

    if not os.path.isfile(pickle_name):
        dict_data = {}

        hdr = fits.getheader(files[i])
        wave = fits.getdata(files[i])

        dict_data['hdr'] = hdr
        dict_data['wave'] = wave

        save_pickle(pickle_name,dict_data)
    else:
        dict_data = read_pickle(pickle_name)
        hdr = dict_data['hdr']
        wave = dict_data['wave']


    mjdmid.append(hdr['MJDMID'])
    if i == 0:
        wave0 = np.zeros([len(files),wave.shape[0],wave.shape[1]])
    wave0[i] = wave

mjdmid = np.array(mjdmid)

ord = np.argsort(mjdmid)
wave0 = wave0[ord]
files = files[ord]
mjdmid = mjdmid[ord]

# Normalize MJD values for color mapping
norm = plt.Normalize(vmin=np.min(mjdmid), vmax=np.max(mjdmid))
cmap = cm.inferno

# Create a plot
fig, ax = plt.subplots()

pix1 = 1024
pix2 = 3072
step = 512

wave_mean = np.nanmean(wave0,axis=0)
# For the first file, store the wavelength and gradient data for normalization
grad = np.gradient(wave_mean, axis=1) / wave_mean
grad = grad[:, pix1:pix2:step]
wave_mean = wave_mean[:, pix1:pix2:step]
grad0 = np.array(grad)

wave_mean[wave_mean > 1800] = np.nan
wave_mean[(wave_mean > 1350) * (wave_mean < 1450)] = np.nan


# Loop through each file with a progress bar
for i, file in enumerate(tqdm(files, desc="Processing files")):
    # Read the wavelength data from the FITS file
    wave = wave0[i]

    # Calculate the gradient of the wavelength data
    grad = np.gradient(wave, axis=1) / wave
    grad = grad[:, pix1:pix2:step]

    # Determine the color for the current file based on its MJD value
    color = cmap(norm(mjdmid[i]))



    # Plot the normalized gradient data
    for ord in range(wave_mean.shape[0]):
        ax.plot(wave_mean[ord], (grad/grad0)[ord],'.-' , label=file, alpha=0.2, color=color)

# Add a colorbar to the plot to represent the MJD values
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('MJD')

# Display the plot
plt.show()