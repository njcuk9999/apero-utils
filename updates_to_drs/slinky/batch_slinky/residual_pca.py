import numpy as np
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from etienne_tools import lowpassfilter, doppler, save_pickle, read_pickle
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from tqdm import tqdm
from wpca import WPCA, EMPCA
from astropy.table import Table
import os
import fitsio
from scipy.optimize import minimize
import warnings
from scipy.special import erf

# Suppress the specific RuntimeWarning for all-NaN slices
warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN slice encountered")

key_vsys = 'BERV'

import fitsio

global temporary_files
temporary_files = []

def whoami():
    """
    Find the user name of the session
    """
    import os
    return os.popen('whoami').read().strip()

def fiber_setup(files):
    dd = read_t(files[0])
    if 'FluxAB' in dd:
        return 'AB'
    if 'FluxA' in dd:
        return 'A'
    else:
        return 'X'

def write_t(data, file, temporary=True):
    """
    Write a dictionary of data to a Multi-Extension FITS (MEF) file.

    This function opens a FITS file for reading and writing ('rw' mode), and allows overwriting if the file already exists.
    It iterates over each key in the data dictionary, retrieves the data and its associated header, and writes them
    to the FITS file as new extensions.

    Parameters:
    data (dict): A dictionary containing the data and headers to be written to the FITS file.
    file (str): The path to the FITS file to be written.

    Returns:
    None
    """

    if temporary:
        pickle_name = file.replace('.fits', '.pkl')
        temporary_files.append(pickle_name)


    # Open the FITS file for reading and writing ('rw' mode), and allow overwriting if the file already exists (clobber=True).
    with fitsio.FITS(file, 'rw', clobber=True) as outfile:
        # Iterate over each key in the data dictionary.
        for key in data.keys():
            # Check if the key does not contain '_header' to avoid processing header keys.
            if '_header' not in key:
                # Retrieve the data associated with the current key.
                data_to_write = data[key]
                
                if key+'_header' not in data:
                    header = fitsio.FITSHDR()
                else:
                    header = data[key + '_header']
                
                # Set the 'EXTNAME' keyword in the header to the current key. This names the extension in the FITS file.
                header['EXTNAME'] = key
                
                # Write the data and its header to the FITS file as a new extension.
                # The extname parameter explicitly sets the extension name in the FITS file.
                outfile.write(data_to_write, header=header, extname=key)
    
    if temporary:
        save_pickle(pickle_name,data)

def read_t(file,keys = None, temporary=True):
    """
    Read a Multi-Extension FITS (MEF) file and create a dictionary containing all the extensions and their headers.

    This function opens a FITS file for reading and iterates over each Header Data Unit (HDU) in the file.
    It reads the data and headers from each HDU and stores them in a dictionary, with the extension names as keys.

    Parameters:
    file (str): The path to the FITS file to be read.

    Returns:
    dict: A dictionary containing the data and headers from all the extensions in the FITS file.
    """
    if temporary:
        pickle_name = file.replace('.fits', '.pkl')
        if os.path.exists(pickle_name):
            return read_pickle(pickle_name)
        keys = None


    # Open the FITS file for reading.
    with fitsio.FITS(file) as infile:
        # Initialize an empty dictionary to store the data and headers.
        data = dict()
        
        # Iterate over each Header Data Unit (HDU) in the FITS file.
        for hdu in infile:
            # Get the extension name of the current HDU.
            key = hdu.get_extname()
            if (keys is not None) and (key not in keys):
                continue
            
            # Read the data from the current HDU and store it in the dictionary with the extension name as the key.
            data[key] = hdu.read()
            
            # Read the header from the current HDU and store it in the dictionary with the key appended by '_header'.
            data[key + '_header'] = hdu.read_header()
    
    if temporary:
        temporary_files.append(pickle_name)
        save_pickle(pickle_name,data)

    # Return the dictionary containing all the data and headers.
    return data

def sigma(x):
    """
    Calculate the standard deviation (sigma) of the input array x.

    This function calculates the standard deviation by finding the difference between the percentiles
    corresponding to -1 and +1 sigma cuts, which are calculated using the error function (erf), and then dividing by 2.

    Parameters:
    x (array-like): Input array for which the standard deviation is to be calculated.

    Returns:
    float: The calculated standard deviation (sigma) of the input array.
    """
    
    # Calculate the percentile corresponding to +1 sigma using the error function.
    upper_percentile = 100 * (1 + erf(1 / np.sqrt(2))) / 2
    
    # Calculate the percentile corresponding to -1 sigma using the error function.
    lower_percentile = 100 * (1 + erf(-1 / np.sqrt(2))) / 2
    
    # Calculate the upper value at the +1 sigma percentile.
    upper_value = np.nanpercentile(x, upper_percentile)
    
    # Calculate the lower value at the -1 sigma percentile.
    lower_value = np.nanpercentile(x, lower_percentile)
    
    # Calculate the standard deviation (sigma) by finding the difference between the upper and lower values,
    # and then dividing by 2.
    return (upper_value - lower_value) / 2

def mkbins(p1, p2, n1, n2):
    """
    Create a 2D binning of data based on two parameters.

    This function sorts the data based on the first parameter and assigns bins accordingly.
    It then further sorts the data within each bin based on the second parameter and assigns sub-bins.
    The result is a 2D binning of the data.

    Parameters:
    p1 (array-like): The first parameter used for binning.
    p2 (array-like): The second parameter used for sub-binning within each bin.
    n1 (int): The number of bins for the first parameter.
    n2 (int): The number of sub-bins for the second parameter within each bin.

    Returns:
    array-like: An array of bin indices representing the 2D binning of the data.
    """
    
    # Sort the data based on the first parameter and assign bins accordingly.
    bin1 = np.array(np.argsort(np.argsort(p1)) / len(p1) * n1, dtype=int)
    
    # Initialize an array to store the sub-bin indices.
    bin2 = np.zeros_like(bin1)

    # Iterate over each unique bin index in bin1.
    for ubin1 in np.unique(bin1):
        # Create a boolean mask for the current bin.
        g = bin1 == ubin1
        
        # Sort the data within the current bin based on the second parameter and assign sub-bins.
        bin2[g] = np.argsort(np.argsort(p2[g])) / np.sum(g) * n2

    # Combine the bin indices from bin1 and bin2 to create a 2D binning.
    bin12 = bin1 * n2 + bin2

    # Return the array of bin indices representing the 2D binning of the data.
    return bin12

def plot_map(map1, map2 = None, map3=None):

    nplots = 1
    if map2 is not None:
        nplots += 1
    if map3 is not None:
        nplots += 1


    if nplots == 1:

        v1 = np.nanpercentile(map1, 1)
        v2 = np.nanpercentile(map1, 99)
        rangeplot = (v2+v1)/2 + np.array([-1.5, 1.5]) * (v2-v1) / 2

        fig, ax = plt.subplots(1, 1, figsize=[10, 5])
        im = ax.imshow(map1, aspect='auto', vmin=rangeplot[0], vmax=rangeplot[1], interpolation='nearest')
        fig.colorbar(im, ax=ax)
        plt.show()
    else:
        fig, ax = plt.subplots(nplots,1,  figsize=[10*nplots, 5], sharex='all', sharey='all')

        v1 = np.nanpercentile(map1, 1)
        v2 = np.nanpercentile(map1, 99)
        rangeplot = (v2+v1)/2 + np.array([-1.5, 1.5]) * (v2-v1) / 2


        im = ax[0].imshow(map1, aspect='auto', vmin=rangeplot[0], vmax=rangeplot[1], interpolation='nearest')
        fig.colorbar(im, ax=ax[0])
        if map2 is not None:
            v1 = np.nanpercentile(map2, 1)
            v2 = np.nanpercentile(map2, 99)
            rangeplot = (v2+v1)/2 + np.array([-1.5, 1.5]) * (v2-v1) / 2


            im = ax[1].imshow(map2, aspect='auto', vmin=rangeplot[0], vmax=rangeplot[1], interpolation='nearest')
            fig.colorbar(im, ax=ax[1])
        if map3 is not None:
            v1 = np.nanpercentile(map3, 1)
            v2 = np.nanpercentile(map3, 99)
            rangeplot = (v2+v1)/2 + np.array([-1.5, 1.5]) * (v2-v1) / 2

            im = ax[2].imshow(map3, aspect='auto', vmin=rangeplot[0], vmax=rangeplot[1], interpolation='nearest')
            fig.colorbar(im, ax=ax[2])
        plt.show()

def mk_residual_maps(files, obj, outpath, doplot=False):
    """
    Create residual maps for a given object using PCA correction.

    This function processes a list of FITS files, extracts relevant header information, and creates residual maps
    by applying PCA correction. The resulting maps are saved as pickle files.

    Parameters:
    files (list): List of paths to the FITS files to be processed.
    obj (str): The object name to be used in the output file names.
    outpath (str): The path to the directory where the output files will be saved.

    Returns:
    None
    """
    
    # Initialize arrays to store header information.
    bervs = np.zeros(len(files))
    mjds = np.zeros(len(files))
    TLPEH2O = np.zeros(len(files))
    TLPEOTR = np.zeros(len(files))  
    fiber = fiber_setup(files)
    
    # Extract header information from each file.
    for i in tqdm(range(len(files))):
        #hdr = fits.getheader(files[i],ext=1)
        dd = read_t(files[i])
        hdr = dd[f'Flux{fiber}_header']
        bervs[i] = hdr[key_vsys]
        mjds[i] = hdr['MJD-OBS']
        TLPEH2O[i] = hdr['TLPEH2O']
        TLPEOTR[i] = hdr['TLPEOTR']

    # Determine the number of bins based on the number of files.
    if len(files) > 50:
        nbin = 5, 5
    else:
        nbin = 2, 2

    # Create 2D binning based on TLPEH2O and TLPEOTR.
    bin_tellu = mkbins(TLPEH2O, TLPEOTR, nbin[0], nbin[1])
    ubin_tellu = np.unique(bin_tellu)

    if doplot:
        for ubin_tellu in np.unique(bin_tellu):
            gg = bin_tellu == ubin_tellu
            plt.plot(TLPEH2O[gg], TLPEOTR[gg], 'o')
        plt.xlabel('TLPEH2O')
        plt.ylabel('TLPEOTR')
        plt.show()

    # Round MJD values to the nearest integer and find unique values.
    jd_round = np.array(mjds, dtype=int)

    # Get the reference wavelength from the first file.
    dd = read_t(files[0])
    waveref = dd[f'Wave{fiber}']

    # Loop over each order in the wavelength reference.
    for iord in range(0, waveref.shape[0]):

        outname = f'{outpath}/residual_pca_tcorr_{obj}_{iord}.fits'
        print(outname)


        grad_log_wave = np.gradient(np.log(waveref[iord]))

        # Check if the output file already exists.
        if not os.path.isfile(outname):
            # Initialize arrays to store the 2D map and other data.
            map2d = np.zeros([len(files), 4088])
            sp_amp = np.zeros_like(map2d)

            # Process each file.
            for i in tqdm(range(len(files)), leave=False):
                # Get the wavelength and flux data for the current order.
                dd = read_t(files[i], keys = [f'Wave{fiber}', f'Flux{fiber}'])
                wave = dd[f'Wave{fiber}'][iord]
                tmp = dd[f'Flux{fiber}'][iord]

                # Calculate the amplitude of the spectrum.
                with warnings.catch_warnings():
                    sp_amp[i] = np.nanmedian(tmp)

                # Check for finite values and interpolate.
                g = np.isfinite(tmp)
                if np.sum(g) > 20:
                    # Apply Doppler shift correction and interpolate the data.
                    tmp2 = ius(doppler(wave[g], -1e3 * bervs[i]), tmp[g], k=1, ext=1)(wave)
                    mask = ius(doppler(wave, -1e3 * bervs[i]), np.array(g, dtype=float), k=1, ext=1)(wave) > 0.9

                    # Mask invalid values.
                    tmp2[tmp2 == 0] = np.nan
                    tmp2[~mask] = np.nan
                    map2d[i] = tmp2

            recon2 = np.zeros_like(map2d.copy(), dtype = float)
            Nite = 3
            for ite in range(Nite):
                # Initialize the residual map.
                with warnings.catch_warnings():
                    sp_star = np.nanmedian((map2d - recon2)/sp_amp, axis=0)
                    grad_velo_sp_star = np.gradient(sp_star)/grad_log_wave
                    
                
                # Calculate the residual map.
                residual_map = []
                for i in range(len(files)):
                    tmp = (map2d[i] - sp_star*sp_amp[i])
                    g = np.isfinite(tmp)

                    amp_grad = np.nansum(grad_velo_sp_star[g] * tmp[g]) / np.nansum(grad_velo_sp_star[g]**2)
                    tmp -= amp_grad * grad_velo_sp_star
                    residual_map.append(tmp)
                    #plt.plot(mjds[i], amp_grad, 'o')
                #plt.show()

                residual_map = np.array(residual_map)
                
                # Interpolate and filter the residual map.
                for i in range(len(files)):
                    tmp = residual_map[i]
                    g = np.isfinite(tmp)
                    if np.sum(g) > 20:
                        # Apply Doppler shift correction and interpolate the data.
                        tmp2 = ius(doppler(wave[g], 1e3 * bervs[i]), tmp[g], k=1, ext=1)(waveref[iord])
                        mask = ius(doppler(wave, 1e3 * bervs[i]), np.array(g, dtype=float), k=1, ext=1)(waveref[iord]) > 0.99

                        # Mask invalid values and apply low-pass filter.
                        tmp2[tmp2 == 0] = np.nan
                        tmp2[~mask] = np.nan
                        tmp2 -= lowpassfilter(tmp2, 51)
                        residual_map[i] = tmp2

                # Bin the residual map.
                residual_map_bin = np.zeros([len(np.unique(bin_tellu)), 4088])
                for i, ubin_tellu in enumerate(np.unique(bin_tellu)):
                    g = bin_tellu == ubin_tellu
                    with warnings.catch_warnings():
                        residual_map_bin[i] = np.nanmedian(residual_map[g], axis=0)



                # Construct a PCA basis and the corresponding weights for bins
                sigmap_bin = np.ones_like(residual_map_bin)
                print(f'We determine the PCA weights for bin map {ite}/{Nite}')
                for i in range(residual_map_bin.shape[0]):
                    sigmap_bin[i] = sigma(residual_map_bin[i])

                weights_bin = 1 / sigmap_bin
                bad_bin = ~np.isfinite(residual_map_bin + weights_bin)

                residual_map_bin2 = residual_map_bin.copy()
                residual_map_bin2[bad_bin] = 0
                weights_bin2 = weights_bin.copy()
                weights_bin2[bad_bin] = 0

                if (ite != Nite-1):
                    # Construct a PCA basis and the corresponding weights.
                    sigmap = np.ones_like(residual_map)
                    print(f'We determine the PCA weights full pixel map {ite}/{Nite}')
                    for i in range(residual_map_bin.shape[0]):
                        sigmap[i] = sigma(residual_map[i])

                    weights = 1 / sigmap
                    bad = ~np.isfinite(residual_map + weights)

                    residual_map2 = residual_map.copy()
                    residual_map2[bad] = 0
                    weights[bad] = 0

                    wpca = EMPCA(n_components=2)
                    print(f'Fitting PCA {ite}/{Nite} on binned data')
                    wpca.fit(residual_map_bin2, weights_bin2)
                    print(f'Reconstructing full pixel map {ite}/{Nite}')
                    try:
                        recon = wpca.reconstruct(residual_map2, weights)
                    except:
                        recon = np.zeros_like(residual_map)

                    # unshift the recon back to the input grid
                    recon2 = np.zeros_like(residual_map)
                    for i in range(len(files)):
                        tmp = recon[i]
                        g = np.isfinite(tmp)
                        if np.sum(g) > 20:
                            tmp2 = ius(waveref[iord], tmp[g], k=1, ext=1)(doppler(waveref[iord], 1e3 * bervs[i]))
                            mask = ius(waveref[iord], np.array(g, dtype=float), k=1, ext=1)(doppler(waveref[iord], 1e3 * bervs[i])) > 0.9
                            tmp2[tmp2 == 0] = np.nan
                            tmp2[~mask] = np.nan
                            recon2[i] = tmp2

            # Mask invalid values in the residual map and weights.
            print('We construct the PCA basis')

            # Create a dictionary to store the residuals and other data.
            dict_residuals = dict()
            dict_residuals['wave'] = waveref[iord]
            dict_residuals['residual_map'] = residual_map_bin2
            dict_residuals['weights'] = weights_bin2
            dict_residuals['sp_star'] = sp_star

            # Save the dictionary as a pickle file.
            write_t(dict_residuals, outname)


def apply_residual_corr(files, obj, outpath, replace_path):
    """
    Apply residual correction to a list of FITS files using PCA.

    This function processes a list of FITS files, applies residual correction using PCA, and saves the corrected
    spectra to new FITS files.

    Parameters:
    files (list): List of paths to the FITS files to be processed.
    obj (str): The object name to be used in the output file names.
    outpath (str): The path to the directory where the output files will be saved.

    Returns:
    None
    """
    
    # Find the number of orders by getting the wavelength reference from the first file.
    dd = read_t(files[0])
    waveref = dd['WaveA']
    fiber = fiber_setup(files)

    # Initialize a dictionary to store all residual tables.
    all_residual_tables = dict()

    # Loop over each order in the wavelength reference.
    for iord in range(0, waveref.shape[0]):
        outname = f'{outpath}/residual_pca_tcorr_{obj}_{iord}.fits'
        
        # Check if the PCA correction file exists.
        if not os.path.isfile(outname):
            print('File not found', outname)
            continue

        # Read the PCA correction data from the pickle file.
        dict_pca = read_t(outname)
        residual_map_bin = dict_pca['residual_map']
        weights = dict_pca['weights']
        
        # Perform weighted PCA on the residual map.
        wpca = EMPCA(n_components=2)
        wpca.fit(residual_map_bin, weights)

        print(f'Loading {outname} for spline')
        wave = dict_pca['wave']
        sp_star = dict_pca['sp_star']

        # Store the PCA and spline data in the dictionary.
        all_residual_tables[iord] = dict()
        all_residual_tables[iord]['wave'] = wave
        g = np.isfinite(sp_star)
        if np.sum(g) > 20:
            all_residual_tables[iord]['spl_star'] = ius(wave[g], sp_star[g], k=3, ext=1)
        else:
            all_residual_tables[iord]['spl_star'] = ius(wave, np.zeros_like(wave), k=3, ext=1)
        all_residual_tables[iord]['wpca'] = wpca

    # Loop over each file to apply the residual correction.
    for i in range(len(files)):
        outname_full_correction = files[i].replace(replace_path[0], replace_path[1])
        outname_slinky = files[i].replace(replace_path[0], replace_path[2])
        
        # get folder name for outname
        folder = os.path.dirname(outname_full_correction)
        if not os.path.exists(folder):
            print(f'Creating folder {folder}')
            os.makedirs(folder)
        # same for slinky
        folder = os.path.dirname(outname_slinky)
        if not os.path.exists(folder):
            print(f'Creating folder {folder}')
            os.makedirs(folder)


        # Skip if the output file already exists.
        if os.path.isfile(outname):
            continue
        
        # Read the flux and wavelength data from the file.

        dd = read_t(files[i])
        sp = dd[f'Flux{fiber}']
        wave = dd[f'Wave{fiber}']
        hdr = dd[f'Flux{fiber}_header']
        berv = hdr[key_vsys]

        # fetch the slinky-corrected wavelength file
        wavelength_calib = dd[f'Wave{fiber}_header']['WFP_FILE']

        # Loop over each order in the wavelength reference.
        for iord in tqdm(range(waveref.shape[0]), leave=False):
            if iord not in all_residual_tables:
                continue

            # Apply Doppler shift correction to the star spectrum.
            sp_star = all_residual_tables[iord]['spl_star'](doppler(wave[iord], -1e3 * berv))
            sp_star[sp_star == 0] = np.nan

            # Check for finite values in the spectrum and star spectrum.
            valid = np.isfinite(sp[iord]) & np.isfinite(sp_star)
            if np.sum(valid) < 20:
                continue

            # Calculate the median ratio of the spectrum to the star spectrum.
            med = np.nanmedian(sp[iord] / sp_star)

            # Calculate the residuals.
            residu1 = sp[iord] - med * sp_star
            
            # Interpolate the residuals to the PCA wavelength grid.
            g = np.isfinite(residu1)
            residu2 = ius(wave[iord][g], residu1[g], k=1, ext=1)(all_residual_tables[iord]['wave'])
            weights2 = np.ones_like(residu2)
            bad = ~np.isfinite(residu2) | (residu2 == 0)
            residu2[bad] = 0
            weights2[bad] = 0
            residu2 = residu2.reshape(1, len(residu2))
            weights2 = weights2.reshape(1, len(weights2))
            
            # Reconstruct the residuals using PCA.
            recon = all_residual_tables[iord]['wpca'].reconstruct(residu2, weights2).flatten()
            
            # Interpolate the reconstructed residuals back to the original wavelength grid.
            recon = ius(all_residual_tables[iord]['wave'], recon, k=1, ext=1)(wave[iord])

            # Optionally plot the residuals and reconstruction for debugging.
            if iord == -1:
                fig, ax = plt.subplots(2, 1, figsize=[10, 5], sharex=True)
                ax[0].plot(residu1, alpha=0.9)
                ax[0].plot(recon, alpha=0.5)
                ax[1].plot(residu1 - recon)
                plt.show()

            # Subtract the reconstructed residuals from the spectrum.
            sp[iord] -= recon

        # Read the original FITS file data.
        dd = read_t(files[i])
        
        # Update the flux data with the corrected spectrum.
        dd[f'Flux{fiber}'] = sp
        
        # Write the corrected data to a new FITS file.
        print(f'[{i + 1}/{len(files)}] Writing corrected file: {outname_full_correction}')
        write_t(dd, outname_full_correction)


def mk_residual_plot(files, obj, iord=60):
    """
    Create a residual plot for a given object and order.

    This function processes a list of FITS files, extracts relevant header information, and creates a residual plot
    for a specified order by applying PCA correction.

    Parameters:
    files (list): List of paths to the FITS files to be processed.
    obj (str): The object name to be used in the plot title.
    iord (int): The order to be plotted.

    Returns:
    None
    """
    
    # Initialize an array to store BERV values.
    bervs = np.zeros(len(files))
    fiber = fiber_setup(files)
    

    # Extract BERV values from each file's header.
    for i in tqdm(range(len(files))):
        dd = read_t(files[i])
        hdr = dd[f'Flux{fiber}_header']
        bervs[i] = hdr[key_vsys]
    
    # Sort the files based on BERV values.
    ord = np.argsort(bervs)
    bervs = bervs[ord]
    files = np.array(files)[ord]

    # Get the reference wavelength from the first file.
    waveref = read_t(files[0])['WaveA']

    # Initialize arrays to store the 2D map and other data.
    map2d = np.zeros([len(files), 4088])
    sp_star = np.ones(waveref.shape[1], dtype=float)
    sp_amp = np.zeros(len(files))

    # Process each file.
    for i in tqdm(range(len(files)), leave=False):

        dd = read_t(files[i])
        # Get the wavelength and flux data for the current order.
        wave = dd[f'Wave{fiber}'][iord]
        tmp = dd[f'Flux{fiber}'][iord]

        # Calculate the amplitude of the spectrum.
        with warnings.catch_warnings():
            sp_amp[i] = np.nanmedian(tmp / sp_star)
        tmp = tmp / sp_amp[i]

        # Check for finite values and interpolate.
        g = np.isfinite(tmp)
        if np.sum(g) > 20:
            # Apply Doppler shift correction and interpolate the data.
            tmp2 = ius(doppler(wave[g], -1e3 * bervs[i]), tmp[g], k=1, ext=1)(wave)
            mask = ius(doppler(wave, -1e3 * bervs[i]), np.array(g, dtype=float), k=1, ext=1)(wave) > 0.9

            # Mask invalid values.
            tmp2[tmp2 == 0] = np.nan
            tmp2[~mask] = np.nan
            map2d[i] = tmp2

    # Initialize the residual map.
    residual_map = np.zeros_like(map2d)
    with warnings.catch_warnings():
        sp_star = np.nanmedian(map2d, axis=0)
    
    # Calculate the residual map.
    for i in range(len(files)):
        residual_map[i] = (map2d[i] - sp_star) * sp_amp[i]
    
    # Interpolate and filter the residual map.
    for i in range(len(files)):
        tmp = residual_map[i]
        g = np.isfinite(tmp)
        if np.sum(g) > 20:
            # Apply Doppler shift correction and interpolate the data.
            tmp2 = ius(doppler(wave[g], 1e3 * bervs[i]), tmp[g], k=1, ext=1)(waveref[iord])
            mask = ius(doppler(wave, 1e3 * bervs[i]), np.array(g, dtype=float), k=1, ext=1)(waveref[iord]) > 0.9

            # Mask invalid values and apply low-pass filter.
            tmp2[tmp2 == 0] = np.nan
            tmp2[~mask] = np.nan
            tmp2 -= lowpassfilter(tmp2, 51)
            residual_map[i] = tmp2
    
    # Determine the minimum and maximum values for the color scale.
    vmin = np.nanpercentile(residual_map, 1)
    vmax = np.nanpercentile(residual_map, 99)

    # Create the residual plot.
    fig, ax = plt.subplots(1, 1, figsize=[10, 5])
    im = ax.imshow(residual_map, aspect='auto', vmin=vmin, vmax=vmax, interpolation='nearest')
    fig.colorbar(im, ax=ax)
    plt.savefig(f'{obj}_residuals_{iord}.png')
    if whoami() == 'eartigau':
        plt.show()
    plt.close()



def residual_pca(files, outpath, replace_path):
    if not files:
        print("No files provided for PCA.")
        return
    fiber = fiber_setup(files)
    obj = read_t(files[0])[f'Flux{fiber}_header']['DRSOBJN']

    mk_residual_maps(files, obj, outpath)
    apply_residual_corr(files, obj, outpath, replace_path)


#objs = ['TOI3397','TOI2952','TOI210','TOI4552','TOI1078','TOI4666']


def run_batch(batch_id = None, params = None):

    objs = params['object_of_interest']
    outpath = params['pca_mef_dir']
    
    replace_path = params['batchname'], params['batchname_corr'], params['batchname_slinky']

    file_path = params['search_tcorr_path'].replace('science/*/','science/{}/')

    for obj in objs:
        files = glob.glob(file_path.format(obj))
        if files:
            residual_pca(files, outpath, replace_path)
        else:
            print(f"No files found for object: {obj}")

