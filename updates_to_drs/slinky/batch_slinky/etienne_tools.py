import glob
import os
import pickle
import warnings
from datetime import datetime
from io import StringIO

import astropy.units as uu
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
import requests
import statsmodels.api as sm
from astropy import constants
from astropy.io import fits
from astropy.table import Table
from astropy.table import vstack
from astroquery.simbad import Simbad
from numba import jit
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.optimize import curve_fit
from tqdm import tqdm
from typing import Optional
from scipy.constants import c


def rot_broad(wvl: np.ndarray, flux: np.ndarray, epsilon: float, vsini: float,
              eff_wvl: Optional[float] = None) -> np.ndarray:
    """
    **********************************************************************
    ***** THIS FUNCTION IS COPIED FROM PyAstronomy/pyasl/rotBroad.py *****
    ***** AND MODIFIED TO USE THE SAME CONVENTIONS AS THE LBL CODE.  *****
    ***** AND AVOID DEPENDENCIES ON OTHER PyAstronomy MODULES.       *****
    **********************************************************************

    Apply rotational broadening using a single broadening kernel.
    The effect of rotational broadening on the spectrum is
    wavelength dependent, because the Doppler shift depends
    on wavelength. This function neglects this dependence, which
    is weak if the wavelength range is not too large.
    .. note:: numpy.convolve is used to carry out the convolution
              and "mode = same" is used. Therefore, the output
              will be of the same size as the input, but it
              will show edge effects.
    Parameters
    ----------
    wvl : array
        The wavelength
    flux : array
        The flux
    epsilon : float
        Linear limb-darkening coefficient
    vsini : float
        Projected rotational velocity in km/s.
    eff_wvl : float, optional
        The wavelength at which the broadening
        kernel is evaluated. If not specified,
        the mean wavelength of the input will be
        used.
    Returns
    -------
    Broadened spectrum : array
        The rotationally broadened output spectrum.
    """
    # Wavelength binsize
    dwl = wvl[1] - wvl[0]
    # deal with no effective wavelength
    if eff_wvl is None:
        eff_wvl = np.mean(wvl)
    # The number of bins needed to create the broadening kernel
    binn_half = int(np.floor(((vsini / (c / 1000)) * eff_wvl / dwl))) + 1
    gwvl = (np.arange(4 * binn_half) - 2 * binn_half) * dwl + eff_wvl
    # Create the broadening kernel
    dl = gwvl - eff_wvl
    # -------------------------------------------------------------------------
    # this bit is from _Gdl.gdl
    #    Calculates the broadening profile.
    # -------------------------------------------------------------------------
    # set vc
    vc = vsini / (c / 1000)
    # set eps (make sure it is a float)
    eps = float(epsilon)
    # calculate the max vc
    dlmax = vc * eff_wvl
    # generate the c1 and c2 parameters
    c1 = 2 * (1 - eps) / (np.pi * dlmax * (1 - eps / 3))
    c2 = eps / (2 * dlmax * (1 - eps / 3))
    # storage for the output
    bprof = np.zeros(len(dl))
    # Calculate the broadening profile
    xvec = dl / dlmax
    indi0 = np.where(np.abs(xvec) < 1.0)[0]
    bprof[indi0] = c1 * np.sqrt(1 - xvec[indi0] ** 2) + c2 * (1 - xvec[indi0] ** 2)
    # Correct the normalization for numeric accuracy
    # The integral of the function is normalized, however, especially in the
    # case of mild broadening (compared to the wavelength resolution), the
    # discrete  broadening profile may no longer be normalized, which leads to
    # a shift of the output spectrum, if not accounted for.
    bprof /= (np.sum(bprof) * dwl)
    # -------------------------------------------------------------------------
    # Remove the zero entries
    indi = np.where(bprof > 0.0)[0]
    bprof = bprof[indi]
    # -------------------------------------------------------------------------
    result = np.convolve(flux, bprof, mode="same") * dwl
    return result


def latex_quantities(quantity_name, value, filename):
    # if a number is present in the 'quantity_name', we should have an error raised
    if any(char.isdigit() for char in quantity_name):
        raise ValueError('quantity_name should not contain any number')

    # if file does not exist, we create it
    if not os.path.exists(filename):
        with open(filename, 'w') as file:
            file.write('')

    # first we read the full file
    with open(filename, 'r') as file:
        data = file.readlines()
    # then we remove the line with the quantity if it exists.
    for i in range(len(data)):
        # \def\nstarsharps{50}
        if '\\def\\' + quantity_name + '{' in data[i]:
            data.pop(i)
            break
    # then we add the new line
    data.append('\\def\\' + quantity_name + '{' + str(value) + '}\n')
    # then we write the new file
    with open(filename, 'w') as file:
        file.writelines(data)


def mjd_to_matplotlib_date(mjd):
    return mjd - 40587.50


def gaucho_table(no_update=False):
    MAIN_URL = ('https://docs.google.com/spreadsheets/d/'
                '13PDmoAirm6p6zL1WcC_JMM6NqQD-C4cq2VBL4whMWIQ/'
                'export?format=csv&gid=0')

    # open main table
    main_dataframe = pd.read_csv(StringIO(requests.get(MAIN_URL).text))

    # if you prefer an astropy table
    gaucho_table = Table.from_pandas(main_dataframe)

    if not no_update:
        date = datetime.isoformat(datetime.now()).split('T')[0]
        nirps_file = '/Users/eartigau/.hidden_files/nirps_object_table_{}.csv'.format(date)
        if not os.path.exists(nirps_file):
            os.system(
                'scp nirps-client@maestria:/cosmos99/nirps/apero-data/nirps_he_online/other/ari/nirps_he_online_udem'
                '/object_table.csv {}'.format(nirps_file))
        spirou_file = '/Users/eartigau/.hidden_files/spirou_object_table_{}.csv'.format(date)
        if not os.path.exists(spirou_file):
            os.system(
                'scp nirps-client@maestria:/cosmos99/spirou/apero-data/spirou_offline/other/ari/spirou_offline_udem'
                '/object_table.csv {}'.format(spirou_file))
    else:
        nirps_files = np.array(glob.glob('/Users/eartigau/.hidden_files/nirps_object_table_*.csv'))
        nirps_files = nirps_files[np.argsort(nirps_files)]
        nirps_file = nirps_files[-1]
        spirou_files = np.array(glob.glob('/Users/eartigau/.hidden_files/spirou_object_table_*.csv'))
        spirou_files = spirou_files[np.argsort(spirou_files)]
        spirou_file = spirou_files[-1]

    tbl_spirou = Table.read(spirou_file)
    tbl_nirps = Table.read(nirps_file)
    tbl_planet = Table.read('/Users/eartigau/.hidden_files/exoplanets.csv')

    keys_float = ['Distance pc', 'MASS', 'H mag', 'RA', 'DEC']
    for key in keys_float:
        v = gaucho_table[key]
        v = np.array(v, dtype='U999')
        v[v == ''] = 'nan'
        v = [vv.replace(',', '.') for vv in v]
        gaucho_table[key] = np.array(v, dtype=float)

    for i in range(len(gaucho_table)):
        ra = gaucho_table['RA'][i]
        dec = gaucho_table['DEC'][i]

        rad_spirou = np.sqrt(
            (tbl_spirou['RA [Deg]'] - ra) ** 2 / np.cos(dec / (180 / np.pi)) + (tbl_spirou['Dec [Deg]'] - dec) ** 2)
        if np.min(rad_spirou) < 0.01:
            imin = np.argmin(rad_spirou)
            if gaucho_table[i]['NVISIT_SPIROU'] != tbl_spirou[imin]['raw files']:
                print('Matched {} with {}'.format(gaucho_table['APERO NAME'][i], tbl_spirou['OBJNAME'][imin]))
                print('Nvisit SPIROU --> {}'.format(tbl_spirou[imin]['raw files']))

        rad_nirps = np.sqrt(
            (tbl_nirps['RA [Deg]'] - ra) ** 2 / np.cos(dec / (180 / np.pi)) + (tbl_nirps['Dec [Deg]'] - dec) ** 2)
        if np.min(rad_nirps) < 0.01:
            imin = np.argmin(rad_nirps)
            if gaucho_table[i]['NVISIT_NIRPS'] != tbl_nirps[imin]['raw files']:
                print('Matched {} with {}'.format(gaucho_table['APERO NAME'][i], tbl_nirps['OBJNAME'][imin]))
                print('Nvisit NIRPS --> {}'.format(tbl_nirps[imin]['raw files']))
                print(tbl_nirps[imin])

    tbl_kirkpatrick = Table.read('/Users/eartigau/.hidden_files/Kirkpatrick2024_20pcCensus_Table4_Parameters.csv')

    for i in range(len(gaucho_table)):
        ra = gaucho_table['RA'][i]
        dec = gaucho_table['DEC'][i]
        rad = np.sqrt((tbl_kirkpatrick['RA (deg)'] - ra) ** 2 / np.cos(dec / (180 / np.pi)) + (
                tbl_kirkpatrick['Dec (deg)'] - dec) ** 2)
        imin = np.argmin(rad)
        if not np.isfinite(gaucho_table['Distance pc'][i]):
            print('\nMatched {} with {}'.format(gaucho_table['APERO NAME'][i], tbl_kirkpatrick['DefaultName'][imin]))
            print('\t\tDistance --> {}'.format(1000 / tbl_kirkpatrick[imin]['Parallax (mas)']))

        if not np.isfinite(gaucho_table['H mag'][i]):
            print('\nMatched {} with {}'.format(gaucho_table['APERO NAME'][i], tbl_kirkpatrick['DefaultName'][imin]))
            print('\t\tH mag --> {}'.format(tbl_kirkpatrick[imin]['H (mag)']))

        if not np.isfinite(gaucho_table['MASS'][i]):
            print('\nMatched {} with {}'.format(gaucho_table['APERO NAME'][i], tbl_kirkpatrick['DefaultName'][imin]))
            print('\t\tMass [Msol] --> {}'.format(tbl_kirkpatrick[imin]['AdoptedInitialMass (Msun)']))

    for i in range(len(gaucho_table)):
        ra = gaucho_table['RA'][i]
        dec = gaucho_table['DEC'][i]
        rad = np.sqrt((tbl_planet['ra'] - ra) ** 2 / np.cos(dec / (180 / np.pi)) + (tbl_planet['dec'] - dec) ** 2)
        imin = np.argmin(rad)
        if (np.min(rad) * 3600 < 120) * (gaucho_table['KNOWN_PLANETS'].mask[i]):
            print('There''s a planet around {} --> {}'.format(gaucho_table['APERO NAME'][i], tbl_planet['name'][imin]))

    return gaucho_table


def get_gaucho(make_smart_template=True, no_update=False):
    tbl = gaucho_table(no_update=no_update)
    tbl = tbl[~tbl['REJECTED']]

    search_string = '*e2dsff_tcorr_A.fits'
    search_string2 = '*_pp_s1d_v_tcorr_A.fits'

    for name in np.array(tbl['APERO NAME']):
        local = '/Users/eartigau/lbl_NIRPS_HE/science/{}'.format(name)
        if not os.path.isdir(local):
            os.mkdir(local)

        server = 'nirps-client@rali:/cosmos99/nirps/apero-data/nirps_he_online/objects/{}/{}'.format(name,
                                                                                                     search_string)
        print(name)
        RSYNC_CMD = 'rsync --copy-links -avu -e "ssh -oport=5822" {} {}'.format(server, local)

        print(RSYNC_CMD)
        os.system(RSYNC_CMD)

        local = '/Users/eartigau/smart_template/{}'.format(name)
        if not os.path.isdir(local):
            os.mkdir(local)

        server = 'nirps-client@rali:/cosmos99/nirps/apero-data/nirps_he_online/objects/{}/{}'.format(name,
                                                                                                     search_string2)
        print(name)
        RSYNC_CMD = 'rsync --copy-links -avu -e "ssh -oport=5822" {} {}'.format(server, local)

        print(RSYNC_CMD)
        os.system(RSYNC_CMD)

    if make_smart_template:
        from smart_template import mk_smart

        for name in np.array(tbl['APERO NAME']):
            local = '/Users/eartigau/smart_template/{}/'.format(name)
            if len(glob.glob(local + search_string2)):
                template_name, tbl = mk_smart(name)
                template_name = os.path.abspath(template_name)
                file_name = template_name.split('/')[-1]
                lbl_path = '/Users/eartigau/lbl_NIRPS_HE/templates/'
                outname = lbl_path + file_name
                os.system('rm ' + outname)
                cmd = 'ln -s {} {}'.format(template_name, outname)
                print(cmd)
                os.system(cmd)
                cmd = 'ln -s {} {}'.format(template_name, outname.replace('-VETTED', ''))
                os.system(cmd)

    return


def smart_download():
    print()
    print('copy-paste the path on maestria : \n\n')
    inp = '-'
    path = []
    while len(inp) > 0:
        inp = input()
        path.append(inp)

    for i in range(len(path)):
        if len(path[i]) < 5:
            continue
        cmd = 'scp spirou@maestria:{} .'.format(path[i])
        print(cmd)
        os.system(cmd)


def get_ari():
    # Copy-paste the name of all files for the ARI interface "file_list" document.
    # Press enter twice when done and files get copied to the current directory.
    # Note that unless you have an ssh key, you will be prompted for your password
    # for each file. This is a feature, not a bug. Look at this website for more
    # information on how to set these keys: https://www.ssh.com/ssh/copy-id

    # empty list to store user inputs
    user_inputs = []
    user = ' '  # initialize user input
    i = 1  # initialize counter
    # empty list to store user inputs
    user_inputs = []
    user = ' '  # initialize user input
    i = 1  # initialize counter
    print('Copy-paste the name of all files. Press enter twice when done.\n')
    while user != '':
        user = input()
        user_inputs.append(user)
        i += 1

    # make a numpy array
    user_inputs = np.array(user_inputs)

    # remove entries that do not contain '.fits'
    valid = ['.fits' in user_input for user_input in user_inputs]
    user_inputs = user_inputs[valid]

    # remove entries that do not contain '.fits'
    valid = ['.fits' in user_input for user_input in user_inputs]
    user_inputs = user_inputs[valid]

    for i, user_input in enumerate(user_inputs):
        file_name = user_input.split('/')[-1]
        if os.path.exists(file_name):
            print('File {} already exists. Skipping download.'.format(file_name))
            continue
        print('Downloading [{}/{}] {}'.format(i + 1, len(user_inputs), user_input))
        cmd = 'scp -r spirou@maestria:{} .'.format(user_input)
        # printout for the user
        print(cmd)
        os.system(cmd)


def bandpass_spline(bandpass):
    url = 'http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id={}'

    bandpassnames = dict()
    bandpassnames['TESS'] = 'TESS/TESS.Red'
    bandpassnames['u'] = 'SLOAN/SDSS.u'
    bandpassnames['g'] = 'SLOAN/SDSS.g'
    bandpassnames['r'] = 'SLOAN/SDSS.r'
    bandpassnames['i'] = 'SLOAN/SDSS.i'
    bandpassnames['z'] = 'SLOAN/SDSS.z'
    bandpassnames['Y'] = 'CFHT/Wircam.Y'
    bandpassnames['J'] = 'CFHT/Wircam.J'
    bandpassnames['H'] = 'CFHT/Wircam.H'
    bandpassnames['K'] = 'CFHT/Wircam.Ks'
    bandpassnames['J2M'] = '2MASS/2MASS.J'
    bandpassnames['H2M'] = '2MASS/2MASS.J'
    bandpassnames['K2M'] = '2MASS/2MASS.J'
    bandpassnames['SPIRou'] = 'SPIROU'

    if bandpass not in bandpassnames.keys():
        key_up = np.array([key.upper() for key in bandpassnames.keys()])
        key = np.array([key for key in bandpassnames.keys()])
        if bandpass.upper() in key_up:
            bandpass = key[np.where(key_up == bandpass.upper())[0][0]]
        else:
            raise ValueError('bandpass {} not recognized'.format(bandpass))

    path = hidden_path() + '/filter_' + bandpass + '.txt'
    if not os.path.exists(path):
        r = requests.get(url.format(bandpassnames[bandpass]))
        print('Downloading {}'.format(bandpass))
        with open(path, 'wb') as f:
            f.write(r.content)

    tbl = Table.read(path, format='ascii')
    # convert to float arras, import as some values may be integers
    tbl['col1'] = np.array(tbl['col1'], dtype=float)
    tbl['col2'] = np.array(tbl['col2'], dtype=float)
    # col1 -> wavelength in nm
    tbl['col1'] /= 10.0  # form Angstrom to nm
    # col2 -> transmission
    tbl['col2'][0] = 0  # for spline of out-of-band, set first value to 0
    tbl['col2'][-1] = 0  # for spline of out-of-band, set last value to 0

    tbl['col1'] = np.array(tbl['col1'], dtype=float)
    tbl['col2'] = np.array(tbl['col2'], dtype=float)
    return ius(tbl['col1'], tbl['col2'], k=1, ext=1)


def get_limb_spectrum(mu, teff, mask=None, logg=4.5, intensity=False, wavelength=None):
    """

    :param mu: 1D or 2D map of mu, should range between 0 and 1, can include NaNs
    :param teff: Stellar model temperature
    :param mask: True for valid, False for invalid
    :param logg: Stellar model logg
    ---->
    :return: wavelength and flux
    """

    limb = get_specific_intensity(teff, logg=logg)

    if wavelength is not None:
        keep = (limb['wavelength_nm'] > wavelength[0]) & (limb['wavelength_nm'] < wavelength[1])
        limb['wavelength_nm'] = limb['wavelength_nm'][keep]
        limb['intensity'] = limb['intensity'][:, keep]

    spectrum_out = np.zeros_like(limb['wavelength_nm'])

    if mask is None:
        mask = np.isfinite(mu)
    else:
        mask &= np.isfinite(mu)

    mu2 = np.append(np.append(0, limb['mu']), 99)
    mu_masked = mu[mask]
    int2 = np.zeros_like(mu_masked)

    if intensity:
        intensity_map = np.zeros_like(mu)

    for i in range(len(limb['mu'])):
        low = mu2[i + 1] - (mu2[i + 1] - mu2[i]) / 2
        high = mu2[i + 1] + (mu2[i + 2] - mu2[i + 1]) / 2
        valid = (mu_masked >= low) & (mu_masked < high)

        if not intensity:
            # express in photos
            spectrum_out += limb['intensity'][i] * np.sum(valid) * limb['wavelength_nm']
        else:

            int2[valid] = np.sum(limb['wavelength_nm'] * limb['intensity'][i])

    if not intensity:
        return limb['wavelength_nm'], spectrum_out
    else:
        intensity_map[mask] = int2
        return intensity_map


def get_specific_intensity(teff, logg=4.5):
    loggtxt = str(logg).ljust(4, '0')
    tefftxt = str(int(teff)).rjust(5, '0')
    file = 'lte{}-{}-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits'.format(tefftxt, loggtxt)
    url = ' wget ftp://phoenix.astro.physik.uni-goettingen.de/SpecIntFITS/PHOENIX-ACES-AGSS-COND-SPECINT-2011/Z-0.0' \
          '/{}'.format(file)

    outname = hidden_path() + '/stellar_models/'
    if not os.path.exists(outname):
        os.makedirs(outname)
    outname2 = outname + file
    outname3 = outname + file.replace('.fits', '.pick')
    if not os.path.exists(outname2):
        os.system(url + ' -O {}'.format(outname2))

    if not os.path.isfile(outname3):
        tmp = dict()
        intensity = fits.getdata(outname2, 'PRIMARY')
        h = fits.getheader(outname2, 'PRIMARY')
        mu = fits.getdata(outname2, 'MU')
        wave = np.arange(intensity.shape[1]) * h['CDELT1'] + h['CRVAL1']
        tmp['wavelength_nm'] = wave / 10.
        tmp['intensity'] = intensity
        tmp['mu'] = mu
        save_pickle(outname3, tmp)
        print('Wrote {}'.format(outname3))
    else:
        tmp = read_pickle(outname3)

    return tmp


def model_flux(wavelength, teff, frac_window=0.1, logg=4.5):
    """
    :param wavelength: wavelength with units
    :param teff: effective temperature of star with units
    : optional parameters:
    :    frac_window: fractional window width, default = 10%
    :   logg: surface gravity of star, default = 4.5
    :
    :return: mean flux in domain of interest splined over Goettingen models
    """

    # get the limb darkening coefficients for the star at the wavelength of interest

    if type(teff) is uu.quantity.Quantity:
        teff = teff.to(uu.K).value
    if type(wavelength) is uu.quantity.Quantity:
        wavelength = wavelength.to(uu.nm).value

    wave0 = (1 - frac_window / 2) * wavelength
    wave1 = (1 + frac_window / 2) * wavelength

    if len(teff) == 1:
        teff0 = np.round(teff, -2) - 400
        teff1 = np.round(teff, -2) + 500
    else:
        teff0 = np.round(np.min(teff), -2) - 400
        teff1 = np.round(np.max(teff), -2) + 500

    teffs = range(int(teff0), int(teff1), 100)

    teff_list = []
    flux = []

    for i in tqdm(range(len(teffs))):
        try:
            sp = get_goettingen_model(teff=teffs[i], logg=logg, wave0=wave0,
                                      wave1=wave1, step=200000)
            flux.append(np.mean(sp['flux']))
            teff_list.append(teffs[i])
        except:
            pass
    flux = np.array(flux)
    teff_list = np.array(teff_list)
    return ius(teff_list, flux, k=2, ext=1)(teff)


def snail(iter, desc='', leave=False):
    txt = [' ']
    for i in range(1000):
        txt.append('🐌')
    txt.append('_')
    # ,'🐌','🐌','🐌','🐌','🐌','🐌','🐌','🐌','_']
    return tqdm(iter, leave=leave, desc=desc,
                colour='green', ascii=txt)


def squid(iter, desc='', leave=False):
    txt = [' ']
    for i in range(1000):
        txt.append('🦑')
    txt.append('_')
    return tqdm(iter, leave=leave, desc=desc,
                colour='blue', ascii=txt)


def save_pickle(filename, variable):
    with open(filename, 'wb') as handle:
        pickle.dump(variable, handle)


def read_pickle(filename):
    with open(filename, 'rb') as handle:
        b = pickle.load(handle)
    return b


def which_band(wave):
    """SLOAN/SDSS.u	3556.52	3572.18	3608.04	3055.11	4030.64	540.97	1582.54	3.75e-9	SLOAN	 	SDSS u full transmission
    Filter ID	λref	λmean	λeff	λmin	λmax	Weff	ZPν	ZPλ
    SLOAN/SDSS.u	3556.52	3572.18	3608.04	3055.11	4030.64	540.97	1582.54	3.75e-9	SLOAN
    SLOAN/SDSS.g	4702.50	4750.82	4671.78	3797.64	5553.04	1064.68	4023.57	5.45e-9	SLOAN	 	SDSS g full transmission
    SLOAN/SDSS.r	6175.58	6204.29	6141.12	5418.23	6994.42	1055.51	3177.38	2.5e-9	SLOAN	 	SDSS r full transmission
    SLOAN/SDSS.i	7489.98	7519.27	7457.89	6692.41	8400.32	1102.57	2593.40	1.39e-9	SLOAN	 	SDSS i full transmission
    SLOAN/SDSS.z	8946.71	8992.26	8922.78	7964.70	10873.33	1164.01	2238.99	8.39e-10	SLOAN	 	SDSS z full transmission
    Filter ID	λref	λmean	λeff	λmin	λmax	Weff	ZPν	ZPλ	Obs. Facility	Instrument	Description
    2MASS/2MASS.J	12350.00	12350.00	12350.00	10806.47	14067.97	1624.32	1594.00	3.13e-10	2MASS	 	2MASS J
    2MASS/2MASS.H	16620.00	16620.00	16620.00	14787.38	18231.02	2509.40	1024.00	1.13e-10	2MASS	 	2MASS H
    2MASS/2MASS.Ks	21590.00	21590.00	21590.00	19543.69	23552.40	2618.87	666.80	4.28e-11	2MASS
    """

    band = ['u', 'g', 'r', 'i', 'z', 'J', 'H', 'K']
    # wave_min = [305.511, 379.764, 541.823, 669.241, 796.470, 1080.647, 1478.738, 1954.369]
    # wave_max = [403.064, 555.304, 699.442, 840.032, 1087.333, 1406.797, 1823.102, 2355.240]
    wave_mean = [360.804, 467.178, 614.112, 745.789, 892.278, 1235.000, 1662.000, 2159.000] * uu.nm

    # u       v       b       y       U       B       V       R       I       J       H       K

    return band[np.argmin(np.abs(wave - wave_mean))]


def limb_darkening(teff=5000 * uu.K, logg=4.5, wave=1600 * uu.nm):
    """
    :param teff: effective temperature in K
    :param band: 'H' or 'K'
    :return: limb darkening coefficient
    """
    url = 'https://cdsarc.cds.unistra.fr/viz-bin/nph-Cat/txt?J/A+A/363/1081/phoenix.dat.gz'

    path = os.path.join(hidden_path(), 'phoenix.dat.gz').replace('.dat.gz', '.fits')
    if not os.path.exists(path):
        r = requests.get(url, allow_redirects=True)
        open(path, 'wb').write(r.content)

        tt = Table.read(path, format='ascii', delimiter='&', guess=False)
        tt = tt[['|' in tt[i][0] for i in range(len(tt))]][1:]
        tt = tt[['-------' not in tt[i][0] for i in range(len(tt))]]
        tt = [tt[i][0] for i in range(len(tt))]
        keys = [key.strip() for key in tt[0].split('|')]
        tbl = Table()
        for key in keys:
            if key == 'eff':
                key_type = 'U4'
            else:
                key_type = float
            tbl[key] = np.zeros(len(tt[1:]), dtype=key_type)
        for i in range(len(tt[1:])):
            for j, key in enumerate(keys):
                tbl[key][i] = tt[i + 1].split('|')[j].strip()

        tbl['Teff'] = tbl['Teff'] * uu.K
        tbl.write(path, format='fits', overwrite=True)

    tbl = Table.read(path, format='fits')

    # select right temperature
    tbl = tbl[tbl['Teff'][np.argmin(np.abs(tbl['Teff'] - teff))] == tbl['Teff']]
    tbl = tbl[tbl['logg'][np.argmin(np.abs(tbl['logg'] - logg))] == tbl['logg']]

    band = which_band(wave)
    if band == 'g':
        band = 'V'
    if band == 'i':
        band = 'I'

    if band not in tbl.keys() and band.upper() in tbl.keys():
        band = band.upper()

    return np.array(tbl[band])


def hidden_path():
    path = os.path.expanduser('~')
    hidden = os.path.join(path, '.etienne_tools')
    if not os.path.exists(hidden):
        os.makedirs(hidden)

    return hidden


def planck(wav, T):
    # define a black body function
    a = 2.0 * constants.h * constants.c ** 2
    b = constants.h * constants.c / (wav * constants.k_B * T)
    intensity = a / ((wav ** 5) * (np.exp(b) - 1.0))
    return intensity


def map2sphere(map_lat_lon, xpix, ypix):
    """

    :param map_lat_lon: A map that is a function of latitude and longitude, does not have to be normalized
    360x180, but may be any 2D size. Will be interpolated to 360x180

    :param xpix: xpixels of the output map, inclide the -1 to +1 domain, will be set to NaN outside x**2+y**2>1
    :param ypix: same as above
    :return: returns the projected map, a 2D with the same size as xpix and ypix
    """

    sphereimage = np.zeros_like(xpix) + np.nan

    resampling_theta = map_lat_lon.shape[0] / 180
    resampling_phi = map_lat_lon.shape[1] / 360

    zpix2 = (1 - xpix ** 2 - ypix ** 2)
    zpix2[zpix2 < 0] = np.nan
    zpix = np.sqrt(zpix2)
    rho = np.sqrt(xpix ** 2 + ypix ** 2 + zpix ** 2)
    phi = np.arctan2(ypix, zpix) / np.pi * 180 + 90
    theta = np.arccos(xpix / rho) / np.pi * 180
    g = np.isfinite(phi) & np.isfinite(theta)

    theta = theta
    phi = phi
    theta[theta > 179] = 179

    theta_int = (theta * resampling_theta).astype(int)
    phi_int = (phi * resampling_phi).astype(int)
    sphereimage[g] = map_lat_lon[theta_int[g], phi_int[g]]

    return sphereimage


def get_goettingen_model(teff=3000, logg=5.0, wave0=700, wave1=2800, step=200000):
    # who am I?
    hidden = hidden_path()
    path_to_models = os.path.join(hidden, 'stellar_models')
    if not os.path.exists(path_to_models):
        os.makedirs(path_to_models)

    binned_model_name = 'teff{}-logg{}-{}-{}-step{}.npy'.format(teff, logg, wave0, wave1, step)
    binned_model_name = os.path.join(path_to_models, binned_model_name)
    if os.path.isfile(binned_model_name):
        tbl = Table(np.load(binned_model_name))
        print('Reading {}'.format(binned_model_name))
        # tbl = Table.read(binned_model_name,format='ascii')
        return tbl

    grid_url = 'ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
    wavegridpath = '{}/wavegrid.fits'.format(path_to_models)
    if not os.path.isfile(wavegridpath):
        os.system(' wget {} -O {}'.format(grid_url, wavegridpath))
    wave = fits.getdata(wavegridpath) / 10

    params = str(int(teff)).zfill(5), logg
    outname = 'lte{}-{:.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(*params)
    if not os.path.isfile('{}/{}'.format(path_to_models, outname)):
        get_url = 'wget ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/' + outname
        cmd = 'wget {} -O {}/{}'.format(get_url, path_to_models, outname)
        print(cmd)
        os.system(cmd)

    sp = fits.getdata('{}/{}'.format(path_to_models, outname))
    # trim the spectrum to the desired wavelength range
    cond1 = (wave > (wave0 * (1 - 20 * step / constants.c.value)))
    cond2 = (wave < (wave1 * (1 + 20 * step / constants.c.value)))
    keep = cond1 & cond2
    wave = wave[keep]
    sp = sp[keep]

    # step in velocity of the sparsed grid step
    step0 = constants.c.value / (np.min(wave / np.gradient(wave)))
    wave2 = get_magic_grid(wave0=wave0, wave1=wave1, dv_grid=step)

    if step0 < step / 5.0:
        # we need at least 5 points per resolution element

        sp2 = np.zeros_like(wave2)
        for i in tqdm(range(len(wave2)), desc='Interpolating spectrum',
                      leave=False):
            mask = (np.abs(wave / wave2[i] - 1) * constants.c.value < step / 2)
            sp2[i] = np.mean(sp[mask])

        tbl = Table([wave2, sp2], names=['wave', 'flux'])
        np.save(binned_model_name, tbl)

    else:
        sp2 = ius(wave, sp, k=2)(wave2)

        tbl = Table([wave2, sp2], names=['wave', 'flux'])
        np.save(binned_model_name, tbl)

    return tbl


def nanpercentile(v, p, axis=None):
    if axis == None:
        return jit_nanpercentile(v, p)
    else:
        return np.nanpercentile(v, p, axis=axis)


@jit(nopython=True)
def jit_nanpercentile(v, p):
    return np.nanpercentile(v, p)


@jit(nopython=True)
def nanstd(v):
    return np.nanstd(v)


@jit(nopython=True)
def nansum(v):
    return np.nansum(v)


@jit(nopython=True)
def sum(v):
    return np.sum(v)


@jit(nopython=True)
def mean(v):
    return np.mean(v)


@jit(nopython=True)
def std(v):
    return np.std(v)


def nanmean(v):
    with warnings.catch_warnings(record=True) as _:
        mean = jitnanmean(v)
    return mean


@jit(nopython=True)
def jitnanmean(v):
    mean = np.nanmean(v)
    return mean


@jit(nopython=True)
def nanmedian(v):
    med = np.nanmedian(v)
    # g = np.isfinite(v)
    return med


@jit(nopython=True)
def median(v):
    return np.median(v)


@jit(nopython=True)
def exp(v):
    return np.exp(v)


# Define Functions Using Numba
# Idea here is to solve ax = b, using least squares, where a represents our coefficients e.g. x**2, x, constants
@jit(nopython=True)
def _coeff_mat(x, deg):
    mat_ = np.zeros(shape=(x.shape[0], deg + 1))
    const = np.ones_like(x)
    mat_[:, 0] = const
    mat_[:, 1] = x
    if deg > 1:
        for n in range(2, deg + 1):
            mat_[:, n] = x ** n
    return mat_


@jit(nopython=True)
def _fit_x(a, b):
    # linalg solves ax = b
    det_ = np.linalg.lstsq(a, b)[0]
    return det_


@jit(nopython=True)
def fit_poly(x, y, deg):
    a = _coeff_mat(x, deg)
    p = _fit_x(a, y)
    # Reverse order so p[0] is coefficient of highest order
    return p[::-1]


@jit(nopython=True)
def eval_polynomial(P, x):
    '''
    Compute polynomial P(x) where P is a vector of coefficients, highest
    order coefficient at P[0].  Uses Horner's Method.
    '''
    result = 0
    for coeff in P:
        result = x * result + coeff
    return result


def art(word, color1='MAGENTA', color2='red'):
    letter = \
        ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
         'W', 'X', 'Y', 'Z']
    length = \
        [3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 4, 3, 4, 3, 3]

    low1 = "┌─┐┌┐ ┌─┐┌┬┐┌─┐┌─┐┌─┐┬ ┬┬  ┬┬┌─┬  ┌┬┐┌┐┌┌─┐┌─┐┌─┐ ┬─┐┌─┐┌┬┐┬ ┬┬  ┬┬ ┬─┐ ┬┬ ┬┌─┐"
    low2 = "├─┤├┴┐│   ││├┤ ├┤ │ ┬├─┤│  │├┴┐│  │││││││ │├─┘│─┼┐├┬┘└─┐ │ │ │└┐┌┘│││┌┴┬┘└┬┘┌─┘"
    low3 = "┴ ┴└─┘└─┘─┴┘└─┘└  └─┘┴ ┴┴└─┘┴ ┴┴─┘┴ ┴┘└┘└─┘┴  └─┘└┴└─└─┘ ┴ └─┘ └┘ └┴┘┴ └─ ┴ └─┘"
    up1 = "╔═╗╔╗ ╔═╗╔╦╗╔═╗╔═╗╔═╗╦ ╦╦  ╦╦╔═╦  ╔╦╗╔╗╔╔═╗╔═╗╔═╗ ╦═╗╔═╗╔╦╗╦ ╦╦  ╦╦ ╦═╗ ╦╦ ╦╔═╗"
    up2 = "╠═╣╠╩╗║   ║║║╣ ╠╣ ║ ╦╠═╣║  ║╠╩╗║  ║║║║║║║ ║╠═╝║═╬╗╠╦╝╚═╗ ║ ║ ║╚╗╔╝║║║╔╩╦╝╚╦╝╔═╝"
    up3 = "╩ ╩╚═╝╚═╝═╩╝╚═╝╚  ╚═╝╩ ╩╩╚═╝╩ ╩╩═╝╩ ╩╝╚╝╚═╝╩  ╚═╝╚╩╚═╚═╝ ╩ ╚═╝ ╚╝ ╚╩╝╩ ╚═ ╩ ╚═╝"

    letter = letter + ['-', ' ', '.', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '[', ']', '?', '!']
    length = length + [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
    low1 = low1 + "         ┌─┐ ┐ ┌─┐┌─┐┌ ┐┌─┐┌─┐┌─┐┌─┐┌─┐ ┌  ┐ ┌─┐ ┐ "
    low2 = low2 + "───      │ │ │ ┌─┘ ─┤└─┤└─┐├─┐  │├─┤└─┤ │  │  ┌┘ │ "
    low3 = low3 + "      ·  └─┘─┴─└─┘└─┘  ┘└─┘└─┘  ┴└─┘└─┘ └  ┘  o  o "

    low_1 = ""
    low_2 = ""
    low_3 = ""

    letter = np.array([l.lower() for l in letter])

    l1 = np.array([(np.cumsum(length))[l.lower() == letter][0] for l in word])
    l2 = np.array([np.array(length)[l.lower() == letter][0] for l in word])
    l0 = l1 - l2

    for i in range(len(l1)):
        if word[i] == word[i].lower():
            low_1 += low1[l0[i]:l1[i]]
            low_2 += low2[l0[i]:l1[i]]
            low_3 += low3[l0[i]:l1[i]]
        else:
            low_1 += up1[l0[i]:l1[i]]
            low_2 += up2[l0[i]:l1[i]]
            low_3 += up3[l0[i]:l1[i]]

    low_0 = color('╔' + '═' * (len(low_1) + 2) + '╗', color1)
    low_4 = color('╚' + '═' * (len(low_1) + 2) + '╝', color1)

    low_1 = color('║ ', color1) + color(low_1, color2) + color(' ║', color1)
    low_2 = color('║ ', color1) + color(low_2, color2) + color(' ║', color1)
    low_3 = color('║ ', color1) + color(low_3, color2) + color(' ║', color1)

    try:
        w = os.get_terminal_size().columns
    except OSError:
        w = 80
    dw = (w - len(low_1) // 2) // 2
    low_0 = ' ' * dw + low_0
    low_1 = ' ' * dw + low_1
    low_2 = ' ' * dw + low_2
    low_3 = ' ' * dw + low_3
    low_4 = ' ' * dw + low_4

    return '\n' + low_0 + '\n' + low_1 + '\n' + low_2 + '\n' + low_3 + '\n' + low_4 + '\n'


def sigma(im_tmp):
    im = np.array(im_tmp).astype('float64')
    p1 = (1 - 0.682689492137086) / 2
    p2 = 1 - p1
    return (nanpercentile(im, p2 * 100) - nanpercentile(im, p1 * 100)) / 2


def running_sigma(v, w):
    # provide a vector and return a running robust dispersion
    # the width (w) is the box size to measure the dispersion
    # pick a w that is wide enough to be representative of the
    # noise but small enough to avoid variations
    ll = len(v)

    sig = np.zeros(ll)
    for i in range(ll):
        i0 = i - w // 2
        i1 = i + w // 2
        if i0 < 0:
            i0 = 0
        if i1 > ll:
            i1 = ll
        sig[i] = sigma(v[i0:i1])

    return sig


def color(message, color):
    COLOURS = dict()
    COLOURS['BLACK'] = '\033[90;1m'
    COLOURS['RED'] = '\033[1;91;1m'
    COLOURS['GREEN'] = '\033[92;1m'
    COLOURS['YELLOW'] = '\033[1;93;1m'
    COLOURS['BLUE'] = '\033[94;1m'
    COLOURS['MAGENTA'] = '\033[1;95;1m'
    COLOURS['CYAN'] = '\033[1;96;1m'
    COLOURS['WHITE'] = '\033[97;1m'
    COLOURS['ENDC'] = '\033[0;0m'

    return COLOURS[color.upper()] + message + COLOURS['ENDC']


def printc(message, msg_type='', print_time=True):
    """
    Print a message with a color
    :param message:
    :param msg_type:
        -> info = green
        -> bad1 = yellow
        -> bad2 = red
        -> bad3 = magenta
        -> number = blue
        -> (other) = white
    :param print_time:
    :return: nothing
    """

    msg_color = "black"

    if msg_type == 'info':
        msg_color = 'green'

    if msg_type == 'bad1':
        msg_color = 'cyan'

    if msg_type == 'bad2':
        msg_color = 'red'

    if msg_type == 'bad3':
        msg_color = 'magenta'

    if msg_type == 'number':
        msg_color = 'blue'

    if print_time:
        time = datetime.now().strftime('%H:%M:%S.%f')[:-4] + '│'
    else:
        time = ''

    if len(message) == 1:
        # get terminal width
        try:
            w = os.get_terminal_size().columns - len(time)
        except:
            w = 80 - len(time)
        message = message[0] * w

    print(color(time + message, msg_color))


def get_magic_grid(wave0=1500, wave1=1800, dv_grid=0.5):
    # default for the function is 500 m/s
    # the arithmetic is a but confusing here, you first find how many
    # elements you have on your grid, then pass it to an exponential
    # the first element is exactely wave0, the last element is NOT
    # exactly wave1, but is very close and is set to get your exact
    # step in velocity
    len_magic = int(np.ceil(np.log(wave1 / wave0) * np.array(constants.c) / dv_grid))
    magic_grid = np.exp(np.arange(len_magic) / len_magic * np.log(wave1 / wave0)) * wave0
    return magic_grid


def doppler(wave, v):
    # velocity expressed in m/s
    # relativistic calculation

    v = np.array(v)
    wave = np.array(wave)
    return wave * np.sqrt((1 - v / constants.c.value) / (1 + v / constants.c.value))


def weighted_median(values, weights):
    keep = np.isfinite(values) * np.isfinite(weights)
    values1 = np.array(values[keep], dtype=float)
    weights1 = np.array(weights[keep], dtype=float)
    weights1 /= np.nansum(weights1)

    ord = np.argsort(values1)
    values1 = values1[ord]
    weights1 = weights1[ord]

    cumsum = np.cumsum(weights1)
    imed = np.min(np.where(cumsum > 0.5))

    return values1[imed]


def wave2wave(e2ds_data_input, wave1, wave2):
    # transform e2ds data from one wavelength grid to another.

    e2ds_data = np.array(e2ds_data_input)

    for iord in range(49):
        keep = np.isfinite(e2ds_data[iord])
        spl = ius(wave1[iord][keep], e2ds_data[iord][keep], k=3, ext=1)
        e2ds_data[iord][keep] = spl(wave2[iord][keep])

    return e2ds_data


def fit_gauss(x, y, p0):
    fit, pcov = curve_fit(gauss, x, y, p0=p0)
    return fit


def get_rough_ccf_rv(wave, sp, wave_mask, weight_line, doplot=False):
    if len(wave.shape) == 2:  # we have the e2ds file, we reshape it
        mask = np.ones_like(wave, dtype=bool)

        for iord in range(1, 49):
            mask[iord] *= (wave[iord - 1, ::-1] < wave[iord])

        for iord in range(0, 48):
            mask[iord] *= (wave[iord] < wave[iord + 1, ::-1])

        mask *= np.isfinite(sp)

        sp2 = sp[mask]
        wave2 = wave[mask]
    else:  # we have the s1d, we just need to get rid of NaNs
        sp2 = sp[np.isfinite(sp)]
        wave2 = wave[np.isfinite(sp)]

    spline_sp = ius(wave2, sp2, k=1, ext=1)

    dvs = np.arange(-1.5e5, 1.5e5, 500)
    ccf = np.zeros_like(dvs)
    print('computing CCF')
    for i in range(len(dvs)):
        ccf[i] = np.nansum(weight_line * spline_sp(doppler(wave_mask, dvs[i])))

    imax = np.argmax(ccf)

    guess = [dvs[imax], 2000, ccf[imax] - nanmedian(ccf), nanmedian(ccf), 0]
    fit, pcov = curve_fit(gauss, dvs, ccf, p0=guess)

    # fit = np.polyfit(dvs[imax - 1:imax + 2], ccf[imax - 1:imax + 2], 2)
    systemic_vel = fit[0]  # -.5 * fit[1] / fit[0]
    print('CCF Velocity : {0:.2f} m/s'.format(-systemic_vel))

    if doplot:
        plt.plot(-dvs / 1000, ccf)
        plt.plot(-dvs / 1000, gauss(dvs, *fit))
        plt.xlabel('RV [km/s]')
        plt.ylabel('Normalized CCF')
        plt.show()

    return systemic_vel


def robust_polyfit(x, y, degree, nsigcut, accept_width=None):
    # if we only want to reject beyond a certain width, we set the accept_width
    # to a given distance to fit.
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    degree = np.array(degree, dtype=int)

    if accept_width is not None:
        nsigcut = 1.0

    keep = np.isfinite(y)
    # set the nsigmax to infinite
    nsigmax = np.inf
    # set the fit as unset at first
    fit = None
    # while sigma is greater than sigma cut keep fitting
    while nsigmax > nsigcut:
        # calculate the polynomial fit (of the non-NaNs)
        fit = fit_poly(x[keep], y[keep], degree)
        # calculate the residuals of the polynomial fit
        res = y - np.polyval(fit, x)
        # work out the new sigma values
        if accept_width is None:
            sig = nanmedian(np.abs(res))
        else:
            nsigcut = 1
            sig = accept_width

        if sig == 0:
            nsig = np.zeros_like(res)
            nsig[res != 0] = np.inf
        else:
            nsig = np.abs(res) / sig
        # work out the maximum sigma
        nsigmax = np.max(nsig[keep])
        # re-work out the keep criteria
        keep = nsig < nsigcut
    # return the fit and the mask of good values
    return fit, keep


def sigma(tmp):
    if type(tmp[0]) != np.float64:
        tmp = np.array(tmp, dtype='float64')
    # return a robust estimate of 1 sigma
    sig1 = 0.682689492137086
    p1 = (1 - (1 - sig1) / 2) * 100
    return (nanpercentile(tmp, p1) - nanpercentile(tmp, 100 - p1)) / 2.0


def gauss(x, cen, ew, amp, zp, slope):
    return np.exp(-0.5 * (x - cen) ** 2 / ew ** 2) * amp + zp + (x - cen) * slope


def lowpassfilter(input_vect, width=101):
    # Computes a low-pass filter of an input vector. This is done while properly handling
    # NaN values, but at the same time being reasonably fast.
    # Algorithm:
    #
    # provide an input vector of an arbitrary length and compute a running NaN median over a
    # box of a given length (width value). The running median is NOT computed at every pixel
    # but at steps of 1/4th of the width value. This provides a vector of points where
    # the nan-median has been computed (ymed) and mean position along the input vector (xmed)
    # of valid (non-NaN) pixels. This xmed/ymed combination is then used in a spline to
    # recover a vector for all pixel positions within the input vector.
    #
    # When there are no valid pixel in a 'width' domain, the value is skipped in the creation
    # of xmed and ymed, and the domain is splined over.

    # indices along input vector
    index = np.arange(len(input_vect))

    # placeholders for x and y position along vector
    xmed = []
    ymed = []

    # loop through the lenght of the input vector
    for i in np.arange(-width // 2, len(input_vect) + width // 2, width // 4):

        # if we are at the start or end of vector, we go 'off the edge' and
        # define a box that goes beyond it. It will lead to an effectively
        # smaller 'width' value, but will provide a consistent result at edges.
        low_bound = i
        high_bound = i + int(width)

        if low_bound < 0:
            low_bound = 0
        if high_bound > (len(input_vect) - 1):
            high_bound = (len(input_vect) - 1)

        pixval = index[low_bound:high_bound]

        if len(pixval) < 3:
            continue

        # if no finite value, skip
        if np.max(np.isfinite(input_vect[pixval])) == 0:
            continue

        # mean position along vector and NaN median value of
        # points at those positions
        xmed.append(nanmean(pixval))
        ymed.append(nanmedian(input_vect[pixval]))

    xmed = np.array(xmed, dtype=float)
    ymed = np.array(ymed, dtype=float)

    # we need at least 3 valid points to return a
    # low-passed vector.
    if len(xmed) < 3:
        return np.zeros_like(input_vect) + np.nan

    if len(xmed) != len(np.unique(xmed)):
        xmed2 = np.unique(xmed)
        ymed2 = np.zeros_like(xmed2)
        for i in range(len(xmed2)):
            ymed2[i] = np.mean(ymed[xmed == xmed2[i]])
        xmed = xmed2
        ymed = ymed2

    # splining the vector
    spline = ius(xmed, ymed, k=2, ext=3)
    lowpass = spline(np.arange(len(input_vect)))

    return lowpass


def file_search(search_string):
    out = np.array(glob.glob(search_string))
    out = out[np.argsort(out)]
    return out


def fits2wave(file_or_header):
    info = """
        Provide a fits header or a fits file
        and get the corresponding wavelength
        grid from the header.

        Usage :
          wave = fits2wave(hdr)
                  or
          wave = fits2wave('my_e2ds.fits')

        Output has the same size as the input
        grid. This is derived from NAXIS 
        values in the header
    """

    # check that we have either a fits file or an astropy header
    if type(file_or_header) == str:
        hdr = fits.getheader(file_or_header, ext=0)
        if 'NAXIS1' not in hdr:
            hdr['NAXIS1'] = 4088
            hdr['NAXIS2'] = 49
    elif str(type(file_or_header)) == "<class 'astropy.io.fits.header.Header'>":
        hdr = file_or_header

        if 'NAXIS1' not in hdr:
            hdr['NAXIS1'] = 4088
            hdr['NAXIS2'] = 49
    else:
        print()
        print('~~~~ wrong type of input ~~~~')
        print()

        print(info)
        return []

    if 'WAVEPOLY' in hdr:
        if hdr['WAVEPOLY'] == 'Chebyshev':
            # get the keys with the wavelength polynomials
            wave_hdr = hdr['WAVE0*']
            # concatenate into a numpy array
            wave_poly = np.array([wave_hdr[i] for i in range(len(wave_hdr))])

            # get the number of orders
            nord = hdr['WAVEORDN']

            # get the per-order wavelength solution
            wave_poly = wave_poly.reshape(nord, len(wave_poly) // nord)

            # get the length of each order (normally that's 4088 pix)
            npix = hdr['NAXIS1']

            # project polynomial coefficiels
            wavesol = [np.polynomial.chebyshev.chebval(np.arange(npix) / npix, wave_poly[i]) for i in range(nord)]
            return np.array(wavesol)
        else:
            raise ValueError('WAVEPOLY not recognized')

    else:
        # get the keys with the wavelength polynomials
        wave_hdr = hdr['WAVE0*']
        # concatenate into a numpy array
        wave_poly = np.array([wave_hdr[i] for i in range(len(wave_hdr))])

        # get the number of orders
        nord = hdr['WAVEORDN']

        # get the per-order wavelength solution
        wave_poly = wave_poly.reshape(nord, len(wave_poly) // nord)

        # get the length of each order (normally that's 4088 pix)
        npix = hdr['NAXIS1']

        # project polynomial coefficiels
        wavesol = [np.polyval(wave_poly[i][::-1], np.arange(npix)) for i in range(nord)]

    # return wave grid
    return np.array(wavesol)


def smart_time(time_in_s):
    if time_in_s > 3600:
        h = str(int(np.floor(time_in_s / 3600))) + 'h'
        time_in_s -= (np.floor(time_in_s / 3600) * 3600)
        flag_h = True
    else:
        h = ''
        flag_h = False

    if time_in_s > 60:
        minimum = str(int(np.floor(time_in_s / 60))).zfill(2) + 'm'
        time_in_s -= (np.floor(time_in_s / 60) * 60)
    else:
        minimum = ''

    if flag_h:
        return h + minimum
    else:
        sec = str(int(np.round(time_in_s))).zfill(2) + 's'
        return minimum + sec


def td_convert(instance):
    if isinstance(instance, Table):
        out = dict()
        for col in instance.keys():
            out[col] = np.array(instance[col])
        return out
    if isinstance(instance, dict):
        out = Table()
        for col in instance.keys():
            out[col] = np.array(instance[col])
        return out


def running_rms(sp1):
    sp1b = np.zeros(4096) + np.nan
    sp1b[4:-4] = sp1
    with warnings.catch_warnings(record=True) as _:
        b1 = np.nanpercentile(np.reshape(sp1b, [16, 256]), [16, 84], axis=1)
    rms = (b1[1] - b1[0]) / 2
    index = np.arange(16) * 256 + 128
    keep = np.isfinite(rms)
    index = index[keep]
    rms = rms[keep]

    return ius(index, rms, k=2, ext=3)(np.arange(len(sp1)))


def sed_ratio(sp1, sp2):
    sp1b = np.zeros(4096) + np.nan
    sp2b = np.zeros(4096) + np.nan
    sp1b[4:-4] = sp1
    sp2b[4:-4] = sp2

    invalid = (np.isfinite(sp1b) * np.isfinite(sp2b)) == False
    sp1b[invalid] = np.nan
    sp2b[invalid] = np.nan

    index = np.arange(128) * 32 + 16
    b1 = np.nansum(np.reshape(sp1b, [128, 32]), axis=1)
    b2 = np.nansum(np.reshape(sp2b, [128, 32]), axis=1)

    invalid = ((b1 != 0) * (b2 != 0)) == False
    b1[invalid] = np.nan
    b2[invalid] = np.nan

    ratio = b1 / b2

    # fit,_ = robust_polyfit(index,ratio,3,3)
    # return np.polyval(fit,np.arange(len(sp1)))

    # plt.plot(index,ratio)
    # plt.plot(index,np.polyval(fit,index))
    # plt.show()

    ratio2 = np.zeros_like(ratio) + np.nan
    for i in range(len(ratio)):
        if np.isfinite(ratio[i]):
            i1 = i - 3
            i2 = i + 4
            ratio2[i] = nanmedian(ratio[i1:i2])

    keep = np.isfinite(ratio2)
    index = index[keep]
    ratio2 = ratio2[keep]

    if len(ratio2) < 4:
        return np.zeros_like(sp1) + nanmedian(ratio2)
    else:
        return ius(index, ratio2, k=2, ext=3)(np.arange(len(sp1)))


@jit(nopython=True)
def odd_ratio_mean(value, err, odd_ratio=1e-4, nmax=10):
    #
    # Provide values and corresponding errors and compute a
    # weighted mean
    #
    #
    # odd_bad -> probability that the point is bad
    #
    # nmax -> number of iterations
    keep = np.isfinite(value) * np.isfinite(err)

    if np.sum(keep) == 0:
        return np.nan, np.nan

    value = value[keep]
    err = err[keep]

    guess = np.nanmedian(value)

    nite = 0
    while (nite < nmax):
        nsig = (value - guess) / err
        gg = np.exp(-0.5 * nsig ** 2)
        odd_bad = odd_ratio / (gg + odd_ratio)
        odd_good = 1 - odd_bad

        w = odd_good / err ** 2

        guess = nansum(value * w) / np.nansum(w)
        nite += 1

    bulk_error = np.sqrt(1 / np.nansum(odd_good / err ** 2))

    return guess, bulk_error


def lin_mini(vector, sample):
    # wrapper function that sets everything for the @jit later
    # In particular, we avoid the np.zeros that are not handled
    # by numba

    # size of input vectors and sample to be adjusted
    sz_sample = sample.shape  # 1d vector of length N
    sz_vector = vector.shape  # 2d matrix that is N x M or M x N

    # define which way the sample is flipped relative to the input vector
    if sz_vector[0] == sz_sample[0]:
        case = 2
    elif sz_vector[0] == sz_sample[1]:
        case = 1
    else:
        emsg = ('Neither vector[0]==sample[0] nor vector[0]==sample[1] '
                '(function = {0})')
        print(emsg)
        raise ValueError(emsg.format(emsg))

    # we check if there are NaNs in the vector or the sample
    # if there are NaNs, we'll fit the rest of the domain
    isnan = (np.sum(np.isnan(vector)) != 0) or (np.sum(np.isnan(sample)) != 0)

    if case == 1:

        if isnan:
            # we create a mask of non-NaN
            keep = np.isfinite(vector) * np.isfinite(np.sum(sample, axis=0))
            # redefine the input vector to avoid NaNs
            vector = vector[keep]
            sample = sample[:, keep]

            sz_sample = sample.shape
            sz_vector = vector.shape

        # matrix of covariances
        mm = np.zeros([sz_sample[0], sz_sample[0]])
        # cross-terms of vector and columns of sample
        v = np.zeros(sz_sample[0])
        # reconstructed amplitudes
        amps = np.zeros(sz_sample[0])
        # reconstruted fit
        recon = np.zeros(sz_sample[1])

    if case == 2:
        # same as for case 1, but with axis flipped
        if isnan:
            # we create a mask of non-NaN
            keep = np.isfinite(vector) * np.isfinite(np.sum(sample, axis=1))
            vector = vector[keep]
            sample = sample[keep, :]

            sz_sample = sample.shape
            sz_vector = vector.shape

        mm = np.zeros([sz_sample[1], sz_sample[1]])
        v = np.zeros(sz_sample[1])
        amps = np.zeros(sz_sample[1])
        recon = np.zeros(sz_sample[0])

    # pass all variables and pre-formatted vectors to the @jit part of the code
    amp_out, recon_out = linear_minimization(vector, sample, mm, v, sz_sample, case,
                                             recon, amps)

    # if we had NaNs in the first place, we create a reconstructed vector
    # that has the same size as the input vector, but pad with NaNs values
    # for which we cannot derive a value
    if isnan:
        recon_out2 = np.zeros_like(keep) + np.nan
        recon_out2[keep] = recon_out
        recon_out = recon_out2

    return amp_out, recon_out


def linear_minimization(vector, sample, mm, v, sz_sample, case, recon, amps):
    # raise ValueError(emsg.format(func_name))
    # ​
    # vector of N elements
    # sample: matrix N * M each M column is adjusted in amplitude to minimize
    # the chi2 according to the input vector
    # output: vector of length M gives the amplitude of each column
    #
    if case == 1:
        # fill-in the co-variance matrix
        for i in range(sz_sample[0]):
            for j in range(i, sz_sample[0]):
                mm[i, j] = np.sum(sample[i, :] * sample[j, :])
                # we know the matrix is symetric, we fill the other half
                # of the diagonal directly
                mm[j, i] = mm[i, j]
            # dot-product of vector with sample columns
            v[i] = np.sum(vector * sample[i, :])
        # if the matrix cannot we inverted because the determinant is zero,
        # then we return a NaN for all outputs
        if np.linalg.det(mm) == 0:
            amps = np.zeros(sz_sample[0]) + np.nan
            recon = np.zeros_like(v)
            return amps, recon

        # invert coveriance matrix
        inv = np.linalg.inv(mm)
        # retrieve amplitudes
        for i in range(len(v)):
            for j in range(len(v)):
                amps[i] += inv[i, j] * v[j]

        # reconstruction of the best-fit from the input sample and derived
        # amplitudes
        for i in range(sz_sample[0]):
            recon += amps[i] * sample[i, :]
        return amps, recon

    if case == 2:
        # same as for case 1 but with axis flipped
        for i in range(sz_sample[1]):
            for j in range(i, sz_sample[1]):
                mm[i, j] = np.sum(sample[:, i] * sample[:, j])
                mm[j, i] = mm[i, j]
            v[i] = np.sum(vector * sample[:, i])

        if np.linalg.det(mm) == 0:
            return amps, recon

        inv = np.linalg.inv(mm)
        for i in range(len(v)):
            for j in range(len(v)):
                amps[i] += inv[i, j] * v[j]

        for i in range(sz_sample[1]):
            recon += amps[i] * sample[:, i]
        return amps, recon


def polyfit_odd_ratio(x, y, order, odd_ratio=1e-3):
    w_prev = np.zeros_like(x)
    w = np.ones_like(x)
    ite = 0

    while (np.abs(np.max(w - w_prev)) > odd_ratio ** 2) and (ite < 20):
        fit = np.polyfit(x, y, order, w=w)
        guess = np.polyval(fit, x)
        diff = y - guess
        err = sigma(diff)

        nsig = diff / err
        gg = np.exp(-0.5 * nsig ** 2)
        odd_bad = odd_ratio / (gg + odd_ratio)
        odd_good = 1 - odd_bad

        w_prev = np.array(w)
        w = odd_good

        ite += 1

    return fit


def lin_mini_errors(y0, yerr0, sample):
    y = y0.ravel()

    if sample.shape[1] == y.shape[0]:
        X = np.array(sample.T)
    else:
        X = np.array(sample)

    errs = yerr0.ravel()  # np.array(err[i]).ravel()

    # weights should be inverse of *square* error
    res_wls = sm.WLS(y, X, weights=1.0 / errs ** 2, missing='drop').fit()

    amps = res_wls.params
    errs = res_wls.bse

    recon = np.zeros_like(y0)
    for i in range(X.shape[1]):
        recon += amps[i] * X[:, i].reshape(y0.shape)

    return amps, errs, recon


def val_cheby(x, fit, domain=[0, 4088]):
    """
    :param x: x value for the y values with fit
    :param fit: output from fit_cheby
    :param domain: domain to be transformed to -1 -- 1. This is important to
    keep the components orthogonal. For SPIRou orders, the default is 0--4088.
    You *must* use the same domain when getting values with fit_cheby
    :return: corresponding y values to the x inputs
    """

    # transform to a -1 to 1 domain
    domain_cheby = 2 * (x - domain[0]) / (domain[1] - domain[0]) - 1

    return np.polynomial.chebyshev.chebval(domain_cheby, fit)


def hdr2wave(h, polytype='Smart'):
    """
    :param h: header of a SPIRou or NIRPS file
    :param polytype:
        polynomial type for the header. Newer data has Chebychev polynomial
        older (pre-Novembre 2022) uses standard polynomials

        Can have three values : "Cheby", "Standard" or "Smart". If the "smart"
        value is used (default), then we check the time of creation of file.
        Files created with versions after APERO 0.7.259 use 'Cheby'

    :return: wavelength map
    """

    if polytype not in ['Standard', 'Cheby', 'Smart']:
        raise ValueError('Polytype = "{}", but should be "Standard", "Cheby" or "Smart"'.format(polytype))

    if polytype == 'Smart':
        if h['VERSION'] < '0.7.259':
            polytype = 'Standard'
        else:
            polytype = 'Cheby'

    WAVEORDN = h['WAVEORDN']  # number of orders
    WAVEDEGN = h['WAVEDEGN']  # order of the fit

    # placeholder for coefficients
    coeffs = np.zeros([WAVEORDN, WAVEDEGN + 1])

    # number of x pixels along the order
    npix = 4088

    # index of pixels in grid
    index = np.arange(npix)

    # wavelength map
    wave = np.zeros([WAVEORDN, 4088])

    ii = 0
    for iord in range(WAVEORDN):
        for icoeff in range(WAVEDEGN + 1):
            key = 'WAVE' + str(ii).zfill(4)
            ii += 1
            coeffs[iord, icoeff] = h[key]

        if polytype == 'Standard':
            wave[iord] = np.polyval(coeffs[iord, :][::-1], index)
        else:
            wave[iord] = val_cheby(index, coeffs[iord, :], domain=[0, 4088])

    return wave


def simbad_teff(obj):
    Simbad.add_votable_fields('fe_h')
    result_table = Simbad.query_object(obj)

    if type(result_table) == type(None):
        return -1, ''

    if ('Fe_H_Teff' in result_table.colnames) and ('Fe_H_bibcode' in result_table.colnames):
        teff = result_table['Fe_H_Teff'].value.data[0]
        bibcode = result_table['Fe_H_bibcode'].value.data[0]

        if teff < 1:
            return -1, ''

        return teff, bibcode

    return -1, ''


def write_to_excel(outname, tbl):
    # create a pandas dataframe with all the info
    df = pd.DataFrame(dict(tbl))
    # save to excel file, first create a sheet with the dataframe
    sheetname = 'mySheet'
    with pd.ExcelWriter(outname) as writer:
        if not df.index.name:
            df.index.name = 'Index'
        df.to_excel(writer, sheet_name=sheetname)
    # then add a table to the sheet
    # open the workbook
    wb = openpyxl.load_workbook(filename=outname)
    tab = openpyxl.worksheet.table.Table(displayName="df", ref=f'A1:{chr(len(df.columns) + 64)}{len(df) + 1}')
    wb[sheetname].add_table(tab)
    wb.save(outname)

    return 0


def get_apero_astrometrics():
    """
    Get apero astrometrics

    Created on 2023-06-14
    Last updated 2023-06-14

    @author: cook
    """

    """
    Get the APERO astrometrics database as a table

    Deals with main and pending table

    Written by Neil Cook

    """
    # get the main table

    # =============================================================================
    # define variables
    # =============================================================================
    # url to download main table
    main_url = ('https://docs.google.com/spreadsheets/d/'
                '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/gviz/'
                'tq?tqx=out:csv&gid=0')
    # url to download pending table
    pend_url = ('https://docs.google.com/spreadsheets/d/'
                '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/gviz/'
                'tq?tqx=out:csv&gid=623506317')

    # object name column
    gl_objcol = 'OBJNAME'

    try:
        rawdata = requests.get(main_url)
        main_table = Table.read(rawdata.text, format='ascii')
    except Exception as _:
        raise ValueError('Cannot get table from {0}'.format(main_url))
    # get the pending table
    try:
        rawdata = requests.get(pend_url)
        pend_table = Table.read(rawdata.text, format='ascii')
    except Exception as _:
        raise ValueError('Cannot get table from {0}'.format(pend_url))
    # deal with overlap in pending table
    for _table in [pend_table]:
        # only do this if this table has some entries
        if len(_table) != 0:
            # make sure we have the object name column
            if gl_objcol in _table.colnames:
                pmask = ~np.in1d(_table[gl_objcol], main_table[gl_objcol])
                # add new columns to main table
                main_table = vstack([main_table, _table[pmask]])
    # return the table
    return main_table
