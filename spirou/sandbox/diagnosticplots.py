#!/usr/bin/env python
# coding: utf-8
"""
Plots for data diagnostics

@author: CharlesCadieux
"""
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PyPDF2 import PdfFileMerger
from weasyprint import HTML
import numpy as np
from astropy.table import Table, Column
from astropy.timeseries import LombScargle
#from PyAstronomy.pyasl import foldAt
from datetime import date
from astropy.io import ascii
import requests
from astropy import units as u
from astropy.coordinates import Angle
today = date.today()

def find_nearest(array, value):
    # Find the index position of the nearest value in an array

    index = (np.abs(array - value)).argmin()
    return index

def semiamp(P, Mp, Rp, Ms, i, e):
    # P (days), Mp (M_earth), Rp (R_earth), Ms (M_sun), i (deg)

    G = 6.67408e-11 # Gravitaional constant
    P *= (24 * 60 * 60) # days to sec
    Mp *= 5.972e24 # M_earth to kg
    Rp *= 6.371e6 # R_earth to m
    Ms *= 1.9885e30 # M_sun to kg

    K = (2 * np.pi * G / P) ** (1 / 3) * Mp * np.sin(np.radians(i)) / (Ms +
         Mp) ** (2 / 3) / np.sqrt(1 - e ** 2) # m/s

    return np.around(K, 2) # 2 decimals

def mass_est(Rp):
    # Chen and Kipping (2017) as implemented by Louie et al. (2018)
    # Rp (R_earth)

    if Rp < 1.23:
        Mp = 0.9718 * Rp**3.58 # Mp (M_earth)
        return np.around(Mp, 2) # 2 decimals
    elif 1.23 <= Rp < 14.26:
        Mp = 1.436 * Rp**1.70 # Mp (M_earth)
        return np.around(Mp, 2) # 2 decimals
    else:
        return np.nan # Degeneracy R vs M

def radius_est(Mp):
    # Chen and Kipping (2017)
    # Mp (M_earth)

    if Mp < 2.04:
        # Terrestrial
        Rp = 1.008 * Mp**0.279 #Rp (R_earth)
        return np.around(Rp, 2) # 2 decimals
    elif 2.04 <= Mp < 131.6:
        # Neptunian
        Rp = 0.808 * Mp**0.589 #Rp (R_earth)
        return np.around(Rp, 2) # 2 decimals
    elif 131.6 <= Mp < 26637.2:
        # Jovian
        Rp = 17.394 * Mp**(-0.04) #Rp (R_earth)
        return np.around(Rp, 2) # 2 decimals
    else:
        # Stellar
        Rp = 1.476e-3 * Mp**0.88 #Rp (R_earth)
        return np.around(Rp, 2) # 2 decimals

def bulkrho(Mp, Rp):
    # Mp (M_earth), Rp (R_earth)

    Mp *= 5.972e24 # M_earth to kg
    Rp *= 6.371e6 # R_earth to m

    rho = Mp / (4 * np.pi / 3 * Rp**3) / 1000 # Density in  g/cm3

    return np.around(rho , 2) # 2 decimals

def atmosig(Teq, Rs, rho, mu):
    # Teq (K), Rs (R_sun), rho (g/cm3), mu (amu)

    G = 6.67408E-11 # Gravitaional constant
    kb = 1.38064852E-23 # Boltzmann constant
    Rs *= 6.956E8 #R_sun to m
    rho *= 1000 #g/cm3 to kg/m3
    mu *= 1.66054E-27 #amu to kg

    atmo = (10 * kb * 3 / 4 / np.pi / G) * Teq / Rs**2 / rho / mu * 1E6 # ppm
    return np.around(atmo, 2) # 2 decimals

def TSM(Rp, Mp, Rs, Teff, a, Jmag):
    # Transmission Spectroscopy Metric (Kempton et al. 2018)
    # Rp (R_earth), Mp (M_earth), Rs (R_sun), Teff (K), a (au), Jmag (mag)

    if Rp < 1.5 : scale_fac = 0.190
    elif 1.5 <= Rp < 2.75 : scale_fac = 1.26
    elif 2.75 <= Rp < 4.0 : scale_fac = 1.28
    elif 4.0 <= Rp < 10.0 : scale_fac = 1.15
    else : scale_fac = 1.0

    Rs_m = Rs * 6.96340e8 # R_sun to m
    a_m = a * 1.496e+11 # au to m
    Teq = Teff * np.sqrt(Rs_m / a_m) * (0.25)**(0.25) # zero albedo Teq

    tsm = scale_fac * Rp**3 * Teq / Mp / Rs**2 * 10**(-Jmag/5) # TSM

    return np.around(tsm, 2) # 2 decimals

def Planck(wav, T):
    # Planck function
    # wav (m), T (K)

    h = 6.62607004e-34 # Planck constant
    c = 2.99792458e+8 # speed of light
    kb = 1.38064852E-23 # Bolztmann constant
    a = 2.0 * h * c**2
    b = h * c / (wav * kb * T)
    intensity = a / ((wav**5) * (np.exp(b) - 1.0))

    return intensity

def ESM(Rp, Rs, Teff, a, Kmag):
    # Emission Spectroscopy Metric (Kempton et al. 2018)
    # Rp (R_earth), Rs (R_sun), Teff (K), a (au), Kmag (mag)

    Rp_m = Rp * 6.371e6 # R_earth to m
    Rs_m = Rs * 6.96340e8 # R_sun to m
    a_m = a * 1.496e+11 # au to m
    Teq = Teff * np.sqrt(Rs_m / a_m) * (0.25)**(0.25) # zero albedo Teq
    Tday = 1.10 * Teq # Day side temperature estimate

    esm = 4.29e6 * Planck(7.5e-6, Tday)/Planck(7.5e-6,Teff) *(Rp_m /
          Rs_m)**2 * 10**(-Kmag/5)

    return np.around(esm, 2) # 2 decimals


def system_info(obj_sci, exoarchive_name = None):
    """
    Automatic system parameters query from Exoplanet Archive
    (exoplanetarchive.ipac.caltech.edu/) or ExoFOP
    (exofop.ipac.caltech.edu/tess/index.php).

    obj_sci: Name of the object (host star).
    exoarchive_name: Name of the object on the Exoplanet Archive.

    Produces a summary page pdf of the stellar and planetary parameters.
    Return the name of the planet, the time of mid-transit or periastron and
    the orbital period.
    """

    if exoarchive_name is None:
        exoarchive_name = obj_sci


    # Query from the masterfile (https://github.com/AntoineDarveau/masterfile)
    info = requests.get('http://www.astro.umontreal.ca/~adb/masterfile.ecsv')
    info = Table.read(info.text, format = 'ascii.ecsv') # Exoplanet Archive data
    info['pl_trandur'] *= 24 # days to hours
    info['pl_trandep'] *= 10 # percent to ppt
    info['st_lum'] = np.around(10**info['st_lum'], 4) # log L_sun to L_sun

    # For a TOI object, look on ExoFOP
    if obj_sci[:3] == 'TOI':
        infoTOI = requests.get(
            "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv")
        names = infoTOI.text[0:824]
        names = list(names.split(","))
        infoTOI = Table.read(
                infoTOI.text, format = "ascii", names = names, data_start = 1)
        index, = np.where(
                np.array(infoTOI['TOI'], dtype = int) == int(obj_sci[4:]))
        nbplanet = len(index)

        # Object not found
        if nbplanet == 0:
            print('No Object Found')
            return [], [], []

        # Loop through planets
        for i in range(nbplanet):
            TICid = infoTOI['TIC ID'][index[i]]
            # Look on the TIC individual pages
            ExoFop = requests.get(
                "https://exofop.ipac.caltech.edu/tess/download_target.php?id={0}".format(TICid))

            # Stellar Mass
            str_st1 = ExoFop.text.find("STELLAR PARMAETERS")
            str_st2 = ExoFop.text.find("MAGNITUDES")
            st = ExoFop.text[str_st1:str_st2:]
            st = Table.read(st, format = 'ascii.fixed_width',
                    col_starts=(0, 25, 43, 65, 87, 109, 131, 153, 175, 197, 219,
                    241, 263, 285, 307, 329, 351, 373, 395, 417, 439, 461, 483,
                    505, 530, 552, 574, 599, 621, 643, 665, 687, 709, 731, 749,
                    767), header_start = 1, data_start = 2)

            # Magnitudes
            str_mag1 = ExoFop.text.find("MAGNITUDES")
            str_mag2 = ExoFop.text.find("IMAGING OBSERVATIONS")
            mag = ExoFop.text[str_mag1:str_mag2:]
            mag = Table.read(mag, format = 'ascii.fixed_width',
                    col_starts=(0, 18, 36, 54, 79, 99, 117, 135),
                    header_start = 1, data_start = 2)
            indexV, = np.where(mag['Band'] == 'V')
            indexG, = np.where(mag['Band'] == 'Gaia')
            indexJ, = np.where(mag['Band'] == 'J')
            indexH, = np.where(mag['Band'] == 'H')
            indexK, = np.where(mag['Band'] == 'K')

            # Mass estimate
            masse = mass_est(infoTOI['Planet Radius (R_Earth)'][index[i]])

            # Bulk density estimate
            rho = bulkrho(masse.copy(),
                          infoTOI['Planet Radius (R_Earth)'][index[i]])

            # RV semi-amplitude estimate
            K = semiamp(infoTOI['Period (days)'][index[i]],
                        masse.copy(),
                        infoTOI['Planet Radius (R_Earth)'][index[i]],
                        float(st['Mass (M_Sun)'][0]),
                        90,
                        0)

            # Add row to info table
            info.add_row()
            info[-1]['pl_hostname'] = 'TIC {0}'.format(infoTOI['TIC ID'][index[i]])
            info[-1]['ra_str'] = str(infoTOI["RA"][index[i]])
            info[-1]['dec_str'] = str(infoTOI["Dec"][index[i]])
            info[-1]['ra'] = Angle(infoTOI["RA"][index[i]], unit=u.hour).deg
            info[-1]['dec'] = Angle(infoTOI["Dec"][index[i]], unit = u.deg).deg
            info[-1]['st_rad'] = infoTOI['Stellar Radius (R_Sun)'][index[i]]
            info[-1]['st_mass'] = st['Mass (M_Sun)'][0]
            info[-1]['st_lum'] = st['Luminosity'][0]
            info[-1]['st_teff'] = infoTOI['Stellar Eff Temp (K)'][index[i]]
            info[-1]['st_spstr'] = '0'
            info[-1]['st_dist'] = infoTOI['Stellar Distance (pc)'][index[i]]
            if len(indexV) == 1 : info[-1]['st_vj'] = mag['Value'][indexV]
            if len(indexG) == 1 : info[-1]['gaia_gmag'] = mag['Value'][indexG]
            if len(indexJ) == 1 : info[-1]['st_j'] = mag['Value'][indexJ]
            if len(indexH) == 1 : info[-1]['st_h'] = mag['Value'][indexH]
            if len(indexK) == 1 : info[-1]['st_k'] = mag['Value'][indexK]
            info[-1]['pl_name'] = 'TOI {0}'.format(infoTOI['TOI'][index[i]])
            info[-1]['pl_bmassprov'] = 'Mass Estimate'
            info[-1]['pl_masse'] = masse
            info[-1]['pl_rade'] = infoTOI['Planet Radius (R_Earth)'][index[i]]
            info[-1]['pl_dens'] = rho
            info[-1]['pl_tranmid'] = infoTOI['Epoch (BJD)'][index[i]]
            info[-1]['pl_orbper'] = infoTOI['Period (days)'][index[i]]
            info[-1]['pl_trandur'] = infoTOI['Duration (hours)'][index[i]]
            info[-1]['pl_orbsmax'] = 0
            info[-1]['pl_eqt'] = infoTOI['Planet Equil Temp (K)'][index[i]]
            info[-1]['pl_insol'] = infoTOI['Planet Insolation (Earth Flux)'][index[i]]
            info[-1]['pl_rvamp'] = K
            info[-1]['pl_masselim'] = 0
            info[-1]['pl_tranflag'] = 1
            info[-1]['pl_trandep'] = infoTOI['Depth (ppm)'][index[i]] / 1000
            info[-1]['pl_facility'] = 'Transiting Exoplanet Survey Satellite (TESS)'

        info = info[-nbplanet:] # Select only the TOI planets
        mass_colors = nbplanet * ['#ff0000']
        rho_colors = nbplanet * ['#ff0000']
        K_colors = nbplanet * ['#ff0000']

    # Look on Exoplanet Archive
    else:
        info = info[info['pl_hostname'] == '{0}'.format(exoarchive_name)]
        nbplanet = len(info)

        # Object not found
        if nbplanet == 0:
            print('No Object Found')
            return [], [], []

        # Radius estimate for non-transiting planet
        info['pl_masse'].mask = False
        info['pl_rade'].mask = False
        # Check where the mass is known, but not the radius
        index = np.intersect1d(
                (info['pl_masse'].data > 0).nonzero(),
                (info['pl_rade'].data == 0).nonzero())
        for i in index:
            info['pl_rade'][i] = radius_est(info['pl_masse'][i])
        # Black if known, red if estimated
        mass_colors = nbplanet * ['#000000']
        rade_colors = nbplanet * ['#000000']
        mass_colors = np.array(mass_colors)
        rade_colors = np.array(rade_colors)
        rade_colors[index] = '#ff0000'

        # Mass estimate for transiting planet
        info['pl_masse'].mask = False
        info['pl_rade'].mask = False
        # Check where the radius is known, but not the mass
        index = np.intersect1d(
                (info['pl_masse'].data == 0).nonzero(),
                (info['pl_rade'].data > 0).nonzero())
        for i in index:
            info['pl_masse'][i] = mass_est(info['pl_rade'][i])
        # Black if known, red if estimated
        mass_colors[index] = '#ff0000'


    # Fill in missing spectral type
    # Teff to SpT prescription from Table 5 of Pecaut & Mamajek (2013)
    spectral = requests.get(
    "http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt")
    sptype = Table.read(spectral.text, format = 'ascii.csv', delimiter = ' ',
             comment = '&', header_start = 22, data_start = 23, data_end = 141)
    del sptype['#SpT_1']
    sptype.rename_column('#SpT', 'SpT')
    info['st_teff'].mask = False
    info['st_spstr'].mask = False
    index, = (info['st_spstr'].data == '0').nonzero()
    index.sort()
    spt_colors = nbplanet * ['#000000']
    spt_colors = np.array(spt_colors)
    for i in index:
        if info['st_teff'][i] > 0:
            info['st_spstr'][i] = sptype['SpT'][find_nearest(sptype['Teff'],
                                                info['st_teff'][i])]
            spt_colors[i] = '#ff0000'
        else:
            continue
    index, = np.where(info['st_spstr'].mask == 1)
    info['st_spstr'][index] = 0


    # Fill in missing semi-major axis
    info['pl_orbsmax'].mask = False
    sma = np.around((info['st_mass'].data *
                    (info['pl_orbper'].data / 365.25)**2)**(1/3), 4)
    index, = np.logical_or(np.isnan(info['pl_orbsmax']),
                           info['pl_orbsmax'] == 0).nonzero()
    info['pl_orbsmax'][index] = sma[index]
    semima_colors = nbplanet * ['#000000']
    semima_colors = np.array(semima_colors)
    semima_colors[index] = '#ff0000'
    index, = np.where(info['pl_orbsmax'].mask == 1)
    info['pl_orbsmax'][index] = 0


    # Fill in missing Teq
    info['pl_eqt'].mask = False
    teq = np.around(info['st_teff'].data * np.sqrt(info['st_rad'].data / 215
            / 2 / info['pl_orbsmax'].data) * (1. - 0.3) **0.25, 2) # 0.3 albedo
    index, = np.logical_or(np.isnan(info['pl_eqt']),
                           info['pl_eqt'] == 0).nonzero()
    info['pl_eqt'][index] = teq[index]
    teq_colors = nbplanet * ['#000000']
    teq_colors = np.array(teq_colors)
    teq_colors[index] = '#ff0000'
    index, = np.where(info['pl_eqt'].mask == 1)
    info['pl_eqt'][index] = 0


    # Fill in missing planet density
    info['pl_dens'].mask = False
    index, = np.logical_or(np.isnan(info['pl_dens']),
                           info['pl_dens'] == 0).nonzero()
    density = bulkrho(info['pl_masse'].data.copy(), info['pl_rade'].data.copy())
    info['pl_dens'][index] = density[index]
    rho_colors = nbplanet * ['#000000']
    rho_colors = np.array(rho_colors)
    rho_colors[index] = '#ff0000'
    index, = np.where(info['pl_dens'].mask == 1)
    info['pl_dens'][index] = 0


    # Fill in missing insolation
    info['pl_insol'].mask = False
    index, = np.logical_or(np.isnan(info['pl_insol']),
                           info['pl_insol'] == 0).nonzero()
    insol = np.around(info['st_rad']**2 * (info['st_teff']
                      / 5777)**4 / info['pl_orbsmax']**2, 2)
    info['pl_insol'][index] = insol[index]
    insol_colors = nbplanet * ['#000000']
    insol_colors = np.array(insol_colors)
    insol_colors[index] = '#ff0000'
    index, = np.where(info['pl_insol'].mask == 1)
    info['pl_insol'][index] = 0


    # Fill in missing K
    info['pl_rvamp'].mask = False
    index, = np.logical_or(np.isnan(info['pl_rvamp']),
                           info['pl_rvamp'] == 0).nonzero()
    K = semiamp(info['pl_orbper'].data.copy(),
                info['pl_masse'].data.copy(),
                info['pl_rade'].data.copy(),
                info['st_mass'].data.copy(),
                90,
                0)
    info['pl_rvamp'][index] = K[index]
    K_colors = nbplanet * ['#000000']
    K_colors = np.array(K_colors)
    K_colors[index] = '#ff0000'
    index, = np.where(info['pl_rvamp'].mask == 1)
    info['pl_rvamp'][index] = 0


    # Create atmospheric signal column
    atmosigcolumn = Column(name = 'pl_atmosig',
                           data = np.zeros(nbplanet))
    info.add_column(atmosigcolumn, 21)
    muarray = np.zeros(nbplanet)
    for i in range(nbplanet):
        if info['pl_rade'][i] <= 2:
            muarray[i] = 28.97 # Earth atmosphere mu
        elif 2 < info['pl_rade'][i] <= 6:
            muarray[i] = 2.61 # Neptune atmosphere mu
        elif info['pl_rade'][i] > 6:
            muarray[i] = 2.22 # Jupiter atmosphere mu
    info['pl_atmosig'] = atmosig(info['pl_eqt'].data.copy(),
                                 info['st_rad'].data.copy(),
                                 info['pl_dens'].data.copy(),
                                 muarray.copy())
    red_colors = nbplanet * ['#ff0000']
    red_colors = np.array(red_colors)
    index, = np.where(info['pl_atmosig'].mask == 1)
    info['pl_atmosig'][index] = 0

    # Create a TSM column
    tsmcolumn = Column(name = 'pl_tsm',
                       data = np.zeros(nbplanet))
    info.add_column(tsmcolumn, 21)
    for i in range(nbplanet):
        info['pl_tsm'][i] = TSM(info['pl_rade'][i],
                                info['pl_masse'][i],
                                info['st_rad'][i],
                                info['st_teff'][i],
                                info['pl_orbsmax'][i],
                                info['st_j'][i])


    # Create a ESM column
    esmcolumn = Column(name = 'pl_esm',
                       data = np.zeros(nbplanet))
    info.add_column(esmcolumn, 21)
    info['pl_esm'] = ESM(info['pl_rade'].data.copy(),
                         info['st_rad'].data.copy(),
                         info['st_teff'].data.copy(),
                         info['pl_orbsmax'].data.copy(),
                         info['st_k'].data.copy())


    # Fill in missing transit depth
    info['pl_trandep'].mask = False
    index, = np.logical_or(np.isnan(info['pl_trandep']),
                           info['pl_trandep'] == 0).nonzero()
    depth = np.around((6.371e6 / 6.956e8)**2 * (info['pl_rade']
                      / info['st_rad'])**2 * 1000, 2) # ppt
    info['pl_trandep'][index] = depth[index]
    depth_colors = nbplanet * ['#000000']
    depth_colors = np.array(depth_colors)
    depth_colors[index] = '#ff0000'
    index, = np.where(info['pl_trandep'].mask == 1)
    info['pl_trandep'][index] = 0


    # Fill in missing luminosity
    info['st_lum'].mask = False
    index, = np.logical_or(np.isnan(info['st_lum']),
                           info['st_lum'] == 0).nonzero()
    lum = np.around((info['st_rad'].data)**2 * (
                     info['st_teff'].data / 5777)**4, 4)
    info['st_lum'][index] = lum[index]
    lum_colors = nbplanet * ['#000000']
    lum_colors = np.array(lum_colors)
    lum_colors[index] = '#ff0000'
    index, = np.where(info['st_lum'].mask == 1)
    info['st_lum'][index] = 0


    def findNA(x):
        # Formating: Leave a blank space if parameter value is
        # non available (NA)
        index, = np.logical_or(np.isnan(x), x == 0).nonzero()
        listx = list(x)
        for i in index:
            listx[i] = ' '
        return listx

    stname = info['pl_hostname'].data # string
    name = info['pl_name'].data # string
    ra_str = info['ra_str'].data # string
    dec_str = info['dec_str'].data # string
    ra = findNA(info['ra'].data)
    dec = findNA(info['dec'].data)
    stmass = findNA(info['st_mass'].data)
    stradius = findNA(info['st_rad'].data)
    stlum = findNA(info['st_lum'].data)
    teff = findNA(info['st_teff'].data)
    spt = info['st_spstr'].data # string
    dist = findNA(info['st_dist'].data)
    mass = []
    for i in range(nbplanet):
        if info['pl_masse'][i] > 0:
            mass.append(info['pl_masse'][i])
        else:
            mass.append(info['pl_msinie'][i])
    mass = findNA(np.array(mass))
    rade = findNA(info['pl_rade'].data)
    rho = findNA(info['pl_dens'].data)
    tp = findNA(info['pl_orbtper'].data)
    t0 = findNA(info['pl_tranmid'].data)
    period = findNA(info['pl_orbper'].data)
    depth = findNA(info['pl_trandep'].data)
    dur = findNA(info['pl_trandur'].data)
    semima = findNA(info['pl_orbsmax'].data)
    teq = findNA(info['pl_eqt'].data)
    insol = findNA(info['pl_insol'].data)
    K = findNA(info['pl_rvamp'].data)
    atmosignal = findNA(info['pl_atmosig'].data)
    tsm = findNA(info['pl_tsm'].data)
    esm =  findNA(info['pl_esm'].data)
    vmag = findNA(info["st_vj"].data)
    gmag = findNA(info["gaia_gmag"].data)
    jmag = findNA(info["st_j"].data)
    hmag = findNA(info["st_h"].data)
    kmag = findNA(info["st_k"].data)


    for i in range(nbplanet):
    # 1 HTML document per planet
        if obj_sci[:3] == 'TOI':
            html_text = f"""
            <!DOCTYPE html>
            <html>
              <head>
                <meta charset="UTF-8" />
                <style>
                @page {{
                    size: 10in 8in;
                    margin: 1in 1in 1in 1in;
             }}
            </style>
              </head>
            <body id="main">
            <table>
              <tr>
                <th colspan="2">Stellar Parameters</th>
                <th colspan="2">Planet Parameters</th>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Host name</span></td>
                <td>{stname[i]}</td>
                <td><span style="color: #2874a6;">Name</span></td>
                <td>{name[i]}</td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">R.A.</span></td>
                <td>{ra_str[i]}</td>
                <td><span style="color: #2874a6;">Mass (M⊕)</span></td>
                <td><span style="color: {mass_colors[i]};">{mass[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Dec.</span></td>
                <td>{dec_str[i]}</td>
                <td><span style="color: #2874a6;">Radius (R⊕)</span></td>
                <td>{rade[i]}</td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Mass (M⊙)</span></td>
                <td>{stmass[i]}</td>
                <td><span style="color: #2874a6;">Density (g/cm³)</span></td>
                <td><span style="color: {rho_colors[i]};">{rho[i]}</span></td>
              </tr>
                <tr>
                <td><span style="color: #2874a6;">Radius (R⊙)</span></td>
                <td>{stradius[i]}</td>
                <td><span style="color: #2874a6;">Period (days)</span></td>
                <td>{period[i]}</td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Lum. (L☉)</span></td>
                <td><span style="color: {lum_colors[i]};">{stlum[i]}</span></td>
                <td><span style="color: #2874a6;">Transit depth (ppt)</span></td>
                <td><span style="color: {depth_colors[i]};">{depth[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Teff. (K)</span></td>
                <td>{teff[i]}</td>
                <td><span style="color: #2874a6;">Transit dur. (hr)</span></td>
                <td>{dur[i]}</td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Spect. type</span></td>
                <td><span style="color: {spt_colors[i]};">{spt[i]}</span></td>
                <td><span style="color: #2874a6;">Orb. dist. (AU)</span></td>
                <td><span style="color: {semima_colors[i]};">{semima[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Dist. (pc)</span></td>
                <td>{dist[i]}</td>
                <td><span style="color: #2874a6;">Teq. (K)</span></td>
                <td><span style="color: {teq_colors[i]};">{teq[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">V</span></td>
                <td>{vmag[i]}</td>
                <td><span style="color: #2874a6;">Insol. (Earth flux)</span></td>
                <td><span style="color: {insol_colors[i]};">{insol[i]}</span></td>
              </tr>
                <tr>
                <td><span style="color: #2874a6;">Gaia</span></td>
                <td>{gmag[i]}</td>
                <td><span style="color: #2874a6;">RV K (m/s)</span></td>
                <td><span style="color: {K_colors[i]};">{K[i]}</span></td>

              </tr>
              <tr>
                <td><span style="color: #2874a6;">J</span></td>
                <td>{jmag[i]}</td>
                <td><span style="color: #2874a6;">Atmo. sig. (ppm)</span></td>
                <td><span style="color: {red_colors[i]};">{atmosignal[i]}</span></td>
              </tr>
                <tr>
                <td><span style="color: #2874a6;">H</span></td>
                <td>{hmag[i]}</td>
                <td><span style="color: #2874a6;">TSM</span></td>
                <td><span style="color: {red_colors[i]};">{tsm[i]}</span></td>

              </tr>
                <tr>
                <td><span style="color: #2874a6;">K</span></td>
                <td>{kmag[i]}</td>
                <td><span style="color: #2874a6;">ESM</span></td>
                <td><span style="color: {red_colors[i]};">{esm[i]}</span></td>
              </tr>
            </table>
            <br>
            <a href="https://exofop.ipac.caltech.edu/tess/target.php?id={TICid}">Link to ExoFOP</a><br>
            </body>
            </html>
            """
        elif info['pl_tranflag'][i] == 1:
            html_text = f"""
            <!DOCTYPE html>
            <html>
              <head>
                <meta charset="UTF-8" />
                <style>
                @page {{
                    size: 10in 8in;
                    margin: 1in 1in 1in 1in;
             }}
            </style>
              </head>
            <body id="main">
            <table>
              <tr>
                <th colspan="2">Stellar Parameters</th>
                <th colspan="2">Planet Parameters</th>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Host name</span></td>
                <td>{stname[i]}</td>
                <td><span style="color: #2874a6;">Name</span></td>
                <td>{name[i]}</td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">R.A.</span></td>
                <td>{ra_str[i]}</td>
                <td><span style="color: #2874a6;">Mass (M⊕)</span></td>
                <td><span style="color: {mass_colors[i]};">{mass[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Dec.</span></td>
                <td>{dec_str[i]}</td>
                <td><span style="color: #2874a6;">Radius (R⊕)</span></td>
                <td>{rade[i]}</td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Mass (M⊙)</span></td>
                <td>{stmass[i]}</td>
                <td><span style="color: #2874a6;">Density (g/cm³)</span></td>
                <td><span style="color: {rho_colors[i]};">{rho[i]}</span></td>
              </tr>
                <tr>
                <td><span style="color: #2874a6;">Radius (R⊙)</span></td>
                <td>{stradius[i]}</td>
                <td><span style="color: #2874a6;">Period (days)</span></td>
                <td>{period[i]}</td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Lum. (L☉)</span></td>
                <td><span style="color: {lum_colors[i]};">{stlum[i]}</span></td>
                <td><span style="color: #2874a6;">Transit depth (ppt)</span></td>
                <td><span style="color: {depth_colors[i]};">{depth[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Teff. (K)</span></td>
                <td>{teff[i]}</td>
                <td><span style="color: #2874a6;">Transit dur. (hr)</span></td>
                <td>{dur[i]}</td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Spect. type</span></td>
                <td><span style="color: {spt_colors[i]};">{spt[i]}</span></td>
                <td><span style="color: #2874a6;">Orb. dist. (AU)</span></td>
                <td><span style="color: {semima_colors[i]};">{semima[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Dist. (pc)</span></td>
                <td>{dist[i]}</td>
                <td><span style="color: #2874a6;">Teq. (K)</span></td>
                <td><span style="color: {teq_colors[i]};">{teq[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">V</span></td>
                <td>{vmag[i]}</td>
                <td><span style="color: #2874a6;">Insol. (Earth flux)</span></td>
                <td><span style="color: {insol_colors[i]};">{insol[i]}</span></td>
              </tr>
                <tr>
                <td><span style="color: #2874a6;">Gaia</span></td>
                <td>{gmag[i]}</td>
                <td><span style="color: #2874a6;">RV K (m/s)</span></td>
                <td><span style="color: {K_colors[i]};">{K[i]}</span></td>

              </tr>
              <tr>
                <td><span style="color: #2874a6;">J</span></td>
                <td>{jmag[i]}</td>
                <td><span style="color: #2874a6;">Atmo. sig. (ppm)</span></td>
                <td><span style="color: {red_colors[i]};">{atmosignal[i]}</span></td>
              </tr>
                <tr>
                <td><span style="color: #2874a6;">H</span></td>
                <td>{hmag[i]}</td>
                <td><span style="color: #2874a6;">TSM</span></td>
                <td><span style="color: {red_colors[i]};">{tsm[i]}</span></td>

              </tr>
                <tr>
                <td><span style="color: #2874a6;">K</span></td>
                <td>{kmag[i]}</td>
                <td><span style="color: #2874a6;">ESM</span></td>
                <td><span style="color: {red_colors[i]};">{esm[i]}</span></td>
              </tr>
            </table>
            <br>
            <a href="https://exoplanetarchive.ipac.caltech.edu/overview/{obj_sci}">Link to NASA Exoplanet Archive</a><br>
            <a href="http://simbad.u-strasbg.fr/simbad/sim-id?Ident={obj_sci}">Link to Simbad</a>
            </body>
            </html>
            """
        elif info['pl_masse'][i] > 0:
            html_text = f"""
            <!DOCTYPE html>
            <html>
              <head>
                <meta charset="UTF-8" />
                <style>
                @page {{
                    size: 10in 8in;
                    margin: 1in 1in 1in 1in;
             }}
            </style>
              </head>
            <body id="main">
            <table>
              <tr>
                <th colspan="2">Stellar Parameters</th>
                <th colspan="2">Planet Parameters</th>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Host name</span></td>
                <td>{stname[i]}</td>
                <td><span style="color: #2874a6;">Name</span></td>
                <td>{name[i]}</td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">R.A.</span></td>
                <td>{ra_str[i]}</td>
                <td><span style="color: #2874a6;">Mass (M⊕)</span></td>
                <td><span style="color: {mass_colors[i]};">{mass[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Dec.</span></td>
                <td>{dec_str[i]}</td>
                <td><span style="color: #2874a6;">Radius (R⊕)</span></td>
                <td><span style="color: {rade_colors[i]};">{rade[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Mass (M⊙)</span></td>
                <td>{stmass[i]}</td>
                <td><span style="color: #2874a6;">Density (g/cm³)</span></td>
                <td><span style="color: {rho_colors[i]};">{rho[i]}</span></td>
              </tr>
                <tr>
                <td><span style="color: #2874a6;">Radius (R⊙)</span></td>
                <td>{stradius[i]}</td>
                <td><span style="color: #2874a6;">Period (days)</span></td>
                <td>{period[i]}</td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Lum. (L☉)</span></td>
                <td><span style="color: {lum_colors[i]};">{stlum[i]}</span></td>
                <td><span style="color: #2874a6;">Orb. dist. (AU)</span></td>
                <td><span style="color: {semima_colors[i]};">{semima[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Teff. (K)</span></td>
                <td>{teff[i]}</td>
                <td><span style="color: #2874a6;">Teq. (K)</span></td>
                <td><span style="color: {teq_colors[i]};">{teq[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Spect. type</span></td>
                <td><span style="color: {spt_colors[i]};">{spt[i]}</span></td>
                <td><span style="color: #2874a6;">Insol. (Earth flux)</span></td>
                <td><span style="color: {insol_colors[i]};">{insol[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Dist. (pc)</span></td>
                <td>{dist[i]}</td>
                <td><span style="color: #2874a6;">RV K (m/s)</span></td>
                <td><span style="color: {K_colors[i]};">{K[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">V</span></td>
                <td>{vmag[i]}</td>
                <td></td>
                <td></td>
              </tr>
                <tr>
                <td><span style="color: #2874a6;">Gaia</span></td>
                <td>{gmag[i]}</td>
                <td></td>
                <td></td>

              </tr>
              <tr>
                <td><span style="color: #2874a6;">J</span></td>
                <td>{jmag[i]}</td>
                <td></td>
                <td></td>
              </tr>
                <tr>
                <td><span style="color: #2874a6;">H</span></td>
                <td>{hmag[i]}</td>
                <td></td>
                <td></td>

              </tr>
                <tr>
                <td><span style="color: #2874a6;">K</span></td>
                <td>{kmag[i]}</td>
                <td></td>
                <td></td>
              </tr>
            </table>
            <br>
            <a href="https://exoplanetarchive.ipac.caltech.edu/overview/{obj_sci}">Link to NASA Exoplanet Archive</a><br>
            <a href="http://simbad.u-strasbg.fr/simbad/sim-id?Ident={obj_sci}">Link to Simbad</a>
            </body>
            </html>
            """
        else:
            html_text = f"""
            <!DOCTYPE html>
            <html>
              <head>
                <meta charset="UTF-8" />
                <style>
                @page {{
                    size: 10in 8in;
                    margin: 1in 1in 1in 1in;
             }}
            </style>
              </head>
            <body id="main">
            <table>
              <tr>
                <th colspan="2">Stellar Parameters</th>
                <th colspan="2">Planet Parameters</th>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Host name</span></td>
                <td>{stname[i]}</td>
                <td><span style="color: #2874a6;">Name</span></td>
                <td>{name[i]}</td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">R.A.</span></td>
                <td>{ra_str[i]}</td>
                <td><span style="color: #2874a6;">M sin i (M⊕)</span></td>
                <td><span style="color: {mass_colors[i]};">{mass[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Dec.</span></td>
                <td>{dec_str[i]}</td>
                <td><span style="color: #2874a6;">Radius (R⊕)</span></td>
                <td>{rade[i]}</td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Mass (M⊙)</span></td>
                <td>{stmass[i]}</td>
                <td><span style="color: #2874a6;">Density (g/cm³)</span></td>
                <td><span style="color: {rho_colors[i]};">{rho[i]}</span></td>
              </tr>
                <tr>
                <td><span style="color: #2874a6;">Radius (R⊙)</span></td>
                <td>{stradius[i]}</td>
                <td><span style="color: #2874a6;">Period (days)</span></td>
                <td>{period[i]}</td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Lum. (L☉)</span></td>
                <td><span style="color: {lum_colors[i]};">{stlum[i]}</span></td>
                <td><span style="color: #2874a6;">Orb. dist. (AU)</span></td>
                <td><span style="color: {semima_colors[i]};">{semima[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Teff. (K)</span></td>
                <td>{teff[i]}</td>
                <td><span style="color: #2874a6;">Teq. (K)</span></td>
                <td><span style="color: {teq_colors[i]};">{teq[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Spect. type</span></td>
                <td><span style="color: {spt_colors[i]};">{spt[i]}</span></td>
                <td><span style="color: #2874a6;">Insol. (Earth flux)</span></td>
                <td><span style="color: {insol_colors[i]};">{insol[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">Dist. (pc)</span></td>
                <td>{dist[i]}</td>
                <td><span style="color: #2874a6;">RV K (m/s)</span></td>
                <td><span style="color: {K_colors[i]};">{K[i]}</span></td>
              </tr>
              <tr>
                <td><span style="color: #2874a6;">V</span></td>
                <td>{vmag[i]}</td>
                <td></td>
                <td></td>
              </tr>
                <tr>
                <td><span style="color: #2874a6;">Gaia</span></td>
                <td>{gmag[i]}</td>
                <td></td>
                <td></td>

              </tr>
              <tr>
                <td><span style="color: #2874a6;">J</span></td>
                <td>{jmag[i]}</td>
                <td></td>
                <td></td>
              </tr>
                <tr>
                <td><span style="color: #2874a6;">H</span></td>
                <td>{hmag[i]}</td>
                <td></td>
                <td></td>

              </tr>
                <tr>
                <td><span style="color: #2874a6;">K</span></td>
                <td>{kmag[i]}</td>
                <td></td>
                <td></td>
              </tr>
            </table>
            <br>
            <a href="https://exoplanetarchive.ipac.caltech.edu/overview/{obj_sci}">Link to NASA Exoplanet Archive</a><br>
            <a href="http://simbad.u-strasbg.fr/simbad/sim-id?Ident={obj_sci}">Link to Simbad</a>
            </body>
            </html>
            """

        with open('{0}_{1}.html'.format(obj_sci, i + 1), 'w') as f:
            f.write(html_text)
        f.close()

        # HTML to PDF
        HTML('{0}_{1}.html'.format(obj_sci, i + 1)).write_pdf('{0}_{1}.pdf'.format(obj_sci, i + 1))
        os.system('rm {0}_{1}.html'.format(obj_sci, i + 1))
    if info['pl_tranflag'][i] == 1:
        return name, t0, period
    else:
        return name, tp, period


def pdfplots(obj_sci,
             obj_template = None,
             exoarchive_name = None,
             doplot_debug = False,
             nMD_cut = -5):
    """
    Line-by-Line Radial Velocity: Automatic Diagnostics Plots
    obj_sci: Name of the object.
    obj_template: Name of the template.
    exoarchive_name: Name of the object on the Exoplanet Archive
    doplot_debug: If True, will show the plots one by one.
    """

    if obj_template is None:
        obj_template = obj_sci

    if exoarchive_name is None:
        exoarchive_name = obj_sci

    if obj_sci[:3] == 'TOI':
        plname, t_phase, per = system_info(obj_sci)
        nbplanet = len(plname)
    else:
        plname, t_phase, per = system_info(obj_sci, exoarchive_name)
        nbplanet = len(plname)

    tbl = Table.read('lbl_{0}_{1}.rdb'.format(obj_sci, obj_template),
                     format = 'rdb')
    tbl_bin = Table.read('lbl2_{0}_{1}.rdb'.format(obj_sci, obj_template),
                     format = 'rdb')

    # General format
    xspace = 0.05 * (tbl_bin['rjd'].max() - tbl_bin['rjd'].min())
    xrange = [int(tbl_bin['rjd'].min() - xspace),
              int(tbl_bin['rjd'].max() + xspace)]

    # Figures
    pdf = PdfPages('plots.pdf')


    # Title page
    fig = plt.figure(figsize = (10,8))
    plt.axis('off')
    plt.text(0.5, 0.6, 'SPIRou Line-by-Line Radial Velocity',
             ha = 'center', va = 'center', fontsize = 16, fontweight = 'bold')
    plt.text(0.5, 0.5, 'Object: {0}, Template: {1}'.format(obj_sci, obj_template),
             ha = 'center', va = 'center', fontsize = 14)
    plt.text(0.5, 0.4, 'Diagnostic Plots',
             ha = 'center', va = 'center', fontsize = 14)
    plt.text(0.5, 0.1, today.strftime("%Y-%m-%d"),
             ha = 'center', va = 'center', fontsize = 12)
    pdf.savefig(fig)
    if doplot_debug:
        plt.show()
    plt.close(fig)


    # SNR
    # Plot in red data points with SNR Normalized Median Deviation < nMD_cut
    nMD = (tbl['EXTSN035'] - np.nanmedian(tbl['EXTSN035']))/np.nanmedian(
           np.abs(tbl['EXTSN035'] - np.nanmedian(tbl['EXTSN035'])))
    low_SNR = nMD < nMD_cut
    nMD_minus5_value = nMD_cut * np.nanmedian(np.abs(tbl['EXTSN035'] -
            np.nanmedian(tbl['EXTSN035']))) + np.nanmedian(tbl['EXTSN035'])
    fig, ax = plt.subplots(2, 1, figsize = (10,8),
                           gridspec_kw={'height_ratios': [2.5, 1]})
    ax[0].set_title('Signal-to-noise ratio', fontsize = 14)
    ax[0].scatter(tbl['rjd'][~low_SNR], tbl['EXTSN035'][~low_SNR],
                  10, color = 'k', zorder = 3)
    if sum(low_SNR) > 0:
        ax[0].scatter(tbl['rjd'][low_SNR], tbl['EXTSN035'][low_SNR],
                      10, color = 'r', zorder = 3)
    ax[0].axhline(y = np.median(tbl['SNRGOAL']), color = 'k', linestyle = '--',
                  zorder = 0, label = 'SNR Goal')
    ax[0].axhline(y = nMD_minus5_value, color = 'r', linestyle = '--',
                  zorder = 0, label = 'nMD = {0}'.format(nMD_cut))
    ax[0].fill_between(np.linspace(xrange[0], xrange[1], 2), 0, 10, color = 'k',
                       alpha = 0.2, label = 'SNR Threshold = 10')
    ax[0].set_ylabel("SNR for order 35 (around 1.6 $\mu$m)", fontsize = 14)
    ax[0].set_xlabel("Time (BJD - 2400000)", fontsize = 14)
    ax[0].set_xlim([xrange[0], xrange[1]])
    ymax = np.amax([int(1.05 * np.amax(tbl['EXTSN035'])),
                    int(1.05 * np.median(tbl['SNRGOAL']))])
    ax[0].set_ylim([0, ymax])
    ax[0].legend(loc="best", edgecolor = 'k')
    ax[0].minorticks_on()
    ax[0].tick_params(axis="both", labelsize=12, direction = 'in', length = 5)
    ax[0].tick_params(axis="both", which = 'minor', direction = 'in', length = 3)
    ax[0].yaxis.set_ticks_position('both')
    ax[0].xaxis.set_ticks_position('both')
    ax[1].set_title('Histogram', fontsize = 14)
    ax[1].hist(tbl['EXTSN035'], bins = np.histogram_bin_edges(tbl['EXTSN035'],
               'fd'), color = 'k', alpha = 0.8, edgecolor='white',
               linewidth = 1, zorder = 3)
    ax[1].axvline(x = np.median(tbl['SNRGOAL']), color = 'k', linestyle = '--',
                  zorder = 0, label = 'SNR Goal')
    if sum(low_SNR) > 0:
        ax[1].axvline(x = nMD_minus5_value, color = 'r', linestyle = '--',
                      zorder = 0, label = 'nMD = {0}'.format(nMD_cut))
    ax[1].set_ylabel("Nb. of obs. per bin", fontsize = 14)
    ax[1].set_xlabel("SNR for order 35 (around 1.6 $\mu$m)", fontsize = 14)
    ax[1].legend(loc="best", edgecolor = 'k')
    ax[1].minorticks_on()
    ax[1].tick_params(axis="both", labelsize=12, direction = 'in', length = 5)
    ax[1].tick_params(axis="both", which = 'minor', direction = 'in', length = 3)
    ax[1].yaxis.set_ticks_position('both')
    ax[1].xaxis.set_ticks_position('both')
    plt.tight_layout()
    pdf.savefig(fig)
    if doplot_debug:
        plt.show()
    plt.close(fig)


    # RV curves (binned, unbinned)
    fig, ax = plt.subplots(2, 1, figsize = (10,8))
    ax[0].set_title('RV Curve', fontsize = 14)
    ax[0].errorbar(tbl['rjd'][~low_SNR], tbl['vrad'][~low_SNR],
                   yerr = tbl['svrad'][~low_SNR], color = 'k',
                   linestyle = "none", marker = "o", capsize = 2, zorder = 1,
                   label = 'Median error = {0:.2f} m/s'.format(np.median(
                                                        tbl['svrad'])))
    if sum(low_SNR) > 0:
        ax[0].errorbar(tbl['rjd'][low_SNR], tbl['vrad'][low_SNR],
                       yerr = tbl['svrad'][low_SNR], color = 'r',
                       alpha = 0.25, linestyle = "none", marker = "o",
                       capsize = 2, zorder = 1)
    ax[0].set_ylabel("RV (m/s)", fontsize = 14)
    ax[0].set_xlim([xrange[0], xrange[1]])
    ax[0].legend(loc = 'best', edgecolor = 'k')
    ax[0].minorticks_on()
    ax[0].tick_params(axis="both", labelsize=12, direction = 'in',
                      length = 5)
    ax[0].tick_params(axis="both", which = 'minor', direction = 'in',
                      length = 3)
    ax[0].yaxis.set_ticks_position('both')
    ax[0].xaxis.set_ticks_position('both')
    ax[1].set_title('RV Curve (binned per night)', fontsize = 14)
    ax[1].errorbar(tbl_bin['rjd'], tbl_bin['vrad'], yerr = tbl_bin['svrad'],
                   color = 'k', linestyle = "none", marker = "o", capsize = 2,
                   zorder = 1, label = 'Median per night error  = {0:.2f} m/s'.format(
                   np.median(tbl_bin['svrad'])))
    ax[1].set_ylabel("RV (m/s)", fontsize = 14)
    ax[1].set_xlabel("Time (BJD - 2400000)", fontsize = 14)
    ax[1].set_xlim([xrange[0], xrange[1]])
    ax[1].legend(loc = 'best', edgecolor = 'k')
    ax[1].minorticks_on()
    ax[1].tick_params(axis="both", labelsize=12, direction = 'in',
                      length = 5)
    ax[1].tick_params(axis="both", which = 'minor', direction = 'in',
                      length = 3)
    ax[1].yaxis.set_ticks_position('both')
    ax[1].xaxis.set_ticks_position('both')
    plt.tight_layout()
    pdf.savefig(fig)
    if doplot_debug:
        plt.show()
    plt.close(fig)

    # RV, DDV, DDDV and periodograms
    fig, ax = plt.subplots(3, 2, figsize = (20,12))
    ax[0,0].set_title('RV Curve', fontsize = 14)
    ax[0,0].errorbar(tbl['rjd'][~low_SNR], tbl['vrad'][~low_SNR],
                     yerr = tbl['svrad'][~low_SNR], color = 'k',
                     linestyle = "none", marker = "o", capsize = 2, zorder = 1)
    ax[0,0].set_ylabel("RV (m/s)", fontsize = 14)
    ax[0,0].set_xlim([xrange[0], xrange[1]])
    ax[0,0].minorticks_on()
    ax[0,0].tick_params(axis="both", labelsize=12, direction = 'in',
                        length = 5)
    ax[0,0].tick_params(axis="both", which = 'minor', direction = 'in',
                        length = 3)
    ax[0,0].yaxis.set_ticks_position('both')
    ax[0,0].xaxis.set_ticks_position('both')
    ax[1,0].set_title('Line-by-Line 2nd Derivative', fontsize = 14)
    ax[1,0].errorbar(tbl['rjd'][~low_SNR], tbl['per_epoch_DDV'][~low_SNR],
                     tbl['per_epoch_DDVRMS'][~low_SNR], color = 'k',
                     linestyle = "none", marker = "o", capsize = 2, zorder = 1)
    ax[1,0].set_ylabel("DDV", fontsize = 14)
    ax[1,0].set_xlim([xrange[0], xrange[1]])
    ax[1,0].minorticks_on()
    ax[1,0].tick_params(axis="both", labelsize=12, direction = 'in',
                        length = 5)
    ax[1,0].tick_params(axis="both", which = 'minor', direction = 'in',
                        length = 3)
    ax[1,0].yaxis.set_ticks_position('both')
    ax[1,0].xaxis.set_ticks_position('both')
    ax[2,0].set_title('Line-by-Line 3rd Derivative', fontsize = 14)
    ax[2,0].errorbar(tbl['rjd'][~low_SNR], tbl['per_epoch_DDDV'][~low_SNR],
                     tbl['per_epoch_DDDVRMS'][~low_SNR], color = 'k',
                     linestyle = "none", marker = "o", capsize = 2, zorder = 1)
    ax[2,0].set_ylabel("DDDV", fontsize = 14)
    ax[2,0].set_xlabel("Time (BJD - 2400000)", fontsize = 14)
    ax[2,0].set_xlim([xrange[0], xrange[1]])
    ax[2,0].minorticks_on()
    ax[2,0].tick_params(axis="both", labelsize=12, direction = 'in',
                        length = 5)
    ax[2,0].tick_params(axis="both", which = 'minor', direction = 'in',
                        length = 3)
    ax[2,0].yaxis.set_ticks_position('both')
    ax[2,0].xaxis.set_ticks_position('both')

    twilight = plt.get_cmap('twilight_shifted')
    colors = iter(twilight(np.linspace(0,1,10)))
    ls = LombScargle(tbl['rjd'][~low_SNR], tbl['vrad'][~low_SNR],
                     tbl['svrad'][~low_SNR])
    frequency, power = ls.autopower()
    fap = ls.false_alarm_level(0.001)
    ax[0,1].set_title('Lomb-Scargle Periodogram', fontsize = 14)
    ax[0,1].plot(1 / frequency, power, color="k",
                 label = "Lomb-Scargle periodogram")
    ax[0,1].axhline(y=fap,linestyle="--",color="k",label="99.9% confidence")
    wf = LombScargle(tbl['rjd'][~low_SNR], np.ones(len(tbl['rjd'][~low_SNR])),
                     fit_mean = False, center_data = False)
    frequency, power = wf.autopower()
    ax[0,1].plot(1/frequency, power, color="r", label = "Window function")
    for i in range(nbplanet):
        ax[0,1].axvline(x = per[i], linestyle = '--', color = next(colors),
                        label = '{0} ($P$ = {1:.2f} days)'.format(plname[i],
                                                                  per[i]))
    ax[0,1].set_ylabel("Power", fontsize=14)
    ax[0,1].legend(loc="upper left", edgecolor = 'k')
    ax[0,1].set_xlim([0.5,500])
    ax[0,1].set_ylim([0,1])
    ax[0,1].set_xscale('log')
    ax[0,1].tick_params(axis="both", labelsize=12, direction = 'out',
                        length = 5)
    ax[0,1].tick_params(axis="both", which = 'minor', direction = 'out',
                        length = 3)
    twilight = plt.get_cmap('twilight_shifted')
    colors = iter(twilight(np.linspace(0,1,10)))
    ls = LombScargle(tbl['rjd'][~low_SNR], tbl['per_epoch_DDV'][~low_SNR],
                     tbl['per_epoch_DDVRMS'][~low_SNR])
    frequency, power = ls.autopower()
    fap = ls.false_alarm_level(0.001)
    ax[1,1].set_title('Lomb-Scargle Periodogram', fontsize = 14)
    ax[1,1].plot(1 / frequency, power, color="k",
                 label = "Lomb-Scargle periodogram")
    ax[1,1].axhline(y=fap,linestyle="--",color="k",label="99.9% confidence")
    wf = LombScargle(tbl['rjd'][~low_SNR], np.ones(len(tbl['rjd'][~low_SNR])),
                     fit_mean = False, center_data = False)
    frequency, power = wf.autopower()
    ax[1,1].plot(1/frequency, power, color="r", label = "Window function")
    for i in range(nbplanet):
        ax[1,1].axvline(x = per[i], linestyle = '--', color = next(colors),
                        label = '{0} ($P$ = {1:.2f} days)'.format(plname[i],
                                                                  per[i]))
    ax[1,1].set_ylabel("Power", fontsize=14)
    ax[1,1].legend(loc="upper left", edgecolor = 'k')
    ax[1,1].set_xlim([0.5,500])
    ax[1,1].set_ylim([0,1])
    ax[1,1].set_xscale('log')
    ax[1,1].tick_params(axis="both", labelsize=12, direction = 'out',
                        length = 5)
    ax[1,1].tick_params(axis="both", which = 'minor', direction = 'out',
                        length = 3)
    twilight = plt.get_cmap('twilight_shifted')
    colors = iter(twilight(np.linspace(0,1,10)))
    ls = LombScargle(tbl['rjd'][~low_SNR], tbl['per_epoch_DDDV'][~low_SNR],
                     tbl['per_epoch_DDDVRMS'][~low_SNR])
    frequency, power = ls.autopower()
    fap = ls.false_alarm_level(0.001)
    ax[2,1].set_title('Lomb-Scargle Periodogram', fontsize = 14)
    ax[2,1].plot(1 / frequency, power, color="k",
                 label = "Lomb-Scargle periodogram")
    ax[2,1].axhline(y=fap,linestyle="--",color="k",label="99.9% confidence")
    wf = LombScargle(tbl['rjd'][~low_SNR], np.ones(len(tbl['rjd'][~low_SNR])),
                     fit_mean = False, center_data = False)
    frequency, power = wf.autopower()
    ax[2,1].plot(1/frequency, power, color="r", label = "Window function")
    for i in range(nbplanet):
        ax[2,1].axvline(x = per[i], linestyle = '--', color = next(colors),
                        label = '{0} ($P$ = {1:.2f} days)'.format(plname[i],
                                                                  per[i]))
    ax[2,1].set_ylabel("Power", fontsize=14)
    ax[2,1].legend(loc="upper left", edgecolor = 'k')
    ax[2,1].set_xlim([0.5,500])
    ax[2,1].set_ylim([0,1])
    ax[2,1].set_xscale('log')
    ax[2,1].tick_params(axis="both", labelsize=12, direction = 'out',
                        length = 5)
    ax[2,1].tick_params(axis="both", which = 'minor', direction = 'out',
                        length = 3)
    plt.tight_layout()
    pdf.savefig(fig)
    if doplot_debug:
        plt.show()
    plt.close(fig)


    # RV H band, RV Y-H, RV J-H, RV H-K
    fig, ax = plt.subplots(3, 2, figsize = (20,12))
    ax[0,0].set_title('RV Curve (H band only)', fontsize = 14)
    ax[0,0].errorbar(tbl['rjd'][~low_SNR], tbl['vrad_H'][~low_SNR],
                     yerr = tbl['svrad_H'][~low_SNR], color = 'k',
                     linestyle = "none", marker = "o", capsize = 2, zorder = 1,
                     label = 'Median error = {0:.2f} m/s'.format(np.median(
                                                          tbl['svrad_H'])))
    ax[0,0].set_ylabel("RV (m/s)", fontsize = 14)
    ax[0,0].set_xlim([xrange[0], xrange[1]])
    ax[0,0].legend(loc = 'best', edgecolor = 'k')
    ax[0,0].minorticks_on()
    ax[0,0].tick_params(axis="both", labelsize=12, direction = 'in',
                        length = 5)
    ax[0,0].tick_params(axis="both", which = 'minor', direction = 'in',
                        length = 3)
    ax[0,0].yaxis.set_ticks_position('both')
    ax[0,0].xaxis.set_ticks_position('both')
    ax[1,0].set_title('RV Curve (K band only)', fontsize = 14)
    ax[1,0].errorbar(tbl['rjd'][~low_SNR], tbl['vrad_K'][~low_SNR],
                     yerr = tbl['svrad_K'][~low_SNR], color = 'k',
                     linestyle = "none", marker = "o", capsize = 2, zorder = 1,
                     label = 'Median error = {0:.2f} m/s'.format(np.median(
                                                          tbl['svrad_K'])))
    ax[1,0].set_ylabel("Residuals (m/s)", fontsize = 14)
    ax[1,0].set_xlim([xrange[0], xrange[1]])
    ax[1,0].legend(loc = 'best', edgecolor = 'k')
    ax[1,0].minorticks_on()
    ax[1,0].tick_params(axis="both", labelsize=12, direction = 'in',
                        length = 5)
    ax[1,0].tick_params(axis="both", which = 'minor', direction = 'in',
                        length = 3)
    ax[1,0].yaxis.set_ticks_position('both')
    ax[1,0].xaxis.set_ticks_position('both')
    ax[2,0].set_title('RV Curve (H - K)', fontsize = 14)
    ax[2,0].errorbar(tbl['rjd'][~low_SNR], tbl['vrad_H'][~low_SNR] - tbl['vrad_K'][~low_SNR],
                     yerr = np.sqrt(tbl['svrad_H'][~low_SNR]**2 + tbl['svrad_K'][~low_SNR]**2),
                     color = 'k',linestyle = "none", marker = "o", capsize = 2,
                     zorder = 1, label = 'Residuals RMS = {0:.2f} m/s'.format(
                             np.std(tbl['vrad_H'] - tbl['vrad_K'])))
    ax[2,0].axhline(y = 0, linestyle = '--', color = 'k', zorder = 0)
    ax[2,0].set_ylabel("Residuals (m/s)", fontsize = 14)
    ax[2,0].set_xlabel("Time (BJD - 2400000)", fontsize = 14)
    ax[2,0].set_xlim([xrange[0], xrange[1]])
    ax[2,0].legend(loc = 'best', edgecolor = 'k')
    ax[2,0].minorticks_on()
    ax[2,0].tick_params(axis="both", labelsize=12, direction = 'in',
                        length = 5)
    ax[2,0].tick_params(axis="both", which = 'minor', direction = 'in',
                        length = 3)
    ax[2,0].yaxis.set_ticks_position('both')
    ax[2,0].xaxis.set_ticks_position('both')

    twilight = plt.get_cmap('twilight_shifted')
    colors = iter(twilight(np.linspace(0,1,10)))
    ls = LombScargle(tbl['rjd'][~low_SNR], tbl['vrad_H'][~low_SNR],
                     tbl['svrad_H'][~low_SNR])
    frequency, power = ls.autopower()
    fap = ls.false_alarm_level(0.001)
    ax[0,1].set_title('Lomb-Scargle Periodogram', fontsize = 14)
    ax[0,1].plot(1 / frequency, power, color="k",
                 label = "Lomb-Scargle periodogram")
    ax[0,1].axhline(y=fap,linestyle="--",color="k",label="99.9% confidence")
    wf = LombScargle(tbl['rjd'][~low_SNR], np.ones(len(tbl['rjd'][~low_SNR])),
                     fit_mean = False, center_data = False)
    frequency, power = wf.autopower()
    ax[0,1].plot(1/frequency, power, color="r", label = "Window function")
    for i in range(nbplanet):
        ax[0,1].axvline(x = per[i], linestyle = '--', color = next(colors),
                        label = '{0} ($P$ = {1:.2f} days)'.format(plname[i],
                                                                  per[i]))
    ax[0,1].set_ylabel("Power", fontsize=14)
    ax[0,1].legend(loc="upper left", edgecolor = 'k')
    ax[0,1].set_xlim([0.5,500])
    ax[0,1].set_ylim([0,1])
    ax[0,1].set_xscale('log')
    ax[0,1].tick_params(axis="both", labelsize=12, direction = 'out',
                        length = 5)
    ax[0,1].tick_params(axis="both", which = 'minor', direction = 'out',
                        length = 3)
    twilight = plt.get_cmap('twilight_shifted')
    colors = iter(twilight(np.linspace(0,1,10)))
    ls = LombScargle(tbl['rjd'][~low_SNR], tbl['vrad_K'][~low_SNR],
                     tbl['svrad_K'][~low_SNR])
    frequency, power = ls.autopower()
    fap = ls.false_alarm_level(0.001)
    ax[1,1].set_title('Lomb-Scargle Periodogram', fontsize = 14)
    ax[1,1].plot(1 / frequency, power, color="k",
                 label = "Lomb-Scargle periodogram")
    ax[1,1].axhline(y=fap,linestyle="--",color="k",label="99.9% confidence")
    wf = LombScargle(tbl['rjd'][~low_SNR], np.ones(len(tbl['rjd'][~low_SNR])),
                     fit_mean = False, center_data = False)
    frequency, power = wf.autopower()
    ax[1,1].plot(1/frequency, power, color="r", label = "Window function")
    for i in range(nbplanet):
        ax[1,1].axvline(x = per[i], linestyle = '--', color = next(colors),
                        label = '{0} ($P$ = {1:.2f} days)'.format(plname[i],
                                                                  per[i]))
    ax[1,1].set_ylabel("Power", fontsize=14)
    ax[1,1].legend(loc="upper left", edgecolor = 'k')
    ax[1,1].set_xlim([0.5,500])
    ax[1,1].set_ylim([0,1])
    ax[1,1].set_xscale('log')
    ax[1,1].tick_params(axis="both", labelsize=12, direction = 'out',
                        length = 5)
    ax[1,1].tick_params(axis="both", which = 'minor', direction = 'out',
                        length = 3)
    twilight = plt.get_cmap('twilight_shifted')
    colors = iter(twilight(np.linspace(0,1,10)))
    ls = LombScargle(tbl['rjd'][~low_SNR], tbl['vrad_H'][~low_SNR] - tbl['vrad_K'][~low_SNR],
                     np.sqrt(tbl['svrad_H'][~low_SNR]**2 + tbl['svrad_K'][~low_SNR]**2))
    frequency, power = ls.autopower()
    fap = ls.false_alarm_level(0.001)
    ax[2,1].set_title('Lomb-Scargle Periodogram', fontsize = 14)
    ax[2,1].plot(1 / frequency, power, color="k",
                 label = "Lomb-Scargle periodogram")
    ax[2,1].axhline(y=fap,linestyle="--",color="k",label="99.9% confidence")
    wf = LombScargle(tbl['rjd'][~low_SNR], np.ones(len(tbl['rjd'][~low_SNR])),
                     fit_mean = False, center_data = False)
    frequency, power = wf.autopower()
    ax[2,1].plot(1/frequency, power, color="r", label = "Window function")
    for i in range(nbplanet):
        ax[2,1].axvline(x = per[i], linestyle = '--', color = next(colors),
                        label = '{0} ($P$ = {1:.2f} days)'.format(plname[i],
                                                                  per[i]))
    ax[2,1].set_ylabel("Power", fontsize=14)
    ax[2,1].legend(loc="upper left", edgecolor = 'k')
    ax[2,1].set_xlim([0.5,500])
    ax[2,1].set_ylim([0,1])
    ax[2,1].set_xscale('log')
    ax[2,1].tick_params(axis="both", labelsize=12, direction = 'out',
                        length = 5)
    ax[2,1].tick_params(axis="both", which = 'minor', direction = 'out',
                        length = 3)
    plt.tight_layout()
    pdf.savefig(fig)
    if doplot_debug:
        plt.show()
    plt.close(fig)

    # PDF close
    pdf.close()

    # System description page(s)
    merger = PdfFileMerger()
    merger.append('plots.pdf')
    for i in range(nbplanet):
        merger.merge(position = i + 1, fileobj = '{0}_{1}.pdf'.format(obj_sci,
                                                                      i + 1))
    output = open('lbl_{0}_{1}_plots.pdf'.format(obj_sci, obj_template), 'wb')
    merger.write(output)
    output.close()
    os.system('rm plots.pdf')
    for i in range(nbplanet):
        os.system('rm {0}_{1}.pdf'.format(obj_sci, i + 1))
