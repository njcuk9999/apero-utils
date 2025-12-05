from astropy.table import Table
import numpy as np
from pixtools import snail, sigma, read_t, write_t
import glob
from scipy.interpolate import InterpolatedUnivariateSpline as ius
import os
from scipy.constants import c

path_to_calibs = '/Volumes/courlan/lbl_NIRPS_HE/calib/'
path = '/Volumes/courlan/lbl_NIRPS_HE/science/PROXIMA_TELLU05/*TELLU05_t.fits'
cavity_files = '/Volumes/courlan/lbl_NIRPS_HE/calib/cavity_slope.csv'
replace_segment = 'TELLU05', 'SLINKY'

tbl_cavity = Table.read(cavity_files)

files = glob.glob(path)

for file in snail(files):

    if os.path.isfile(file.replace('TELLU05_t', 'slinky_t')):
        continue

    data_dict = read_t(file)
    h = data_dict['FluxA_header']
    wavefile = h['WAVEFILE']

    g = tbl_cavity['WAVEFILE'] == wavefile
    if np.sum(g) != 1:
        print('no cavity found, skipping')
        continue
    slope = tbl_cavity[g]['slope'][0]
    intercept = tbl_cavity[g]['offset_1600'][0]

    slinky_name = path_to_calibs + h['WAVEFILE'].replace('.fits', '_sliky.csv')
    if not os.path.isfile(slinky_name):
        print('slinky not found, skipping')
        continue

    wave = np.array(data_dict['WaveA'])

    slinky = Table.read(slinky_name)
    corr = np.array(slinky['corr'])
    wave_corr = np.array(slinky['wavelength'])
    valid = np.isfinite(corr + wave_corr)
    corr = corr[valid]
    wave_corr = wave_corr[valid]
    spl = ius(wave_corr, corr, ext=1, k=1)

    off = (wave - 1600) * slope + intercept
    data_dict['WaveA'] *= (1 - off / c)

    write_t(data_dict, file.replace('TELLU05_t', 'slope_t'))

    data_dict['WaveA'] *= (1 - spl(wave) / c)

    write_t(data_dict, file.replace('TELLU05_t', 'slinky_t'))
