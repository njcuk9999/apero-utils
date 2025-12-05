import getpass
import glob
import os
import pickle
import warnings
from datetime import datetime
from astropy.io import fits

import numpy as np
import yaml
from astropy import constants
from scipy import constants
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from tqdm import tqdm


import fitsio

def write_t(data, file):
    # write the dictionary to a MEF file
    with fitsio.FITS(file, 'rw', clobber=True) as outfile:
        for key in data.keys():
            if '_header' not in key:
                data_to_write = data[key]
                header = data[key+'_header']
                header['EXTNAME'] = key  # set the extension name in the header
                outfile.write(data_to_write, header=header,extname=key)
def read_t(file):
    # read the file that is a MEF and makes a dictionary
    # of all extensions
    with fitsio.FITS(file) as infile:
        data = dict()
        for hdu in infile:
            key = hdu.get_extname()
            data[key] = hdu.read()
            data[key + '_header'] = hdu.read_header()
    return data


def load_yaml(yaml_file):
    with open(yaml_file, 'r') as file:
        params = yaml.safe_load(file)
    return params


from astropy.io import fits


def dict2mef(outname, dict_mef):
    """

    :param outname:
    :param dict:

    :return:
    """

    # we check that the folder exists and create it if not
    # we loop from the base folder up
    tmp = outname.split('/')

    if len(tmp) > 3:
        for i in range(len(tmp) - 2, 1, -1):
            folder = '/'.join(tmp[:-i])
            if os.path.exists(folder):
                continue
            else:
                printc('Creating folder {}'.format(folder), 'info')
                os.mkdir(folder)

    h = fits.Header()
    new_hdul = fits.HDUList()

    keys_mef = [np.ndarray]
    for key in dict_mef.keys():
        if type(dict_mef[key]) not in keys_mef:
            with warnings.catch_warnings() as _:
                h[key] = dict_mef[key]
        else:
            new_hdul.append(fits.ImageHDU(data=dict_mef[key], header=h, name=key))

    new_hdul.writeto(outname, overwrite=True)


def mef2dict(filename):
    printc('Reading {}'.format(filename), 'info')
    dict_mef = {}
    with fits.open(filename) as hdu_list:
        h = hdu_list[0].header

        for i in range(1, len(hdu_list)):
            dict_mef[hdu_list[i].name] = hdu_list[i].data

        for key in h.keys():
            dict_mef[key] = h[key]

    remove_keys = ['SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'EXTEND',
                   'COMMENT', 'HISTORY', 'EXTNAME']
    for key in remove_keys:
        if key in dict_mef:
            del dict_mef[key]

    return dict_mef


def color(message, color_txt):
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

    return COLOURS[color_txt.upper()] + message + COLOURS['ENDC']


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


def doppler(wave, v):
    # velocity expressed in m/s
    # relativistic calculation

    v = np.array(v)
    wave = np.array(wave)
    return wave * np.sqrt((1 - v / constants.c) / (1 + v / constants.c))


def lowpassfilter(input_vect, width=501):
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
        xmed.append(np.nanmean(pixval))
        ymed.append(np.nanmedian(input_vect[pixval]))

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


def sigma(im_tmp):
    im = np.array(im_tmp).astype('float64')
    p1 = (1 - 0.682689492137086) / 2
    p2 = 1 - p1
    return (np.nanpercentile(im, p2 * 100) - np.nanpercentile(im, p1 * 100)) / 2


def save_pickle(filename, variable):
    with open(filename, 'wb') as handle:
        pickle.dump(variable, handle)


def read_pickle(filename):
    with open(filename, 'rb') as handle:
        b = pickle.load(handle)
    return b


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


def snail(iter, desc='', leave=False):
    txt = [' ']
    for i in range(1000):
        txt.append('🐌')
    txt.append('_')
    # ,'🐌','🐌','🐌','🐌','🐌','🐌','🐌','🐌','_']
    return tqdm(iter, leave=leave, desc=desc,
                colour='green', ascii=txt)
