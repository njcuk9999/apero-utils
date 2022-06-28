#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-02-10

@author: cook
"""
from astropy.coordinates import EarthLocation
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import time
import warnings

import barycorrpy_berv
import oli_berv
import pyasl_berv
import astropy_berv

# =============================================================================
# Define variables
# =============================================================================
# Note coordinates, pm and plx must be from same time (SIMBAD does not provide
#   this)

# Here we use all Gaia EDR3
INPUT_RA = 269.448502525438 # deg
INPUT_DEC = 4.73942005111241 # deg
INPUT_PMRA = -801.55097836847 # mas / yr
INPUT_PMDE = 10362.3942065465 # mas / yr
INPUT_PLX = 546.975939730948   # mas
INPUT_COORD_TIME = 2457388.5   # JD
# ----------------------------------------------------------------------------
loc = EarthLocation.of_site("Paranal")
#loc = EarthLocation.of_site("Canada-France-Hawaii Telescope")
INPUT_LONG = loc.lon.value     # in deg
INPUT_LAT = loc.lat.value      # in deg
INPUT_ALT = loc.height.value   # in m
# ----------------------------------------------------------------------------
START_DATE = Time('2020-01-01 00:00:00', format='iso')
END_DATE = Time('2021-01-01 00:00:00', format='iso')
STEPS = 100
# ----------------------------------------------------------------------------
# ra / dec in degrees (forms a map)
GRID_RA = np.linspace(0, 359, 10)
GRID_DEC = np.linspace(-89, 89, 10)
# pmra in mas / yr
GRID_PMRA = np.linspace(-6000, 6000, 6)
# pmde in mas / yr
GRID_PMDE = np.linspace(-6000, 6000, 6)
# plx in mas
GRID_PLX = 10**(np.linspace(0, 3, 6))

# select mode (simple or grid)
MODE = 'grid'

# ----------------------------------------------------------------------------
# where things are saved
workspace = '/data/nirps_he/misc/berv_comp'
# output pdf
outfile = 'berv_comparison_{0}.pdf'
# output grid
grid_file = 'grid_{0}.fits'



# =============================================================================
# Define functions
# =============================================================================
def grid_plot(combinations, variables, shapes, labels, sky_shape,
              r_astropy, r_oli, r_pyasl, r_barycorrpy, pdf):
    # loop around each combination group
    for variable in combinations:
        # get combination
        combination = combinations[variable]
        comb_name = labels[variable]
        # calculate the delta barycorrpy results
        dr_astropy = abs(r_astropy[variable] - r_barycorrpy[variable]) * 1000
        dr_oli = abs(r_oli[variable] - r_barycorrpy[variable]) * 1000
        dr_pyasl = abs(r_pyasl[variable] - r_barycorrpy[variable]) * 1000

        # force to the shape of sky
        sky_dr_astropy = np.max(dr_astropy, axis=1).reshape(shapes[variable])
        sky_dr_oli = np.max(dr_oli, axis=1).reshape(shapes[variable])
        sky_dr_pyasl = np.max(dr_pyasl, axis=1).reshape(shapes[variable])

        # deal with number of rows
        if shapes[variable] == sky_shape:
            rows = 1
            figsize = (18, 6)
        else:
            rows = shapes[variable][-1]
            figsize = (18, 18)

        # left, right, bottom, top
        extent = [np.min(GRID_RA), np.max(GRID_RA),
                  np.min(GRID_DEC), np.max(GRID_DEC)]
        # set up the imshow plot kwargs
        plotkwargs = dict(extent=extent, aspect='auto',
                          interpolation='None')
        # set up this plot
        fig, frames = plt.subplots(ncols=3, nrows=rows, sharex='all',
                                   sharey='all', figsize=figsize)
        # loop around rows
        for row in range(rows):
            # print progress
            print('GRID_PLOT: {0} row = {1}'.format(comb_name, row))
            # deal with differing shapes
            if shapes[variable] == sky_shape:
                frame = frames
                title = '{0} - barycorrpy'
                im1 = frame[0].imshow(sky_dr_astropy, **plotkwargs)
                im2 = frame[1].imshow(sky_dr_oli, **plotkwargs)
                im3 = frame[2].imshow(sky_dr_pyasl, **plotkwargs)
            else:
                vlabel = variables[variable][row]
                frame = frames[row]
                # get title
                titlestr = ' ({0}={1})'.format(labels[variable], vlabel)
                title = '{0} - barycorrpy' + titlestr
                im1 = frame[0].imshow(sky_dr_astropy[:, :, row], **plotkwargs)
                im2 = frame[1].imshow(sky_dr_oli[:, :, row], **plotkwargs)
                im3 = frame[2].imshow(sky_dr_pyasl[:, :, row], **plotkwargs)

            # -----------------------------------------------------------------
            # set up color bars - frame 0
            divider = make_axes_locatable(frame[0])
            cax1 = divider.append_axes('right', size='5%', pad=0.05)
            cbar1 = fig.colorbar(im1, cax=cax1, orientation='vertical')
            cbar1.ax.get_yaxis().labelpad = 15
            cbar1.ax.set_ylabel('max dBERV (m/s)', rotation=270)
            # -----------------------------------------------------------------
            # set up color bars - frame 1
            divider = make_axes_locatable(frame[1])
            cax2 = divider.append_axes('right', size='5%', pad=0.05)
            cbar2 = fig.colorbar(im2, cax=cax2, orientation='vertical')
            cbar2.ax.get_yaxis().labelpad = 15
            cbar2.ax.set_ylabel('max dBERV (m/s)', rotation=270)
            # -----------------------------------------------------------------
            # set up color bars - frame 2
            divider = make_axes_locatable(frame[2])
            cax3 = divider.append_axes('right', size='5%', pad=0.05)
            cbar3 = fig.colorbar(im3, cax=cax3, orientation='vertical')
            cbar3.ax.get_yaxis().labelpad = 15
            cbar3.ax.set_ylabel('max dBERV (m/s)', rotation=270)
            # -----------------------------------------------------------------
            # set labels and titles
            if row == rows -1:
                frame[0].set_xlabel('ra (deg)')
                frame[1].set_xlabel('ra (deg)')
                frame[2].set_xlabel('ra (deg)')
            frame[0].set_ylabel('dec (deg)')
            frame[0].set_title(title.format('astropy'), y=1.0, pad=-14,
                                     color='white')
            frame[1].set_title(title.format('oli'), y=1.0, pad=-14,
                                     color='white')
            frame[2].set_title(title.format('pyasl'), y=1.0, pad=-14,
                                     color='white')

        fig.subplots_adjust(left=0.05, right=0.95, bottom=0.075, top=0.925)
        fig.suptitle('Variable with {0}'.format(comb_name))
        pdf.savefig()
        plt.close()



def grid_simple_plot(combinations, times, labels,
                     r_astropy, r_oli, r_pyasl, r_barycorrpy, pdf):
    # add some single plots
    for variable in combinations:
        # get combination
        combination = combinations[variable]
        comb_name = labels[variable]
        # calculate the delta barycorrpy results
        dr_astropy = abs(r_astropy[variable] - r_barycorrpy[variable]) * 1000
        dr_oli = abs(r_oli[variable] - r_barycorrpy[variable]) * 1000
        dr_pyasl = abs(r_pyasl[variable] - r_barycorrpy[variable]) * 1000
        # get minimum and maximum value across all values
        pmax_astropy = np.argmax(np.max(dr_astropy, axis=1))
        pmax_oli = np.argmax(np.max(dr_oli, axis=1))
        pmax_pyasl = np.argmax(np.max(dr_pyasl, axis=1))
        positions = [pmax_astropy, pmax_oli, pmax_pyasl]

        ptitle = ['{0} (max {1} dBERV)'.format(comb_name, 'astropy'),
                  '{0} (max {1} dBERV)'.format(comb_name, 'oli'),
                  '{0} (max {1} dBERV)'.format(comb_name, 'pyasl')]
        pscale = ['astropy', 'oli', 'pyasl']
        # loop around positions
        for p_it, pos in enumerate(positions):
            # print progress
            print('SIMPLE PLOT: {0} position = {1}'.format(comb_name, pos))
            # get parameters
            kwargs = dict(coordtime=INPUT_COORD_TIME,
                          ra=combination[pos][0][0],
                          dec=combination[pos][0][1],
                          pmra=combination[pos][1],
                          pmde=combination[pos][2],
                          plx=combination[pos][3],
                          obs_long=INPUT_LONG, obs_lat=INPUT_LAT,
                          obs_alt=INPUT_ALT)
            # simple plot
            simple_plot(times,
                        r_astropy[variable][pos], r_oli[variable][pos],
                        r_pyasl[variable][pos], r_barycorrpy[variable][pos],
                        pdf=pdf, title=ptitle[p_it], scale=pscale[p_it],
                        **kwargs)


def simple_plot(times, r_astropy, r_oli, r_pyasl, r_barycorrpy,
                ra, dec, pmra, pmde, plx, pdf, **kwargs):
    fig = plt.figure(figsize=(16, 16))
    shape = (9, 1)
    frame1 = plt.subplot2grid(shape, (0, 0), rowspan=2)
    frame2 = plt.subplot2grid(shape, (2, 0), rowspan=2)
    frame3 = plt.subplot2grid(shape, (4, 0), rowspan=2)
    frame4 = plt.subplot2grid(shape, (6, 0), rowspan=2)
    frame5 = plt.subplot2grid(shape, (8, 0), rowspan=1)

    frame1.plot(times, r_astropy, label='astropy', color='blue')
    frame1.plot(times, r_oli, label='oli', color='orange')
    frame1.plot(times, r_pyasl, label='pyasl', color='green')
    frame1.plot(times, r_barycorrpy, label='barycorrpy', color='red')
    frame1.legend(loc=0)
    frame1.set(ylabel='BERV [km/s]')
    frame1.set_xticklabels([])

    dr_astropy = (r_astropy - r_barycorrpy) * 1000
    dr_oli = (r_oli - r_barycorrpy) * 1000
    dr_pyasl = (r_pyasl - r_barycorrpy) * 1000

    frame2.plot(times, dr_astropy, label='astropy - barycorrpy', color='blue')
    frame3.plot(times, dr_oli, label='oli - barycorrpy', color='orange')
    frame4.plot(times, dr_pyasl, label='pyasl - barycorrpy', color='green')

    frame2.legend(loc=0)
    frame2.set(ylabel='dBERV [m/s]')
    frame2.set_xticklabels([])
    frame3.legend(loc=0)
    frame3.set(ylabel='dBERV [m/s]')
    frame3.set_xticklabels([])
    frame4.legend(loc=0)
    frame4.set(xlabel='JD', ylabel='dBERV [m/s]')

    # add table
    columns = ['ra', 'dec', 'pmra', 'pmde', 'plx']
    rows = ['']
    data = [[ra, dec, pmra, pmde, plx]]

    frame5.axis('off')
    tbl = frame5.table(cellText=data, rowLabels=rows, colLabels=columns,
                          bbox = [0.1, 0.3, 0.8, 0.4])
    tbl.scale(1, 1.5)

    if 'title' in kwargs:
        fig.suptitle(kwargs['title'])

    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95,
                        hspace=0.2)

    pdf.savefig()
    plt.close()


def simple_test():
    """
    Use the INPUT_XXXX parameters to generate the difference for one
    specific target

    :return:
    """
    # get a range of dates (in JD)
    times = np.linspace(START_DATE.jd, END_DATE.jd, STEPS)

    # first test single coords
    kwargs = dict(coordtime=INPUT_COORD_TIME,
                  ra=INPUT_RA, dec=INPUT_DEC, pmra=INPUT_PMDE,
                  pmde=INPUT_PMDE, plx=INPUT_PLX, obs_long=INPUT_LONG,
                  obs_lat=INPUT_LAT, obs_alt=INPUT_ALT)
    # -------------------------------------------------------------------------
    # calculate bervs: astropy
    # -------------------------------------------------------------------------
    print('Calculating BERVS for astropy')
    start = time.time()
    r_astropy = astropy_berv.get_bervs(times, **kwargs)
    end = time.time()
    time_astropy = end - start
    # -------------------------------------------------------------------------
    # calculate bervs: oli
    # -------------------------------------------------------------------------
    print('Calculating BERVS for oli')
    start = time.time()
    r_oli = oli_berv.get_bervs(times, **kwargs)
    end = time.time()
    time_oli = end - start
    # -------------------------------------------------------------------------
    # calculate bervs: pyasl
    # -------------------------------------------------------------------------
    print('Calculating BERVS for pyasl')
    start = time.time()
    r_pyasl = pyasl_berv.get_bervs(times, **kwargs)
    end = time.time()
    time_pyasl = end - start
    # -------------------------------------------------------------------------
    # calculate bervs: barycorrpy
    # -------------------------------------------------------------------------
    print('Calculating BERVS for pyasl')
    start = time.time()
    r_barycorrpy = barycorrpy_berv.get_bervs(times, **kwargs)
    end = time.time()
    time_barycorrpy = end - start

    # -------------------------------------------------------------------------
    # Print stats
    # -------------------------------------------------------------------------
    print('Astropy: {0} s'.format(time_astropy))
    print('Oli: {0} s'.format(time_oli))
    print('PyASL: {0} s'.format(time_pyasl))
    print('Barycorrpy: {0} s'.format(time_barycorrpy))

    # -------------------------------------------------------------------------
    # Plot
    # -------------------------------------------------------------------------
    # construct full path for pdf
    outpdf = os.path.join(workspace, outfile.format('simple'))
    # plot pdf
    plt.close()
    with PdfPages(outpdf) as pdf:
        simple_plot(times, r_astropy, r_oli, r_pyasl, r_barycorrpy, pdf=pdf,
                    **kwargs)


def grid_test():
    """
    Use the grid parameters to
    :return:
    """
    # vary sky position
    sky_positions = list(itertools.product(GRID_RA, GRID_DEC))
    sky_shape = [len(GRID_RA), len(GRID_DEC)]

    combinations = dict()
    shapes = dict()
    labels = dict()
    variables = dict()
    # -------------------------------------------------------------------------
    combinations1 = list(itertools.product(sky_positions,
                                           [np.mean(GRID_PMRA)],
                                           [np.mean(GRID_PMDE)],
                                           [np.mean(GRID_PLX)]))
    second_axis = None
    combinations[0] = combinations1
    shapes[0] = sky_shape
    labels[0] = 'sky'
    variables[0] = None
    # -------------------------------------------------------------------------
    # pmra with sky position
    combinations2 = list(itertools.product(sky_positions, GRID_PMRA,
                                           [np.mean(GRID_PMDE)],
                                           [np.mean(GRID_PLX)]))

    combinations[1] = combinations2
    shapes[1] = sky_shape + [len(GRID_PMRA)]
    labels[1] = 'pmra'
    variables[1] = GRID_PMRA
    # -------------------------------------------------------------------------
    # pmde with sky position
    combinations3 = list(itertools.product(sky_positions,
                                           [np.mean(GRID_PMRA)],
                                           GRID_PMDE,
                                           [np.mean(GRID_PLX)]))
    combinations[2] = combinations3
    shapes[2] = sky_shape + [len(GRID_PMDE)]
    labels[2] = 'pmde'
    variables[2] = GRID_PMDE
    # # -------------------------------------------------------------------------
    # plx with sky position
    combinations4 = list(itertools.product(sky_positions,
                                           [np.mean(GRID_PMRA)],
                                           [np.mean(GRID_PMDE)],
                                           GRID_PLX))
    combinations[3] = combinations4
    shapes[3] = sky_shape + [len(GRID_PLX)]
    labels[3] = 'plx'
    variables[3] = GRID_PLX
    # -------------------------------------------------------------------------
    # storage for results
    r_astropy = dict()
    r_oli = dict()
    r_pyasl = dict()
    r_barycorrpy = dict()
    # get a range of dates (in JD)
    times = np.linspace(START_DATE.jd, END_DATE.jd, STEPS)
    # -------------------------------------------------------------------------
    # calculate berv
    # -------------------------------------------------------------------------
    for variable in combinations:
        # get combination list for this variable
        combination = combinations[variable]
        comb_name = labels[variable]
        # save the outputs
        r_astropy[variable] = np.zeros([len(combination), STEPS])
        r_oli[variable] = np.zeros([len(combination), STEPS])
        r_pyasl[variable] = np.zeros([len(combination), STEPS])
        r_barycorrpy[variable] = np.zeros([len(combination), STEPS])
        # just measure how long the iteration took
        run_time = 0
        # loop through combinations
        for c_it in range(len(combination)):
            # -----------------------------------------------------------------
            # print progress (and a timing for iterations > 0)
            msg = '{0} Iteration {1} of {2} ({3:.5f}%)'
            if c_it > 0:
             msg += '[Previous iteration took {4:.3f} s]'
            it = c_it + 1
            args = [comb_name, it, len(combination),
                    100 * (it / len(combination)), run_time]
            print(msg.format(*args))
            # -----------------------------------------------------------------
            # set up kwargs
            kwargs = dict(coordtime=INPUT_COORD_TIME,
                          ra=combination[c_it][0][0],
                          dec=combination[c_it][0][1],
                          pmra=combination[c_it][1],
                          pmde=combination[c_it][2],
                          plx=combination[c_it][3],
                          obs_long=INPUT_LONG, obs_lat=INPUT_LAT,
                          obs_alt=INPUT_ALT)
            # ---------------------------------------------------------------------
            # start timer
            start = time.time()
            # ---------------------------------------------------------------------
            # calculate bervs: astropy
            # ---------------------------------------------------------------------
            with warnings.catch_warnings(record=True) as _:
                berv_astropy = astropy_berv.get_bervs(obstime=times, **kwargs)
                r_astropy[variable][c_it] = berv_astropy
            # ---------------------------------------------------------------------
            # calculate bervs: oli
            # ---------------------------------------------------------------------
            berv_oli = oli_berv.get_bervs(obstime=times, **kwargs)
            r_oli[variable][c_it] = berv_oli
            # ---------------------------------------------------------------------
            # calculate bervs: pyasl
            # ---------------------------------------------------------------------
            berv_pyasl = pyasl_berv.get_bervs(obstime=times, **kwargs)
            r_pyasl[variable][c_it] = berv_pyasl
            # ---------------------------------------------------------------------
            # calculate bervs: barycorrpy
            # ---------------------------------------------------------------------
            berv_barycorrpy = barycorrpy_berv.get_bervs(obstime=times, **kwargs)
            r_barycorrpy[variable][c_it] = berv_barycorrpy
            # ---------------------------------------------------------------------
            # end timer
            end = time.time()
            # update run time
            run_time = end - start

    # ---------------------------------------------------------------------
    # save the data
    # ---------------------------------------------------------------------


    # ---------------------------------------------------------------------
    # sky plots
    # ---------------------------------------------------------------------
    # construct grid name
    gridname = 'grid_{0}'.format(np.product(sky_shape))
    # construct full path for pdf
    outpdf = os.path.join(workspace, outfile.format(gridname))
    # plot pdf
    plt.close()
    with PdfPages(outpdf) as pdf:
        grid_plot(combinations, variables, shapes, labels, sky_shape,
                  r_astropy, r_oli, r_pyasl, r_barycorrpy, pdf)

        grid_simple_plot(combinations, times, labels,
                         r_astropy, r_oli, r_pyasl, r_barycorrpy, pdf)


def save_grid(combinations, variables, shapes, labels, sky_shape,
              r_astropy, r_oli, r_pyasl, r_barycorrpy):

    # for variable in combinations:
    #     hdu0 = fits.PrimaryHDU()
    #     hdu1 = fits.BinTableHDU(save_comb_table(combinations[variable]))

    pass

def save_comb_table(combination) -> Table:

    ncomb = np.array(combination, dtype=object)

    table = Table()
    table['coords'] = ncomb[:, 0]
    table['pmra'] = ncomb[:, 1]
    table['pmde'] = ncomb[:, 2]
    table['plx'] = ncomb[:, 3]

    return table


def load_grid():
    pass


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    # simple test with static coordinates
    if MODE == 'simple':
        simple_test()

    # grid test
    if MODE == 'grid':
        grid_test()


# =============================================================================
# End of code
# =============================================================================
