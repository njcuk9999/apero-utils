import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
from scipy.signal import medfilt
import os as os

#
# Code to generate a mask that can be used with the DRS from a MODEL spectrum.
# analoguous to the mk_ccf_mask
#


# Path where models are saved
path_to_models = 'HiResFITS'

# some parameters, don't worry
dv = 0.00  # km/s -- width of the CCF box
c = 2.99792458e5 # speed of light
logg = 4.50

# create directory if needed
if not os.path.isdir(path_to_models):
    os.system('mkdir {0}'.format(path_to_models))

# read wavelength and flux. The wavelength is expressed in Ang, we convert to µm
ftp_link = 'ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/'
wave_file = 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
if not os.path.isfile(path_to_models+'/'+wave_file):
    os.system('wget {0}{1}'.format(ftp_link,wave_file) )
    os.system('mv {0} {1}'.format(wave_file,path_to_models))
w = fits.getdata(path_to_models+'/'+wave_file) / 10

spirou_domain = (w>965)*(w<2500)
w = w[spirou_domain]

# get goettigen models if you don't have them.
for temperature in np.arange(3000, 6100, 100):
    outname = '{0}/lte0{1}-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(path_to_models,temperature)

    print('We are downloading/checking presence of model file {0}'.format(outname))

    outname_mask = '{0}/teff{1:.0f}logg{2:.2f}'.format(path_to_models,temperature, logg)

    if not os.path.isfile(outname):
        os.system(
            'wget {0}PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte0{1}-{2:.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(ftp_link,temperature,logg))

        os.system('mv lte0{1:.0f}-{2:.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits {0}'.format(path_to_models,temperature,logg))
    else:
        print('File {0} exists, we are happy!'.format(outname))

    print('\tReading data and finding derivatives')
    f = fits.getdata(outname)
    f = f[spirou_domain]
    f*=w # express in photons, not in energy as models provide things normally
    f/=np.nanmedian(f)
    # find the first and second derivative of the flux. Derivative as a function of log(lambda) is
    # derivative in the velocity space
    df = np.gradient(f)/np.gradient(np.log(w ))
    ddf = np.gradient(df)

    # lines are regions there is a sign change in the derivative of the flux
    # we also have some checks for NaNs
    line = np.where((np.sign(df[1:]) != np.sign(df[:-1])) &
                    np.isfinite(ddf[1:])
                    & np.isfinite(df[1:])
                    & np.isfinite(df[:-1]))[0]

    print('\tConstructing table')
    # create the output table
    tbl = Table()
    tbl['ll_mask_s'] = np.zeros_like(line, dtype=float)
    tbl['ll_mask_e'] = np.zeros_like(line, dtype=float)
    # the weight is the second derivative of the flux. The sharper the line,
    # the more weight we give it
    tbl['w_mask'] = ddf[line]
    tbl['w_mask'] /=np.nanmedian(np.abs(tbl['w_mask']))
    tbl['value'] = f[line]

    tbl['depth'] = np.zeros_like(tbl['value'])

    tbl['depth'][1:-1] = 1 - tbl['value'][1:-1] / ((tbl['value'][0:-2] + tbl['value'][2:]) / 2)

    print('\tcomputing line positions')
    for i in range(len(line)):
        # we perform a linear interpolation to find the exact wavelength
        # where the derivatives goes to zero
        wave_cen = (np.polyfit(df[line[i]:line[i] + 2], w[line[i]:line[i] + 2], 1))[1]

        # we  subtract half of the line width
        corrv = np.sqrt((1 + (-dv / 2) / c) / (1 - (-dv / 2) / c))
        tbl['ll_mask_s'][i] = wave_cen

        # same but for the upper bound to the line position
        corrv = np.sqrt((1 + (dv / 2) / c) / (1 - (dv / 2) / c))
        tbl['ll_mask_e'][i] = wave_cen

    tbl = tbl[np.isfinite(tbl['w_mask'])]
    # tbl = tbl[tbl['w_mask'] > 0]

    wavelines = (tbl['ll_mask_s'] + tbl['ll_mask_e']) / 2.0
    weight = tbl['w_mask']

    pos_mask = tbl['w_mask'] < 0
    neg_mask = tbl['w_mask'] > 0

    print('\twriting {0}'.format(outname_mask + '_pos.csv'))
    tbl[pos_mask].write(outname_mask + '_pos.csv', format='ascii', overwrite=True)

    print('\twriting {0}'.format(outname_mask + '_pos.csv'))
    tbl[neg_mask].write(outname_mask + '_neg.csv', format='ascii', overwrite=True)

    print('\twriting {0}'.format(outname_mask + '_full.mas'))
    tbl.write(outname_mask + '_full.mas', format='ascii', overwrite=True)

    tbl2 = tbl[tbl['w_mask'] > 0]
    tbl2['w_mask'] /= np.nanmedian(tbl2['w_mask'])
    tbl2['depth'] /= np.nanmedian(tbl2['depth'])
    tbl2['depth'] = np.abs(tbl2['depth'])

    print('\twriting {0}'.format(outname_mask + '_neg.mas'))
    f = open(outname_mask + '_neg.mas', 'w')
    for i in range(len(tbl2)):
        f.write('      ' + '      '.join(
            [str(tbl2['ll_mask_s'][i])[0:14], str(tbl2['ll_mask_e'][i])[0:14], str(tbl2['w_mask'][i])[0:12]]) + '\n')
    f.close()

    print('\twriting {0}'.format(outname_mask + '_neg_depth.mas'))
    f = open(outname_mask + '_neg_depth.mas', 'w')
    for i in range(len(tbl2)):
        f.write('      ' + '      '.join(
            [str(tbl2['ll_mask_s'][i])[0:14], str(tbl2['ll_mask_e'][i])[0:14], str(tbl2['depth'][i])[0:12]]) + '\n')
    f.close()

    tbl2 = tbl[tbl['w_mask'] < 0]
    tbl2['w_mask'] /= np.nanmedian(tbl2['w_mask'])

    print('\twriting {0}'.format(outname_mask + '_pos.mas'))
    f = open(outname_mask + '_pos.mas', 'w')
    for i in range(len(tbl2)):
        f.write('      ' + '      '.join(
            [str(tbl2['ll_mask_s'][i])[0:14], str(tbl2['ll_mask_e'][i])[0:14], str(tbl2['w_mask'][i])[0:12]]) + '\n')
    f.close()
