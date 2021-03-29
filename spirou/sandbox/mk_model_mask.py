import numpy as np
from astropy.io import fits
from astropy.table import Table
import os as os
from tqdm import tqdm
#
# Code to generate a mask that can be used with the DRS from a MODEL spectrum.
# analoguous to the mk_ccf_mask
#

# Path where models are saved
path_to_models = 'HiResFITS'

# some parameters, don't worry
dv = 0.00  # km/s -- width of the CCF box
c = 2.99792458e5 # speed of light

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

good_domain = (w>350)*(w<2500)
w = w[good_domain]

for logg in ([6.0,5.5,5.0,4.5,4.0]):
    # get goettigen models if you don't have them.
    for temperature in np.arange(2700, 6100, 100):
        outname = '{0}/lte0{1}-{2:.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(path_to_models,temperature,logg)

        print('We are downloading/checking presence of model file {0}'.format(outname))

        # we need to remove the '.' in the outname_mask as it messes the name of the ccf file
        outname_mask = outname.split('.fits')[0]+'.csv'
        if os.path.isfile(outname_mask):
            print('File {0} exists, we are happy!'.format(outname_mask))
            continue

        if not os.path.isfile(outname):
            os.system(
                'wget {0}PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte0{1}-{2:.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(ftp_link,temperature,logg))

            os.system('mv lte0{1:.0f}-{2:.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits {0}'.format(path_to_models,temperature,logg))
        else:
            print('File {0} exists, we are happy!'.format(outname))

        print('\tReading data and finding derivatives')
        f = fits.getdata(outname)
        f = f[good_domain]
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
                        & np.isfinite(df[:-1]) & (ddf[1:]>0))[0]

        weight = ddf[line]

        for ite in range(9):
            print(len(line))
            keep = np.ones_like(weight,dtype = bool)

            for i in range(5,len(line)-6):
                # line is less than 1/5th of best neighbour from -5 to +5 in the list
                if weight[i]/np.max(weight[i-5:i+5]) < 0.1:
                    keep[i] = False
            line = line[keep]
            weight = weight[keep]



        print('\tConstructing table')
        # create the output table
        tbl = Table()
        tbl['wavelength'] = np.zeros_like(line, dtype=float)
        tbl['weight'] = weight
        tbl['weight'] /=np.nanmedian(np.abs(tbl['weight']))
        tbl['weight'] = np.round(tbl['weight'],4)

        print('\tcomputing line positions')
        wave_cen = np.zeros_like(tbl,dtype = float)
        for i in tqdm(range(len(line))):
            # we perform a linear interpolation to find the exact wavelength
            # where the derivatives goes to zero
            wave_cen[i] = (np.polyfit(df[line[i]:line[i] + 2], w[line[i]:line[i] + 2], 1))[1]
        tbl['wavelength'] = np.round(wave_cen,6)

        tbl = tbl[np.isfinite(tbl['wavelength'])]



        print('\twriting {0}'.format(outname_mask))
        tbl.write(outname_mask, overwrite=True)

