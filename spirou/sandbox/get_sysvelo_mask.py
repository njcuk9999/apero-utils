import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
from scipy.signal import medfilt
import os as os
from scipy import constants

# get the systematic velocity of a mask by comparing it to a Goettigen model


# Model temperature
temperature = 3500

# Path where models are saved
path_to_models = 'HiResFITS'

# some parameters, don't worry
dv = 0.00  # km/s -- width of the CCF box
c = (constants.c/1000)

# create directory if needed
if not os.path.isdir(path_to_models):
    os.system('mkdir {0}'.format(path_to_models))

# read wavelength and flux. The wavelength is expressed in Ang, we convert to µm
ftp_link = 'ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/'
wave_file = 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
if not os.path.isfile(path_to_models+'/'+wave_file):
    os.system('wget {0}{1}'.format(ftp_link,wave_file) )
    os.system('mv {0} {1}'.format(wave_file,path_to_models))
wave_phoenix = fits.getdata(path_to_models+'/'+wave_file) / 10

# get goettigen models if you don't have them.
for temperature in np.arange(3000, 6100, 100):
    temperature = str(np.int(np.round(temperature, -2)))
    outname = '{0}/lte0{1}-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(path_to_models,temperature)

    if not os.path.isfile(outname):
        os.system(
            'wget {0}PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte0{1}-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(ftp_link,temperature))

        os.system('mv lte0{1}-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits {0}'.format(path_to_models,temperature,))
    else:
        print('File {0} exists, we are happy!'.format(outname))


# read template and header
tbl, hdr = fits.getdata(template, ext=1, header=True)

# round temperature to nearest 100 and get the right model
temperature = str(np.int(np.round(temperature, -2)))


print('Temperature = ', temperature)
outname = '{0}/lte0{1}-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(path_to_models, temperature)
print('Model file = ', outname)
flux_phoenix = fits.getdata(outname)


wavelines = (tbl['ll_mask_s'] + tbl['ll_mask_e']) / 2.0
weight = tbl['w_mask']

# create a spline of the model
model = InterpolatedUnivariateSpline(wave_phoenix, flux_phoenix)

# assume a 0 velocity and search
dv0 = 0
scale = 1.0
for ite in range(2):
    dvs = np.arange(400, dtype=float)
    dvs -= np.mean(dvs)

    dvs *= scale
    dvs += dv0

    cc = np.zeros_like(dvs)
    for i in range(len(dvs)):
        corrv = np.sqrt((1 + dvs[i] / c) / (1 - dvs[i] / c))
        cc[i] = np.sum(model(wavelines / corrv))

    # just centering the cc around one and removing low-f trends
    mini = np.argmin(cc / medfilt(cc, 21))
    dv0 = dvs[mini]
    scale /= 10.0

    plt.plot(dvs, cc)
    plt.title('CCF of model SP with target''s line list\nThis gets you the systemic velocity')
    plt.xlabel('Velocity')
    plt.ylabel('Abritrary flux')
    plt.show()

minpos = np.argmin(cc)
fit = np.polyfit(dvs[minpos - 1:minpos + 2], cc[minpos - 1:minpos + 2], 2)

systemic_velocity = -.5 * fit[1] / fit[0]

print('\n\tsystemic velocity : {0:.2f}km/s\n'.format(systemic_velocity))


