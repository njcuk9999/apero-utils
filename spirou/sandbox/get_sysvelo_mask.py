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
temperature_model = 3200

mask_file = 'GL699_neg.mas'

# We assume that the mask file has 3 columns, start/end of mask domain and weight. The line center is assume
# to be the mean between the first two columns
tbl = np.loadtxt(mask_file)
wavelines =  (tbl[:,0]+tbl[:,1])/2
weight =  tbl[:,2]

tbl = Table()
tbl['ll_mask_s'] = wavelines
tbl['ll_mask_e'] = wavelines

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



# round temperature to nearest 100 and get the right model
temperature = str(np.int(np.round(temperature_model, -2)))


print('Temperature = ', temperature)
outname = '{0}/lte0{1}-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(path_to_models, temperature)
print('Model file = ', outname)
flux_phoenix = fits.getdata(outname)


# create a spline of the model
model = InterpolatedUnivariateSpline(wave_phoenix, flux_phoenix)

# assume a 0 velocity and search
dv0 = 0
scale = 1.0

systemic_velocity = 0

tbl0 = Table(tbl)

for ite in range(3):
    corrv = np.sqrt((1 + systemic_velocity / c) / (1 - systemic_velocity / c))
    tbl['ll_mask_s'] = tbl0['ll_mask_s'] / corrv
    tbl['ll_mask_e'] = tbl0['ll_mask_e'] / corrv

    wavelines = (tbl['ll_mask_s']+tbl['ll_mask_e'])/2.0

    dvs = np.arange(400, dtype=float)
    dvs -= np.mean(dvs)

    dvs *= scale
    #dvs += systemic_velocity

    neg_mask = weight > 0

    weight_tmp = weight[neg_mask]
    wave_tmp = wavelines[neg_mask]

    cc = np.zeros_like(dvs)
    for i in range(len(dvs)):
        corrv = np.sqrt((1 + dvs[i] / c) / (1 - dvs[i] / c))
        cc[i] = np.sum(weight_tmp*model(wave_tmp / corrv))

    # just centering the cc around one and removing low-f trends
    #cc = (cc / medfilt(cc, 21))

    minpos = np.argmin(cc)
    fit = np.polyfit(dvs[minpos - 1:minpos + 2], cc[minpos - 1:minpos + 2], 2)

    plt.plot(dvs+systemic_velocity, cc,alpha = 0.5)


    systemic_velocity += (-.5 * fit[1] / fit[0])
    print(systemic_velocity)
    scale /= 5.0

plt.title('CCF of model SP with target''s line list\nThis gets you the systemic velocity')
plt.xlabel('Velocity')
plt.ylabel('Abritrary flux')
plt.show()



print('\n\tsystemic velocity of mask : {0:.2f}km/s\n'.format(systemic_velocity))


