import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
from scipy.signal import medfilt
import os as os

#
# Code to generate a mask that can be used with the DRS. You need to provide a template file name
# and the code finds features and correlates the mask against a model at the proper temperature
# (from header or 3600K if not provided). This gets you the systemic velocity, so you can offset
# your lines to a zero velocity. The code creates both .csv and .mas files. The csv can be read
# with astropy.table, while the .mas files are in the right format for the DRS. You get masks for
# the negative and positive spectroscopic features (_neg and _pos) masks. You also get a _full mask,
# but the DRS does not handle both positive and negative features yet.
#

# Provide a template file, needs to be a _s1d_v file
template = 'Template_s1d_Gl846_sc1d_v_file_AB.fits'

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

# round temperature in header to nearest 100 and get the right model
if 'OBJTEMP' in hdr:
    temperature = hdr['OBJTEMP']
    if temperature < 3000:
        temperature = 3000
    if temperature > 6000:
        temperature = 6000

    temperature = str(np.int(np.round(temperature, -2)))
else:
    temperature = '3600'

print('Temperature = ', temperature)
outname = '{0}/lte0{1}-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(path_to_models, temperature)
print('Model file = ', outname)
flux_phoenix = fits.getdata(outname)

# get wave and flux vectors
w = np.array(tbl['wavelength'])
f = np.array(tbl['flux'])

f2 = np.array(f)
mask = np.isfinite(f)
f2[~mask] = 0
mask = mask*1.0
f = np.convolve(f2,np.ones(5), mode = 'same')/np.convolve(mask,np.ones(5), mode = 'same')

# find the first and second derivative of the flux
df = np.gradient(f)
ddf = np.gradient(np.gradient(f))

# lines are regions there is a sign change in the derivative of the flux
# we also have some checks for NaNs
line = np.where((np.sign(df[1:]) != np.sign(df[:-1])) &
                np.isfinite(ddf[1:])
                & np.isfinite(df[1:])
                & np.isfinite(df[:-1]))[0]

# create the output table
tbl = Table()
tbl['ll_mask_s'] = np.zeros_like(line, dtype=float)
tbl['ll_mask_e'] = np.zeros_like(line, dtype=float)
# the weight is the second derivative of the flux. The sharper the line,
# the more weight we give it
tbl['w_mask'] = ddf[line]
tbl['value'] = f[line]

tbl['depth'] = np.zeros_like(tbl['value'])

tbl['depth'][1:-1] = 1-tbl['value'][1:-1]/((tbl['value'][0:-2]+tbl['value'][2:])/2)

plt.plot(tbl['w_mask'],tbl['depth'],'g.')
plt.show()

for i in range(len(line)):
    # we perform a linear interpolation to find the exact wavelength
    # where the derivatives goes to zero
    wave_cen = (np.polyfit(df[line[i]:line[i] + 2], w[line[i]:line[i] + 2], 1))[1]

    # we offset that wavelength by the systemic velocity and subtract
    # half of the line width
    corrv = np.sqrt((1 + (-dv / 2) / c) / (1 - (-dv / 2) / c))
    tbl['ll_mask_s'][i] = wave_cen * corrv

    # same but for the upper bound to the line position
    corrv = np.sqrt((1 + (dv / 2) / c) / (1 - (dv / 2) / c))
    tbl['ll_mask_e'][i] = wave_cen * corrv

tbl = tbl[np.isfinite(tbl['w_mask'])]
#tbl = tbl[tbl['w_mask'] > 0]

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

hdr['SYSVELO'] = systemic_velocity, 'meas. systemic velocity (km/s)'
hdr['VELOFILE'] = outname, 'model used for SYSVEL cc'

print('\n\tsystemic velocity : {0:.2f}km/s\n'.format(systemic_velocity))

plt.plot(w,f, 'g-',label = 'input spectrum')
plt.vlines(tbl[tbl['w_mask'] < 0]['ll_mask_s'], np.nanmin(f), np.nanmax(f), 'k',alpha = 0.2,label = 'positive feature')
plt.vlines(tbl[tbl['w_mask'] > 0]['ll_mask_s'], np.nanmin(f), np.nanmax(f), 'r',alpha = 0.2,label = 'negative feature')
plt.legend()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Arbitrary flux')
plt.show()


corrv = np.sqrt((1 + systemic_velocity / c) / (1 - systemic_velocity / c))

# updating the table to account for systemic velocity of star
tbl['ll_mask_s'] = tbl['ll_mask_s'] / corrv
tbl['ll_mask_e'] = tbl['ll_mask_e'] / corrv

# write the output table
fits.writeto(hdr['OBJECT'] + '.fits', tbl, hdr, overwrite=True)

pos_mask = tbl['w_mask']<0
neg_mask = tbl['w_mask']>0

tbl['w_mask']/=np.nanmean(np.abs(tbl['w_mask']))

tbl[pos_mask].write(hdr['OBJECT'] + '_pos.csv', format='ascii', overwrite=True)
tbl[neg_mask].write(hdr['OBJECT'] + '_neg.csv', format='ascii', overwrite=True)
tbl.write(hdr['OBJECT'] + '_full.mas', format='ascii', overwrite=True)


tbl2 = tbl[tbl['w_mask'] > 0]
tbl2['w_mask'] /= np.nanmedian(tbl2['w_mask'])
tbl2['depth'] /= np.nanmedian(tbl2['depth'])
tbl2['depth'] = np.abs(tbl2['depth'])

f = open(hdr['OBJECT'] + '_neg.mas', 'w')
for i in range(len(tbl2)):
    f.write('      ' + '      '.join(
        [str(tbl2['ll_mask_s'][i])[0:14], str(tbl2['ll_mask_e'][i])[0:14], str(tbl2['w_mask'][i])[0:12]]) + '\n')
f.close()

f = open(hdr['OBJECT'] + '_neg_depth.mas', 'w')
for i in range(len(tbl2)):
    f.write('      ' + '      '.join(
        [str(tbl2['ll_mask_s'][i])[0:14], str(tbl2['ll_mask_e'][i])[0:14], str(tbl2['depth'][i])[0:12]]) + '\n')
f.close()



tbl2 = tbl[tbl['w_mask'] < 0]
tbl2['w_mask'] /= np.nanmedian(tbl2['w_mask'])

f = open(hdr['OBJECT'] + '_pos.mas', 'w')
for i in range(len(tbl2)):
    f.write('      ' + '      '.join(
        [str(tbl2['ll_mask_s'][i])[0:14], str(tbl2['ll_mask_e'][i])[0:14], str(tbl2['w_mask'][i])[0:12]]) + '\n')
f.close()
