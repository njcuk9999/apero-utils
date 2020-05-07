import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
from scipy.signal import medfilt

template = 'Template_s1d_Gl699_sc1d_v_file_AB.fits'
# template = 'Template_s1d_Gl15A_sc1d_v_file_AB.fits'
# template = 'Template_s1d_HD189733_sc1d_v_file_AB.fits'

dv = 0.00  # km/s -- width of the CCF box
c = 2.99792458e5

# read wavelength and flux. The wavelength is expressed in Ang, we convert to µm
wave_phoenix = fits.getdata('WAVE_PHOENIX-ACES-AGSS-COND-2011.fits') / 10

# get goettigen models. Use only you don't have the models locally
if False:
    import os as os

    for temperature in np.arange(3000, 6100, 100):
        temperature = str(np.int(np.round(temperature, -2)))
        print(temperature)
        os.system(
            'wget ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte0' + temperature + '-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')

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
model_file = 'lte0' + temperature + '-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
print('Model file = ', model_file)
flux_phoenix = fits.getdata(model_file)

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
    plt.show()

minpos = np.argmin(cc)
fit = np.polyfit(dvs[minpos - 1:minpos + 2], cc[minpos - 1:minpos + 2], 2)

systemic_velocity = -.5 * fit[1] / fit[0]

hdr['SYSVELO'] = systemic_velocity, 'meas. systemic velocity (km/s)'
hdr['VELOFILE'] = model_file, 'model used for SYSVEL cc'

print('systemic velocity : ', systemic_velocity, 'km/s')

plt.plot(w,f, 'g-.')
plt.vlines(tbl[tbl['w_mask'] < 0]['ll_mask_s'], np.nanmin(f), np.nanmax(f), 'k')
plt.vlines(tbl[tbl['w_mask'] > 0]['ll_mask_s'], np.nanmin(f), np.nanmax(f), 'r')
plt.show()


corrv = np.sqrt((1 + systemic_velocity / c) / (1 - systemic_velocity / c))

# updating the table to account for systemic velocity of star
tbl['ll_mask_s'] = tbl['ll_mask_s'] / corrv
tbl['ll_mask_e'] = tbl['ll_mask_e'] / corrv

# ll_mask_s', 'll_mask_e', 'w_mask


# write the output table
fits.writeto(hdr['OBJECT'] + '.fits', tbl, hdr, overwrite=True)

tbl[tbl['w_mask'] < 0].write(hdr['OBJECT'] + '_pos.mas', format='ascii', overwrite=True)
tbl[tbl['w_mask'] > 0].write(hdr['OBJECT'] + '_neg.mas', format='ascii', overwrite=True)
tbl.write(hdr['OBJECT'] + '_full.mas', format='ascii', overwrite=True)




