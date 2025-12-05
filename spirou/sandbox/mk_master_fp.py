import numpy as np
import glob
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from tqdm import tqdm

outname = '/Users/eartigau/lbl/templates/Template_s1d_FP_sc1d_v_file_AB.fits'
maskname = '/Users/eartigau/lbl/templates/FP_pos.csv'


files = np.array(glob.glob('*C.fits'))
WFPDRIFT = np.zeros_like(files, dtype = float)+np.nan

low_bound, high_bound = -0.92, -0.84

for i in range(len(files)):
    try:
        hdr = fits.getheader(files[i],ext = 1)
        WFPDRIFT[i] = hdr['WFPDRIFT']
    except:
        print('err')

keep = (WFPDRIFT>low_bound)*(WFPDRIFT<high_bound)
files = files[keep]

ngood = 0
all = []
for file in tqdm(files):
    try:
        tmp = fits.getdata(file)
    except:
        continue
    flux = tmp['flux']/np.nanpercentile(tmp['flux'],95)

    if ngood == 0:
        cube = np.zeros([11,len(flux)])

    cube[ngood % 11] += flux
    ngood+=1

for ii in range(11):
    cube[ii]/=np.nanpercentile(cube[ii],95)


tbl = Table()
tbl['wavelength'] = tmp['wavelength']
tbl['flux'] = np.nanmedian(cube, axis = 0)

hdu1 = fits.PrimaryHDU()
hdu2 = fits.BinTableHDU(tbl)
new_hdul = fits.HDUList([hdu1, hdu2])
new_hdul.writeto(outname, overwrite=True)






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


for i in range(len(line)):
    # we perform a linear interpolation to find the exact wavelength
    # where the derivatives goes to zero
    wave_cen = (np.polyfit(df[line[i]:line[i] + 2], w[line[i]:line[i] + 2], 1))[1]

    # we offset that wavelength by the systemic velocity and subtract
    # half of the line width
    tbl['ll_mask_s'][i] = wave_cen
    tbl['ll_mask_e'][i] = wave_cen

tbl[tbl['w_mask']>0].write(maskname, format='ascii', overwrite=True)
