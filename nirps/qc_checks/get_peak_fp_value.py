# code to get the peak FP value expressed in ADU (not ADU/s)
from astropy.io import fits
from astropy.table import Table
import numpy as np
import glob
from tqdm import tqdm

path = '/cosmos99/nirps/raw-data/nirps_ha/2023-05-*/*.fits'
files = np.array(glob.glob(path))

for ifile in range(len(files)):
    if 'FAKE' in files[ifile]:
        continue
    h = fits.getheader(files[ifile])
    if 'FP,FP' not in h['OBJECT']:
        continue

    #print(files[ifile])
    im = fits.getdata(files[ifile],ext=1)
    p995 = np.nanpercentile(im,99.9)
    nd = h['HIERARCH ESO INS FILT1 ID']
    print('{:5.2f}, {:8.2f},{}, {}, {}, {}, {}'.format(p995,p995*h['EXPTIME'],nd,h['EXPTIME']+5.5,h['OBJECT'],
                                                       h['DATE'],files[ifile].split(
        '/')[
        -1]))
