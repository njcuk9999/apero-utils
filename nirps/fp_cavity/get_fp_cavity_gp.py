import os
import glob
import numpy as np
from astropy.table import Table
from astropy.io import fits
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.constants import c


outname = '/cosmos99/nirps/apero-data/nirps_he_07276_online/calib/cavity_times.rdb'

# TODO : commands to be run on maestria/rali
if False:
    wave_files = np.array(glob.glob('/cosmos99/nirps/apero-data/nirps_he_07276_online/calib/*wave*A.fits'))

    tbl = Table()
    tbl['FILENAME'] = [os.path.basename(f) for f in wave_files]
    tbl['WAVETIME'] = np.zeros_like(tbl['FILENAME'], dtype=float)
    tbl['WCAV_PED'] = np.zeros_like(tbl['FILENAME'], dtype=float)
    tbl['WCAV000'] = np.zeros_like(tbl['FILENAME'], dtype=float)

    for ifile in tqdm(range(len(wave_files)),leave = False):
        h = fits.getheader(wave_files[ifile])
        tbl['WAVETIME'][ifile] = h['WAVETIME']
        tbl['WCAV_PED'][ifile] = h['WCAV_PED']
        tbl['WCAV000'][ifile] = h['WCAV000']
    tbl['CAVITY'] = tbl['WCAV_PED'] + tbl['WCAV000']
    tbl['DV'] = (tbl['CAVITY']/np.nanmedian(tbl['CAVITY'])-1)*c.value
    tbl = tbl[np.argsort(tbl['WAVETIME'])]

    tbl.write(outname, format='ascii.rdb', overwrite=True)


os.system('scp spirou@maestria:{} .'.format(outname))
tbl = Table.read(os.path.basename(outname))

plt.plot(tbl['WAVETIME'], tbl['DV'], 'k.',alpha = 0.5)
plt.show()