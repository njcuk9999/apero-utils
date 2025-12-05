from astropy.io import fits
from astropy.table import Table
import numpy as np
import glob
from tqdm import tqdm

def sigma(x):
    # difference between -1 and +1 sigma confidence interval
    v1n = np.nanpercentile(x,16)
    v1p = np.nanpercentile(x,84)
    return (v1p-v1n)/2



path = '/cosmos99/nirps/apero-data/nirps_ha_07276_online'
# path to the reduced data
all_files = glob.glob('{}/red/2023-*-*/*_pp_e2dsff_tcorr_A.fits'.format(path))
all_files = np.array(all_files)
all_files = all_files[np.argsort(all_files)]

h = fits.getheader(all_files[0])
wave = fits.getdata('{}/calib/{}'.format(path,h['WAVEFILE']),ext=1)

wave_cen = ['{:.0f}'.format(w) for w in wave[:,2044]]


tbl = Table()
tbl['FILE'] = np.array([f.split('/')[-1] for f in all_files])

# placeholders for the table
tbl['MJDMID'] = np.zeros(len(all_files),dtype = float)
tbl['OBJECT'] = np.zeros(len(all_files),dtype = 'U99')
tbl['DATE-OBS'] = np.zeros(len(all_files),dtype = 'U99')
tbl['JMAG'] = np.zeros(len(all_files),dtype = float)
tbl['EXPTIME'] = np.zeros(len(all_files),dtype = float)
tbl['DELTA_T'] = np.zeros(len(all_files),dtype = float)+np.nan
tbl['OVERHEAD'] = np.zeros(len(all_files),dtype = float)+np.nan

for w in wave_cen:
    tbl['SNR_'+w] = np.zeros(len(all_files),dtype = float)+np.nan
    tbl['SNR_'+w+'_APERO'] = np.zeros(len(all_files),dtype = float)+np.nan
    tbl['PEAK_'+w] = np.zeros(len(all_files),dtype = float)+np.nan


order_Y = 17 # 1.10µm
order_J = 33 # 1.25µm
order_H = 58 # 1.6µm

previous_obj = ''
previous_sp = 0
previous_mjdmid = 0
for i in tqdm(range(len(all_files))):
    h = fits.getheader(all_files[i])
    tbl['DATE-OBS'][i] = h['DATE-OBS']
    tbl['MJDMID'][i] = h['MJDMID']
    tbl['OBJECT'][i] = h['OBJECT']
    tbl['JMAG'][i] = h['HIERARCH ESO OCS TARG JMAG']
    tbl['EXPTIME'][i] = h['EXPTIME']

    sp = fits.getdata(all_files[i],ext=1)
    sp = sp[:,2048-200:2048+200]

    # continuous and same object
    if (h['OBJECT'] == previous_obj) and ((h['MJDMID'] - previous_mjdmid)<1800):
        tbl['DELTA_T'][i] = (h['MJDMID'] - previous_mjdmid)*86400
        overhead = tbl['DELTA_T'][i] - tbl['EXPTIME'][i]

        if overhead>1800:
            continue

        tbl_e2ds = Table.read(all_files[i].replace('_tcorr', ''))

        tbl['OVERHEAD'][i] = overhead

        for iw in range(len(wave_cen)):
            w = wave_cen[iw]

            diff = (sp[iw] / np.nanmedian(sp[iw]) - previous_sp[iw] / np.nanmedian(previous_sp[iw]))
            diff = diff - np.roll(diff, 1)

            tbl['SNR_' + w][iw] = 2/sigma(diff)
            tbl['SNR_' + w + '_APERO'][iw] = np.nanmedian(sp[iw])
            tbl['PEAK_' + w] [iw]=  tbl_e2ds[iw]['SNR']


    previous_obj = h['OBJECT']
    previous_mjdmid = h['MJDMID']
    previous_sp = sp

for iw in range(len(wave_cen)):
    w = wave_cen[iw]
    tbl['SNR_' + w][tbl['SNR_' + w]<0] = np.nan
    tbl['SNR_' + w][tbl['SNR_' + w]>0] = np.nan


tbl.write('all_snr_apero.csv',overwrite=True)


