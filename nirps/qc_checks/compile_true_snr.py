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

# path to the reduced data
all_files = glob.glob('/cosmos99/nirps/apero-data/nirps_he_07276_online/red/2023-*-*/*_pp_e2dsff_tcorr_A.fits')
all_files = np.array(all_files)
all_files = all_files[np.argsort(all_files)]

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

tbl['SNR_Y'] = np.zeros(len(all_files),dtype = float)+np.nan
tbl['SNR_Y_APERO'] = np.zeros(len(all_files),dtype = float)+np.nan
tbl['PEAK_Y'] = np.zeros(len(all_files),dtype = float)+np.nan

tbl['SNR_J'] = np.zeros(len(all_files),dtype = float)+np.nan
tbl['SNR_J_APERO'] = np.zeros(len(all_files),dtype = float)+np.nan
tbl['PEAK_J'] = np.zeros(len(all_files),dtype = float)+np.nan

tbl['SNR_H'] = np.zeros(len(all_files),dtype = float)+np.nan
tbl['SNR_H_APERO'] = np.zeros(len(all_files),dtype = float)+np.nan
tbl['PEAK_H'] = np.zeros(len(all_files),dtype = float)+np.nan

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

        #### Y band ####
        diff = (sp[order_Y] / np.nanmedian(sp[order_Y]) - previous_sp[order_Y] / np.nanmedian(previous_sp[order_Y]))
        diff = diff-np.roll(diff,1)

        # factor 2 as this is a point-to-point difference of a difference of spectra
        tbl['SNR_Y'][i] = 2/sigma(diff)
        tbl['PEAK_Y'][i] = np.nanmedian(sp[order_Y])
        tbl['SNR_Y_APERO'][i] = tbl_e2ds[order_Y]['SNR']

        #### J band ####
        diff = (sp[order_J] / np.nanmedian(sp[order_J]) - previous_sp[order_J] / np.nanmedian(previous_sp[order_J]))
        diff = diff-np.roll(diff,1)

        # factor 2 as this is a point-to-point difference of a difference of spectra
        tbl['SNR_J'][i] = 2/sigma(diff)
        tbl['SNR_J_APERO'][i] = tbl_e2ds[order_J]['SNR']

        tbl['PEAK_J'][i] = np.nanmedian(sp[order_J])

        #### H band ####
        diff = (sp[order_H] / np.nanmedian(sp[order_H]) - previous_sp[order_H] / np.nanmedian(previous_sp[order_H]))
        diff = diff-np.roll(diff,1)

        # factor 2 as this is a point-to-point difference of a difference of spectra
        tbl['SNR_H'][i] = 2/sigma(diff)
        tbl['SNR_H_APERO'][i] = tbl_e2ds[order_H]['SNR']

        tbl['PEAK_H'][i] = np.nanmedian(sp[order_H])

    previous_obj = h['OBJECT']
    previous_mjdmid = h['MJDMID']
    previous_sp = sp

tbl['SNR_Y'][tbl['SNR_Y']<0] = np.nan
tbl['SNR_J'][tbl['SNR_J']<0] = np.nan
tbl['SNR_H'][tbl['SNR_H']<0] = np.nan

tbl['SNR_Y'][tbl['SNR_Y']>1000] = np.nan
tbl['SNR_J'][tbl['SNR_J']>1000] = np.nan
tbl['SNR_H'][tbl['SNR_H']>1000] = np.nan


tbl.write('all_snr_apero.csv',overwrite=True)


