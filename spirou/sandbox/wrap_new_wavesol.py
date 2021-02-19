import glob
from astropy.io import fits
import numpy as np
from get_cavity import *
from astropy.table import Table
from tqdm import tqdm


if False:
    files_hc = np.array(glob.glob('*c_pp_e2dsff_AB_wave_hclines_AB.fits'))
    files_hc = files_hc[np.argsort(files_hc)]

    mjdates_hc = np.zeros_like(files_hc, dtype = float)
    files_fp = np.array(glob.glob('*a_pp_e2dsff_AB_wave_fplines_AB.fits'))
    mjdates_fp = np.zeros_like(files_hc, dtype = float)

    for i in tqdm(range(len(files_fp))):
        hdr = fits.getheader(files_fp[i], ext=1)
        mjdates_fp[i] = hdr['MJDATE']

    cavity_length = np.zeros_like(mjdates_hc)
    SBCFPI_T = np.zeros_like(mjdates_hc)
    for i in range(len(files_hc)):

        hdr = fits.getheader(files_hc[i],ext=1)
        mjdates_hc[i] = hdr['MJDATE']
        SBCFPI_T[i] = hdr['SBCFPI_T']

        id_fp = np.argmin(np.abs(mjdates_hc[i] -mjdates_fp ))

        print(files_hc[i],files_fp[id_fp])
        print(np.abs(mjdates_hc[i] - mjdates_fp[id_fp]))

        hcs = Table.read(files_hc[i])
        fps = Table.read(files_fp[id_fp])

        wave_fit,wave_map,cavity = get_wave_sol(hcs,fps,cavity_table_name = 'cavity.csv')
        cavity_length[i] =cavity

        hdu = fits.PrimaryHDU()
        hdr2 = hdu.header

        flag = False

        for key in hdr.keys():
            if key == 'FILENAME':
                flag = True
            if flag:
                hdr2[key] = hdr[key]

        for ii in range(len(wave_fit.ravel())):
            hdr2['WAVE'+str(ii).zfill(4)] = wave_fit[:,::-1].ravel()[ii]
        hdr2['WAVEORDN'] = wave_fit.shape[0]
        hdr2['WAVEDEGN'] = wave_fit.shape[1]-1

        fits.writeto(files_hc[i].split('_')[0]+'_wavesol.fits', wave_map,hdr2, overwrite = True)


    tbl = Table()
    tbl['mjdates'] = mjdates_hc
    tbl['cavity'] = cavity_length
    tbl['SBCFPI_T'] = SBCFPI_T

    tbl.write('cavity_x.csv',overwrite = True)


    scale_ratio = np.polyval(np.array(Table.read('cavity.csv')['coeffs']),1500)/constants.c

    tbl = Table.read('cavity_x.csv')

    plt.close()
    plt.plot(tbl['mjdates'],tbl['cavity']/scale_ratio,'g.')
    plt.xlabel('date [mjd]')
    plt.ylabel('equivalent drift [m/s]')
    plt.show()


files = np.array(glob.glob('2*tcorr_AB.fits'))
for i in tqdm(range(len(files))):
    if int(files[i][0:3])<245:
        files[i] = ''
files = files[files != '']

for i in tqdm(range(len(files))):

    file = files[i]
    if 'wavesol' in file:
        continue

    sp,hdr = fits.getdata(file,header = True)

    if True:#i ==0:
        wavefile = '/Users/eartigau/nirps_wave/data/'+hdr['WFP_FILE'].split('_')[0] + '_wavesol.fits'
        if not os.path.isfile(wavefile):
            print('err...')
            continue

        hdr2 = fits.getheader(wavefile)
        hdr2 = hdr2['WAVE*']

    for key in hdr2.keys():
        hdr[key] = hdr2[key]

    fits.writeto('tcorr_wavesol001'.join(file.split('tcorr')), sp,hdr,overwrite = True)

