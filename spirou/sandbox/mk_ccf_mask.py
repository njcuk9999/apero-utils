import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import os as os
from scipy import constants
import glob
from tqdm import tqdm
import etienne_tools as et

#
# Code to generate a mask that can be used with the DRS. You need to provide a template file name
# and the code finds features and correlates the mask against a model at the proper temperature
# (from header or 3600K if not provided). This gets you the systemic velocity, so you can offset
# your lines to a zero velocity. The code creates both .csv and .mas files. The csv can be read
# with astropy.table, while the .mas files are in the right format for the DRS. You get masks for
# the negative and positive spectroscopic features (_neg and _pos) masks. You also get a _full mask,
# but the DRS does not handle both positive and negative features yet.
#


def mk_ccf_mask(template, doplot = False):
    # CSV table to contain systemic velocities. Will replace the entry of the same object if it is already present
    # in the table, will create the table if it does not exist
    systemic_velocity_table = 'systemic_velo.csv'

    print(template)

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

    if 'FP' not in template:
        # get goettigen models if you don't have them.
        for temperature in np.arange(3000, 6100, 100):
            temperature = str(np.int(np.round(temperature, -2)))
            outname = '{0}/lte0{1}-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(path_to_models,temperature)

            if not os.path.isfile(outname):
                os.system(
                    'wget {0}PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte0{1}-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(ftp_link,temperature))

                os.system('mv lte0{1}-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits {0}'.format(path_to_models,temperature,))



    # read template and header
    tbl, hdr = fits.getdata(template, ext=1, header=True)

    if 'FP' not in template:
        hdr2 = fits.getheader(template, ext=2)
        nsp_input = hdr2['NAXIS2']
    else:
        nsp_input = 0
        hdr['OBJECT'] = 'FP'

    out_pos_name = hdr['OBJECT'].upper() + '_pos.fits'
    if os.path.isfile(out_pos_name):
        print('File {} exists, we skip'.format(out_pos_name))
        return


    if 'FP' not in template:
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
    tbl = dict()
    tbl['ll_mask_s'] = np.zeros_like(line, dtype=float)
    tbl['ll_mask_e'] = np.zeros_like(line, dtype=float)
    # the weight is the second derivative of the flux. The sharper the line,
    # the more weight we give it
    tbl['w_mask'] = ddf[line]
    tbl['value'] = f[line]

    tbl['depth'] = np.zeros_like(tbl['value'])

    tbl['depth'][1:-1] = 1-tbl['value'][1:-1]/((tbl['value'][0:-2]+tbl['value'][2:])/2)

    for i in tqdm(range(len(line))):
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

    weight = tbl['w_mask']

    systemic_velocity = 0

    if 'FP' not in template:
        # create a spline of the model
        model = InterpolatedUnivariateSpline(wave_phoenix, flux_phoenix)

        # assume a 0 velocity and search
        dv0 = 0
        scale = 1.0


        tbl0 = Table(tbl)

        low_contrast = False

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



            if doplot:


                plt.plot(dvs+systemic_velocity, cc,alpha = 0.5)


            systemic_velocity += (-.5 * fit[1] / fit[0])
            print(systemic_velocity)
            scale /= 5.0

            if np.min(cc)/np.max(cc) > 0.95:
                low_contrast = True
                print('not enough ccf contrast, will end after the plot')

        if doplot:
            plt.title('CCF of model SP with target''s line list\nThis gets you the systemic velocity')
            plt.xlabel('Velocity')
            plt.ylabel('Abritrary flux')
            plt.show()
        if low_contrast:
            return

        hdr['SYSVELO'] = systemic_velocity, 'meas. systemic velocity (km/s)'
        hdr['VELOFILE'] = outname, 'model used for SYSVEL cc'

        print('\n\tsystemic velocity : {0:.2f}km/s\n'.format(systemic_velocity))

        if os.path.isfile(systemic_velocity_table) == False:
            tbl_sysvelo = Table()
            tbl_sysvelo['OBJECT'] = [hdr['OBJECT']]
            tbl_sysvelo['SYSTEMIC_VELOCITY'] = [systemic_velocity]
            tbl_sysvelo['MODEL_FILE'] = [hdr['VELOFILE']]

            print('We create {0}'.format(systemic_velocity_table))
            tbl_sysvelo.write(systemic_velocity_table)

        else:
            tbl_old = Table.read(systemic_velocity_table)
            tbl_sysvelo = Table()

            tbl_sysvelo['OBJECT'] = np.append(hdr['OBJECT'],tbl_old['OBJECT'])
            tbl_sysvelo['SYSTEMIC_VELOCITY'] = np.append(systemic_velocity,tbl_old['SYSTEMIC_VELOCITY'])
            tbl_sysvelo['MODEL_FILE'] = np.append(hdr['VELOFILE'],tbl_old['MODEL_FILE'])

            print('We append {0}'.format(systemic_velocity_table))
            tbl_sysvelo.write(systemic_velocity_table, overwrite = True)

    # convert back to table for manipulation
    tbl = et.td_convert(tbl)

    if 'FP' not in template:
        valid = np.isfinite(f)
        spline = InterpolatedUnivariateSpline(w[valid],f[valid], k = 1, ext=0)
        # DETERMINATION OF H-band FWHM
        # cen, ew, amp, zp, slope
        dvs = np.arange(-50000,50000,500)+systemic_velocity*1000
        cc = np.zeros_like(dvs,dtype = float)

        H = (tbl['ll_mask_s'] > 1500) * (tbl['ll_mask_s'] > 1800) * (tbl['w_mask'] > 0)
        wave_H = np.array(tbl['ll_mask_s'][H])
        weights_H = np.array(tbl['w_mask'][H])
        for i in range(len(dvs)):
            cc[i] = np.sum(weights_H*spline(et.doppler(wave_H,-dvs[i])))

        imin = np.nanargmin(cc)
        p0 = [dvs[imin], 4000, np.nanmin(cc) - np.nanmedian(cc), np.nanmedian(cc), 0]

        fit_gau = et.fit_gauss(dvs, cc, p0)
        gfit = et.gauss(dvs, *fit_gau)

        cc /= np.polyval(fit_gau[[4, 3]], dvs)
        gfit /= np.polyval(fit_gau[[4, 3]], dvs)

        print(fit_gau)


        if doplot:
            plt.plot(dvs/1000, cc, color='black', alpha=0.5, label = 'normalized CCF')
            plt.plot(dvs/1000, gfit,alpha = 0.5,label = 'normalized gaussian fit')
            plt.ylabel('flux')
            plt.xlabel('velocity [km/s]')
            plt.legend()
            plt.show()

        hdr['CCF_FWHM'] = np.sqrt(2*np.log(2))*2*fit_gau[1]/1000,'H-band CCF FWHM in km/s'
        hdr['CCF_CONT'] = 1-np.min(cc),'Fractionnal CCF contrast'

    if doplot:



        plt.plot(w,f, 'g-',label = 'input spectrum')
        plt.vlines(tbl[tbl['w_mask'] < 0]['ll_mask_s'], np.nanmin(f), np.nanmax(f), 'k',alpha = 0.2,label = 'positive feature')
        plt.vlines(tbl[tbl['w_mask'] > 0]['ll_mask_s'], np.nanmin(f), np.nanmax(f), 'r',alpha = 0.2,label = 'negative feature')
        plt.legend()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('Arbitrary flux')
        plt.show()

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

    hdu1 = fits.PrimaryHDU()
    hdu1.header['SYSTVEL'] =systemic_velocity,'Systemic velocity'
    hdu1.header['NSPTEMPL'] =nsp_input,'Number of spectra used for tempalte'

    keys_transfer = ['OBJTEMP','PI_NAME','CCF_FWHM','CCF_CONT']
    for key in keys_transfer:
        if key in hdr.keys():
            hdu1.header[key] = hdr[key]

    hdu2 = fits.BinTableHDU(tbl2)
    # convert back from dictionnary to table and save
    new_hdul = fits.HDUList([hdu1, hdu2])
    new_hdul.writeto(out_pos_name, overwrite=True)

templates = glob.glob('Template_s1d_*_sc1d_v_file_AB.fits')
for template in templates:
        mk_ccf_mask(template,doplot = True)