import os
import glob
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import etienne_tools as et

def get_tapas(request, tapas_path = '/Volumes/courlan/TAPAS/', user = 'artigau@astro.umontreal.ca'):
    # get curret directory
    current_dir = os.getcwd()

    # create folder for run
    if not os.path.isdir(str(request)):
        # create folder for TAPAS run
        os.system('mkdir {0}{1}'.format(tapas_path,request))
        # change to tapas directory
        os.chdir('{0}{1}'.format(tapas_path,request))

        # loop through 12 files in case this is a double request
        for id in range(1,12):
            url = 'ftp://tapas:tapas@ftp.ipsl.fr/{0}/Ether_TAPAS_{1}/tapas_0000{2}.ipac.gz'.format(user,
                                                                                                   request,
                                                                                                   str(id).zfill(2))
            print('wget '+url)
            # get the URL
            os.system('wget {}'.format(url))
        # unzip the ipac files
        os.system('gunzip *.gz')
        # find the ipac files in the folder
        files = glob.glob('tapas_??????.ipac')

        # loop through files and plot the spectra
        for file in files:
            outname = file.split('.')[0]+'.fits'
            print(outname)
            if not os.path.isfile(outname):
                tbl = Table.read(file,format='ascii')
                tbl.write(outname)

                plt.plot(tbl['wavelength'],tbl['transmittance'],label = file)
        plt.xlabel('wavelength')
        plt.ylabel('transmittance')
        plt.legend()
        plt.savefig('{}_transmittance.pdf'.format(request))
        plt.show()

    # back to initial directory
    os.chdir(current_dir)

    return

def tapas2npy():
    # creates the npy and mask files necessary for the telluric preclean

    outname = 'tapas_HARPS.npy'
    dirs = glob.glob('????')

    tapas_wave = []
    trans_water = []
    trans_others = []
    for dir in dirs:
        print(dir)
        for i in tqdm(range(1,7)):
            tbl = Table.read(dir+'/tapas_00000'+str(i)+'.fits')
            if i == 1:
                tapas_wave = np.append(tapas_wave,np.array(tbl['wavelength']))
                trans_water = np.append(trans_water,np.array(tbl['transmittance']))
            if i == 2:
                tmp = np.array(tbl['transmittance'])
            if i >2:
                tmp *= np.array(tbl['transmittance'])

        trans_others = np.append(trans_others,tmp)

    ord = np.argsort(tapas_wave)
    tapas_wave = tapas_wave[ord]
    trans_water = trans_water[ord]
    trans_others = trans_others[ord]

    tmp_tapas = np.zeros([3, len(tapas_wave)])
    tmp_tapas[0] = tapas_wave
    tmp_tapas[1] = trans_others
    tmp_tapas[2] = trans_water
    plt.plot(tapas_wave,trans_others, label = 'others',color = 'orange',alpha = 0.5)
    plt.plot(tapas_wave,trans_water, label = 'water',color = 'blue',alpha = 0.5)
    plt.legend()
    plt.show()
    np.save(outname, tmp_tapas)

    tbl = et.td_convert(mk_tapas_mask(tapas_wave,trans_others))
    tbl.write(outname.split('.')[0]+'_others.csv')
    tbl = et.td_convert(mk_tapas_mask(tapas_wave,trans_water))
    tbl.write(outname.split('.')[0]+'_water.csv')

    return


def mk_tapas_mask(w,f):

    df = np.gradient(f)
    ddf = np.gradient(np.gradient(f))

    line = np.where((np.sign(df[1:]) != np.sign(df[:-1]))
                    & np.isfinite(ddf[1:])
                    & (ddf[1:]>0)
                    & np.isfinite(df[1:])
                    & np.isfinite(df[:-1])
                    & (f[1:]>0.5)
                    & (f[1:]<0.95)
                    )[0]

    # create the output table
    tbl = dict()
    tbl['ll_mask_s'] = np.zeros_like(line, dtype=float)
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
        tbl['ll_mask_s'][i] = wave_cen

    return tbl


