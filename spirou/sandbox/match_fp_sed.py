import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import etienne_tools as et
from tqdm import tqdm
import glob


sp1 = fits.getdata('2571575o_pp_e2dsff_C.fits')

files = glob.glob('25*o_pp_e2dsff_C.fits')

for file in tqdm(files):
    sp2, hdr = fits.getdata(file, header = True)

    outname = '_sedfit_'.join(file.split('_pp_e2dsff_'))
    #print(outname)

    for iord in tqdm(range(sp1.shape[0]),leave = False):
        tmp1 = np.array(sp1[iord])
        tmp2 = np.array(sp2[iord])

        #fig, ax = plt.subplots(nrows = 2, ncols = 1,sharex = True, sharey = True)
        #ax[0].plot(tmp1/np.nanpercentile(tmp1,90))
        #ax[0].plot(tmp2/np.nanpercentile(tmp2,90),alpha = 0.5)

        for ite in range(5):
            ratio = np.nanpercentile(tmp1, 90) / np.nanpercentile(tmp2, 90)
            tmp2 *= ratio

            index = np.arange(len(tmp1))

            continu =et.lowpassfilter(tmp1-tmp2)

            tmp2 += continu

            sed = et.sed_ratio(tmp1,tmp2)
            tmp2*=sed

            #print(np.nanstd(continu),np.nanstd(sed))
        #ax[1].plot(tmp1/np.nanpercentile(tmp1,90))
        #ax[1].plot(tmp2/np.nanpercentile(tmp2,90),alpha = 0.5)

        #plt.show()
        #plt.plot(tmp1,tmp1-tmp2,'g.')
        #plt.show()

        sp2[iord] = tmp2

    fits.writeto(outname,sp2,hdr,overwrite = True)