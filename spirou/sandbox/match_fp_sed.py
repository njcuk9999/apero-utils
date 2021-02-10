import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import etienne_tools as et
from tqdm import tqdm
import glob


sp1 = fits.getdata('2571513o_pp_e2dsff_C.fits')
files = glob.glob('2*o_pp_e2dsff_C.fits')

#2571513o_pp_e2dsff_C.fits          	1.57
#2571514o_pp_e2dsff_C.fits          	3.03

for file in tqdm(files):
    sp2, hdr = fits.getdata(file, header = True)

    outname = '_sedfit_'.join(file.split('_pp_e2dsff_'))
    #print(outname)

    for iord in tqdm(np.arange(0,49),leave = False):
        tmp1 = np.array(sp1[iord])
        tmp2 = np.array(sp2[iord])

        #fig, ax = plt.subplots(nrows = 2, ncols = 1,sharex = True, sharey = True)
        #ax[0].plot(tmp1/np.nanpercentile(tmp1,90))
        #ax[0].plot(tmp2/np.nanpercentile(tmp2,90),alpha = 0.5)

        for ite in range(2):
            sed = et.sed_ratio(tmp1,tmp2)
            tmp2*=sed

            continu = et.lowpassfilter( tmp1-tmp2,width = 101)
            tmp2 += continu
            #print(et.nanstd(continu[np.isfinite(tmp1+tmp2)]),et.nanstd(sed[np.isfinite(tmp1+tmp2)]))

        sp2[iord] = tmp2

    #plt.legend()
    #plt.xlabel('reference flux')
    #plt.ylabel('flux - ref')
    #plt.show()

    fits.writeto(outname,sp2,hdr,overwrite = True)
    #fits.writeto('diff.fits',sp2-sp1,hdr,overwrite = True)
    #stop