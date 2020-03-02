import numpy as np
from astropy.io import fits
from scipy.signal import medfilt
import matplotlib.pyplot as plt
import os
import glob

# master DARK_FP files created by mk_master_DARK_MASTER
master_AB = fits.getdata('master_AB.fits')
master_C = fits.getdata('master_C.fits')

# files to be corrected. Do *not* correct tcorr files, this
# correction beens to happen before telluric correction. Only
# correct files with the "_pp_e2dsff_AB.fits" suffix.
# We may create a similar code for _A and _B fiber setups.
# it would follow the same logic.
all_files = glob.glob('2279068a_pp_e2dsff_AB.fits')

for fic_AB in all_files:
    # check that the file is really what we think it should be

    if "_pp_e2dsff_AB" not in fic_AB:
        print('bad bad bad person! {0} is not the right file type as it does ')
        print('contain the "_pp_e2dsff_AB" string'.format(fic_AB))
        continue

    # reading file
    sci_AB,hdr = fits.getdata(fic_AB,header = True)

    # just checking that the file has not been corrected yet. Correcting
    # more than once would be bad
    if 'CORRDARK' in hdr:
        print('File {0} already corrected'.format(fic_AB))
        continue

    # add a flag to header to know that the output file has been corrected
    hdr['CORRDARK'] = True,'corrected dark_fp'

    # look for the corresponding "C" fiber file. It should be in the same
    # directory or something is from
    fic_C = 'C'.join(fic_AB.split('AB'))

    if os.path.isfile(fic_C) == False:
        print('You do it on purpose, {0} is missing from the'.format(fic_C))
        print('folder where {0} is located, naughty person!'.format(fic_AB))
        continue

    # read C fiber image
    sci_C =  fits.getdata(fic_C)

    # for each order, find the normalization factor that would
    # take to C fiber flux in the reference frame and scale it to the
    # observed C fiber that accompanies the AB science data.
    ratio = np.zeros(49)
    for iord in np.arange(49):
        sci_C[iord] -= np.nanpercentile(sci_C[iord],5)

        # only perform the measurement of the amplitude of the persistence signal
        # on the 1-99th percentile. This allows for a small number of hot/dark pixels
        # along the order. Without this, we end up with some spurious amplitude values
        # in ~10% of the frames (i.e., 1 in ~500 orders)
        # expect numpy to complain as there are some NaNs in the sci_C and we apply
        # a > and a <

        mask = sci_C[iord]>np.nanpercentile(sci_C[iord],1)
        mask &= sci_C[iord]<np.nanpercentile(sci_C[iord],99)
        mask &= master_C[iord]>np.nanpercentile(master_C[iord],1)
        mask &= master_C[iord]<np.nanpercentile(master_C[iord],99)

        # approx ratio. We know that frames were normalized with their 90th percentile prior to
        # median combining
        approx_ratio = 1/np.nanpercentile(sci_C[iord],90)

        # much more accurate ratio from a dot product
        ratio[iord] = np.nansum(master_C[iord][mask]*sci_C[iord][mask])/np.nansum(sci_C[iord][mask]**2)

        # must be within 0.9-1.1 of the approximate but robust ratio to be kept, else we use the
        # approx_ratio
        if (ratio[iord]/approx_ratio) < 0.9 or (ratio[iord]/approx_ratio > 1.1):
            print('QC warning, we have a spurious FP_C ratio and use the robust but approx one for order {0}'.format(iord))
            ratio[iord] = approx_ratio


    for iord in np.arange(49):
        # scale the leakage for that order to the observed amplitude of
        sci_AB[iord] -= master_AB[iord]/ratio[iord]
        
    # to be kept in the header as a measure of FP versus science flux ratio
    ratio_leak = np.zeros(49)
    for iord in np.arange(49):
        ratio_leak[iord] = np.nanpercentile(sci_C[iord],90)/np.nanmedian(sci_AB[iord])

    # overwrite the AB file with the corrected from
    fits.writeto('x_'+fic_AB,sci_AB,hdr,overwrite = True)
