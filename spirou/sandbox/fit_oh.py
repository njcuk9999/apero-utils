import numpy as np
from astropy.io import fits
import os
import glob
import etienne_tools as et

#### TO BE UPDATED #####
# one needs e2dsff_AB and e2dsff_recon_AB in the directory
files = glob.glob('2516853o_pp_e2dsff_AB.fits')
pc_file = '/Users/eartigau/smart_sky/sky_PCs.fits'

# constants
debug = True # more plots, makes it slow
# number of bright OH lines that will be individually adjusted in amplitude. Done
# only on lines that are at an SNR > 1
nbright = 300

#### END OF UPDATES #####


if debug:
    import matplotlib.pyplot as plt


ifile = 0
for file in files:
    print('[{0}/{1}]'.format(ifile,len(files)))
    ifile+=1

    outname = '_tcorr_AB'.join(file.split('_AB'))

    abso_name = '_recon_AB'.join(file.split('_AB'))

    if os.path.isfile(abso_name) == False:
        print('file {0} does not exist'.format(abso_name))
        continue

    if os.path.isfile(outname) == False:
        abso,hdr = fits.getdata(abso_name, header = True)


        # e-width of the region over which we measure the residuals of brighter lines
        # to adjust them
        FWHM_PIXEL_PSF = 2.1
        ew_weight = 2.5*FWHM_PIXEL_PSF
        # region of which we will compute the weight falloff of a bright sky line
        width = np.int(ew_weight*4)


        ### exists in the DRS ###

        # read science data
        sp = fits.getdata(file)
        print(file,hdr['OBJECT'],hdr['DATE'])

        # get sci wavelength
        wave = et.fits2wave(hdr)

        # get the sky PCs
        pcs = fits.getdata(pc_file)

        # reconstruct the PC wavelength solution
        wavepc = np.reshape(pcs[:,0],sp.shape)

        # get only the sky PCs
        pcs = pcs[:,1:]

        sp2 = et.wave2wave(np.array(sp),wave,wavepc)

        # mimicking the DRS behavior
        print('performing reconstruction')
        amps,recon = et.lin_mini(np.gradient(sp2.ravel()),np.gradient(pcs,axis=0))

        skymodel = np.zeros_like(sp.ravel())
        for iamp in range(len(amps)):
            skymodel+=(pcs[:,iamp]*amps[iamp])
        skymodel[skymodel<0] = 0

        sp2 = sp2.ravel()

        ### just for debug
        skymodel0 = np.array(skymodel)
        skymodel0 = et.wave2wave(skymodel0.reshape(sp.shape),wavepc,wave)


        ### NEW STUFF ###

        # sky amplitude correction
        amp_sky = np.ones_like(skymodel)

        # weight vector to have a seamless falloff of the sky weight
        weight = np.exp(-0.5*np.arange(-width+.5,width+.5)**2/ew_weight**2)

        # mask to know where we looked for a bright line
        mask = np.zeros_like(sp2)

        # keep a mask of what has actually been masked
        mask_plot = np.zeros_like(sp2)+np.nan

        # number of masked lines
        masked_lines = 0

        # loop through bright lines
        for i in range(nbright):
            # find brightest sky pixel that has not yet been looked at
            imax = np.nanargmax(skymodel+mask)
            # keep track of where we looked
            mask[imax-width:imax+width] = np.nan

            # segment of science spectrum minus current best guess of sky
            tmp1 = (sp2-skymodel*amp_sky)[imax-width:imax+width]
            # segment of sky sp
            tmp2 = (skymodel*amp_sky)[imax-width:imax+width]

            # find rms of derivative of science vs sky line
            snr_line = (et.nanstd(np.gradient(tmp2))/et.nanstd(np.gradient(tmp1)))

            # if above 1 sigma, we adjust
            if snr_line>1:
                # dot product of derivative vs science sp
                amp = np.nansum(np.gradient(tmp1)*np.gradient(tmp2)*weight**2)/np.nansum(np.gradient(tmp2)**2*weight**2)
                if amp<(-1):
                    amp = 0
                # modify the amplitude of the sky
                amp_sky[imax - width:imax + width] *= (amp*weight+1)
                mask_plot[imax-width:imax+width] = 0
                masked_lines +=1
        # print how many lines were masked
        print('N masked = {0}'.format(masked_lines))

        # update sky model
        skymodel *= amp_sky

        # reshape sky model
        skymodel = np.reshape(skymodel, sp.shape)

        # back to science wavelength grid
        skymodel = et.wave2wave(skymodel,wavepc,wave)

        if debug:
            plt.plot(wave.ravel(),sp.ravel(),label = 'input',color = 'orange',alpha =0.3)
            plt.plot(wave.ravel(),sp.ravel() - skymodel0.ravel(),label = 'PCA sky model',color = 'red',alpha = 0.9)
            plt.plot(wave.ravel(),sp.ravel() - skymodel.ravel(),label = 'PCA+per-line',color = 'green',alpha = 0.9)

            plt.plot(wave.ravel(),mask_plot,color = 'black',label = 'domain line-by-line')
            plt.title(hdr['object'])
            plt.legend()
            plt.show()

        ### END OF NEW STUFF ###
        sp = (sp-skymodel)/abso

        hdr['OHMASKED'] = masked_lines,'Number of bright OH lines adjusted'
        fits.writeto(outname,sp,hdr,overwrite = True)
    else:
        print('File {0} exists '.format(outname))
