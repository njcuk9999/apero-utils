import numpy as np
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline
from tqdm import tqdm
import os
import warnings
from scipy.signal import medfilt
from scipy.optimize import curve_fit
from scipy import constants
import etienne_tools as et
import pwd
from mk_harps_template import *


if pwd.getpwuid(os.getuid()).pw_name != 'spirou':
    import matplotlib.pyplot as plt


#from fits2wave import fits2wave
#from wave2wave import *
#from lin_mini import *

def gauss(x,x0,ew,amp):
    # gaussian for curve_fit to get the absorption velocities
    return np.exp( -0.5*(x-x0)**2/ew**2)*amp

def get_abso_sp(wave, expo_others, expo_water, spl_others, spl_water, ww = 4.95, ex_gau = 2.20, dv_abso = 0,
                wave0 = 965, wave1 = 2500, dv_grid = 1.0):

    # return an absorption spectrum from exponents describing water and 'others' in absorption
    # inputs :
    # wave -> wavelength grid onto which the spectrum is splined
    # expo_others -> optical depth of all species other than water
    # expo_water -> optical depth of water
    # optional:
    #   ww -> gaussian width of the kernel
    #   ex_gau -> exponent of the gaussian, ex_gau = 2 is a gaussian, ex_gau > 2 is boxy

    if (expo_others == 0)*(expo_water==0):
        # for some tests, one may give 0 as exponents and get just a flat 1
        return np.ones_like(wave)

    # define the convolution kernel for the model. This shape factor can be modified if needed
    w = ww / 2.35 # fwhm of a gaussian of exp = 2.0

    # defining the convolution kernel x grid, defined over 4 fwhm
    kernel_width = int(ww*4)
    dd = np.arange(-kernel_width,kernel_width+1,.5,dtype = float)

    # normalization of the kernel
    ker = np.exp(-0.5 * np.abs(dd / w) ** ex_gau)
    ker = ker[ker > 1e-6 * np.max(ker)] # shorten then kernel to keep only pixels that are more than 1e-6 of peak
    ker /= np.sum(ker) # normalize your kernel

    # create a magic grid onto which we spline our transmission, same as for the s1d_v
    c = (constants.c.value/1000)
    len_magic = int(np.ceil(np.log(wave1 / wave0) * c / dv_grid))
    magic_grid = np.exp(np.arange(len_magic) / len_magic * np.log(wave1 / wave0)) * wave0

    sp_others = spl_others(magic_grid)
    sp_water = spl_water(magic_grid)

    # for numerical stability, we may have values very slightly below 0 from the spline above.
    # negative values don't work with fractionnal exponents
    sp_others[sp_others<0] = 0
    sp_water[sp_water<0] = 0

    # applying optical depths
    trans_others = sp_others**expo_others
    trans_water = sp_water**expo_water

    # getting the full absorption at full resolution
    trans = trans_others*trans_water
    # convolving after product (to avoid the infamous commutativity problem
    trans_convolved = np.convolve(trans, ker, mode='same')

    # spline that onto the input grid and allow a velocity shift
    if wave1<1000:
        # we have HARPS data that is computed in air
        wave_spline =  et.doppler(magic_grid,1000*dv_abso)
        wave_spline /= et.air_index(wave_spline,Unit='nm')
    else:
        wave_spline =  et.doppler(magic_grid,1000*dv_abso)

    spl2 = InterpolatedUnivariateSpline( wave_spline  ,trans_convolved)


    # if this is a 2d array from an e2ds, we loop on the orders
    sz = wave.shape

    if len(sz) == 2:
        out = np.zeros(sz)
        for iord in range(sz[0]):
            out[iord] = spl2(wave[iord])
    else:
        out = spl2(wave)

    return out

def get_tapas_spl(params):
    tapas_npy_file = str(params['tapas_npy_file'])
    dv0 = params['dv0']

    # ----------------------------------------------------------------------------------------------------------------
    # read table with TAPAS model for 6 molecules. We may want to keep it as a numpy array to make it faster.
    # will need to be replaced with numpy readouts consistent with the rest of the DRS
    # ----------------------------------------------------------------------------------------------------------------

    tmp_tapas = np.load(tapas_npy_file)
    tapas_wave = tmp_tapas[0]
    trans_others = tmp_tapas[1]
    trans_water = tmp_tapas[2]

    #define spline function for both, optionally shift the grid
    spl_others = InterpolatedUnivariateSpline(et.doppler(tapas_wave,dv0*1000), trans_others, k=1, ext=3)
    spl_water = InterpolatedUnivariateSpline(et.doppler(tapas_wave,dv0*1000), trans_water, k=1, ext=3)

    return spl_others, spl_water

def clean_oh_pca(image_e2ds, wave_e2ds):
    # provide an e2ds image and corresponding wavelength grid. We load the n sky principal components and
    # fit the derivatives to the derivatives of the e2ds image. As the wavenlength grid changes between
    # nights, we spline the PCs [on the master grid] onth that night's grid.

    principalComponents = fits.getdata('sky_PCs.fits')
    n_components = principalComponents.shape[1]-1

    # load the master wavelength grid
    master_wave = principalComponents[:,0].reshape([49,4088])

    principalComponents = principalComponents[:,1:].reshape([49,4088,n_components])

    # replace NaNs in the science data with zeros to avoid problems in the fitting below
    ribbon_e2ds = image_e2ds.ravel()
    ribbon_e2ds[~np.isfinite(ribbon_e2ds)] = 0


    # make the PCs a ribbon that is N_pc * (4088*49)
    ribbons_pcs = np.zeros([n_components,len(ribbon_e2ds)])

    # lead the PCs and transform to the night grid
    for i in range(n_components):
        tmp = wave2wave(principalComponents[:,:,i], master_wave, wave_e2ds )
        ribbons_pcs[i] = tmp.ravel()

    # output for the sky model
    sky_model = np.zeros_like(ribbon_e2ds)

    # --> here we could have a loop, that's why the sky_model is in a difference within the fitting function.
    #     for now, let's play safe and not have a loop.
    # linear minimisation of the ribbon's derivative to the science data derivative
    amps,model = lin_mini(np.gradient(ribbon_e2ds-sky_model),np.gradient(ribbons_pcs,axis=1))



    # reconstruct the sky model with the amplitudes derived above
    for i in range(n_components):
        sky_model+=ribbons_pcs[i]*amps[i]

    # sky cannot be negative
    sky_model[sky_model<0] = 0

    # return the clean image and sky model
    return image_e2ds - sky_model.reshape([49,4088]), sky_model.reshape([49,4088])

"""
def graceful_error_in_tellu(image_e2ds, wave_e2ds, hdr_e2ds):
    # default parameters for a crash-less exit of the fit_tellu function.
    print('\n\n ERROR ... we have the graceful exit!')
    print(hdr_e2ds['FILENAME'], hdr_e2ds['OBJECT'])
    print('\n\n')

    default_water_abso = 4.0 # typical water abso exponent. Compare to values in headers for high-snr targets later

    flag_error = True

    spl_others, spl_water = get_tapas_spl()


    expo_others = hdr_e2ds['AIRMASS']
    expo_water = default_water_abso

    print(expo_others, expo_water)
    abso_e2ds = get_abso_sp(wave_e2ds, expo_others, expo_water, spl_others, spl_water, ww=5.1, ex_gau=2.3)

    mask = abso_e2ds < np.exp(-1)
    corrected_e2ds = image_e2ds / abso_e2ds
    corrected_e2ds[mask] = np.nan

    # return corrected e2ds, mask, absorption, expo_water, expo_others
    return corrected_e2ds, mask, abso_e2ds, sky_model,  expo_water, expo_others, 0, 0
"""


def hybrid_fit_tellu(params, image_e2ds, wave_e2ds, hdr_e2ds,spl_others, spl_water):

    # Pass an e2ds  image and return the telluric-corrected data. This is a rough model fit and
    # we will need to perform PCA correction on top of it.
    #
    # Will fit both water and all dry components of the absorption separately.
    #
    # Set force_airmass = True to have the dry components forced at the level expected from airmass in hdr. In practice
    # this is fine for most/all purposes. Don't forget, this is a first pass of correction and we'll catch residuals
    # left by inaccuracies here.
    #
    # set doplot = True to see ccf values and before/after correction.
    #
    # Underlying idea: We correct with a super naive tapas fit and iterate until the CCF of the telluric absorption
    # falls to zero. We have 2 degrees of freedom, the dry and water components of the atmosphere. The instrument profile
    # is defined by two additional parameters [ww -> fwhm, ex_gau -> kernel shape] that should be defined at the DRS
    # level: FWHM and kernel shape parameter.
    #
    # Again, this is just a cleanup PRIOR to PCA correction, so if the code is not perfect in it's correction, this is
    # fine as we will empirically determine the residuals and fit them in a subsequent step.
    #
    # we set bounds to the limits of the reasonable domain for both parameters. Reaching the limit will raise an error
    # flag.

    doplot = params['doplot']
    force_airmass = params['force_airmass']
    snr_min = params['snr_min']
    water_bounds = params['water_bounds']
    others_bounds = params['others_bounds']
    clean_oh = params['clean_oh']
    ww = params['ww']
    ex_gau = params['ex_gau']
    return_ccf_power = params['return_ccf_power']
    water_ccf_file = params['water_ccf_file']
    others_ccf_file = params['others_ccf_file']
    wave0 = params['wave0']
    wave1 = params['wave1']
    dv0 = params['dv0'] # zeroth velocity of the wavelength grid
    ccf_scan_range = int(params['ccf_scan_range'])
    template_file = params['template_file']
    mask_domain = params['mask_domain']
    default_water_abso = params['default_water_abso']
    BERV = params['BERV']
    recenter_ccf = params['recenter_ccf']
    berv_offset = params['berv_offset']

    if os.path.isfile(template_file):
        template =  Table(fits.getdata(template_file))
        wave_template = template['wavelength']
        flux_template = template['flux']
        g = np.isfinite(flux_template)
        wave_template = wave_template[g]
        flux_template = flux_template[g]

        template_flag = True
        template_spline = InterpolatedUnivariateSpline(wave_template,flux_template,k=1,ext=3)
    else:
        template_flag = False


    image_e2ds_ini = np.array( image_e2ds)

    if clean_oh: # remove OH lines
        image_e2ds, sky_model = clean_oh_pca(np.array(image_e2ds),wave_e2ds)

    # we ravel the wavelength grid to make it a 1d array of increasing wavelength. We will trim the overlapping domain
    # between orders
    keep = np.ones_like(wave_e2ds)
    if image_e2ds.shape[0] > 3:
        for iord in range(1,wave_e2ds.shape[0]-1):
            keep[iord] = (wave_e2ds[iord]>(wave_e2ds[iord-1][::-1]))*(wave_e2ds[iord]<(wave_e2ds[iord+1][::-1]))
        keep[0] = 0 # firts and last orders are rejected
        keep[-1] = 0

    wave = wave_e2ds.ravel()
    keep = keep.ravel() != 0
    sp = image_e2ds.ravel()
    sp_ini = image_e2ds_ini.ravel()

    # keep non-overlapping bits
    sp = sp[keep]
    sp_ini = sp_ini[keep]
    wave = wave[keep]

    if template_flag:
        if not berv_offset:
            template = template_spline(et.doppler(wave,-BERV*1000) )
        else:
            template = template_spline(et.doppler(wave,0) )

    else:
        template = np.ones_like(wave)




    flag_error = False # to be used for QCing

    if hdr_e2ds['instrume'] != 'HARPS':
        snr = np.array([hdr_e2ds[key] for key in hdr_e2ds['EXTSN*']],dtype = float)
        snr[np.isfinite(snr) == False] = 0

        max_snr = np.nanmax(snr)

        #if np.nanmax(snr) < snr_min:
        #    print('We are at SNR<{0} for all orders and will not be able to fit tellurics'.format(snr_min))
        #    return graceful_error_in_tellu(image_e2ds, wave_e2ds, hdr_e2ds)

        # we are going to mask all orders with an snr below a certain threshold. SNR=10 is fine
        mask_low_snr = np.zeros_like(wave)
        for iord in range(49):
            if snr[iord]<snr_min:
                w1 = np.min(wave_e2ds[iord])
                w2 = np.max(wave_e2ds[iord])
                mask_low_snr[ (wave>w1)*(wave<w2) ] = 1

        # mask bad domain to avoid fitting it
        sp[mask_low_snr == 1] = np.nan

    # for numerical stabiility, remove NaNs. Setting to zero biases a bit the CCF, but this should be OK
    # after we converge
    sp[np.isfinite(sp) == False] = 0
    sp[sp<0] = np.nan


    dv_abso = 0.0
    dv_water = 0.0
    dv_others = 0.0

    # define ploting params if doplot ==  True
    if doplot:
        fig, ax = plt.subplots(nrows=1, ncols=2)

    # species => 0 all absorbers other than water. These files are provided outside of the DRS
    mask_others = Table.read(others_ccf_file)
    keep = (mask_others['ll_mask_s']>mask_domain[0])*(mask_others['ll_mask_s']<mask_domain[1])
    mask_others = mask_others[keep]



    mask_water = Table.read(water_ccf_file)
    keep = (mask_water['ll_mask_s']>mask_domain[0])*(mask_water['ll_mask_s']<mask_domain[1])
    mask_water = mask_water[keep]

    if wave1<1000: # we have HARPS data, express things in air
        mask_water['ll_mask_s'] /= et.air_index(mask_water['ll_mask_s'],Unit='nm')
        mask_others['ll_mask_s'] /= et.air_index(mask_others['ll_mask_s'],Unit='nm')

    if berv_offset:
        mask_water['ll_mask_s'] = et.doppler(mask_water['ll_mask_s'],-BERV*1000)
        mask_others['ll_mask_s'] = et.doppler(mask_others['ll_mask_s'],-BERV*1000)


    # start with no correction of abso to get the ccf
    expo_water = 0
    expo_others = 0 # we start at zero to get a velocity measurement even if we may force to airmass

    # keep track of consecutive exponents and test convergence
    expo_water_prev = np.inf
    if ~force_airmass:
        expo_others_prev = np.inf
    dexpo = np.inf

    amps_water = [] # amplitude of CCF
    if ~force_airmass:
        amps_others = []

    expos_water = []
    if ~force_airmass:
        expos_others = []

    # first guess at velocity of absorption is 0 km/s
    dv_abso = 0.0

    # scanning range for the ccf computation
    dd = np.arange(-ccf_scan_range, ccf_scan_range + 1) #+ dv_abso
    # placeholder for the ccfs
    ccf_others = np.zeros_like(dd, dtype=float)
    ccf_water = np.zeros_like(dd, dtype=float)


    ite = 0 # stop at 20th iteration, normally it converges in ~6

    wave_trans = np.array(wave)
    if berv_offset:
        wave_trans = et.doppler(wave_trans, BERV * 1000)

    while dexpo>1e-4 and (ite<20):

        trans = get_abso_sp(wave_trans,expo_others, expo_water, spl_others, spl_water, ww = ww,
                            ex_gau = ex_gau, dv_abso = dv_abso, wave0 = wave0, wave1 = wave1)
        sp_tmp = sp/trans
        sp_tmp/=template

        valid = np.isfinite(sp_tmp)

        # transmission with the exponenent value
        valid &= (trans > np.exp(-1))

        # apply some cuts to very discrepant points. these will be set to zero not to bias the ccf too much
        cut = np.nanmedian(np.abs(sp_tmp))*10 #>10 sigma outliers
        sp_tmp[np.isfinite(sp_tmp) == False] = 0 # not finite
        sp_tmp[sp_tmp > cut]=0
        sp_tmp[sp_tmp<0]=0 # below zero


        # get the ccf of the test spectrum
        spl = InterpolatedUnivariateSpline(wave[valid],sp_tmp[valid], k=1, ext=1)

        all_others = np.zeros([len(mask_others),len(dd)])
        all_water = np.zeros([len(mask_water),len(dd)])
        for i in range(len(dd)):
            # we compute the ccf_others all the time, even when forcing the airmass, just to look at its structure
            # and potential residuals
            tmp = spl(et.doppler(mask_others['ll_mask_s'],dd[i]*1000))*mask_others['w_mask']
            tmp[tmp ==0]=np.nan
            all_others[:,i] = tmp


            tmp = spl(et.doppler(mask_water['ll_mask_s'],dd[i]*1000))*mask_water['w_mask']
            tmp[tmp ==0]=np.nan
            all_water[:,i] = tmp

        #for  i in range(len(mask_water)):
        #    all_water[i] /=np.nanmedian(all_water[i])
        #for  i in range(len(mask_others)):
        #    all_others[i] /=np.nanmedian(all_others[i])

        ccf_water = np.nansum(all_water, axis=0)
        ccf_others = np.nansum(all_others, axis=0)


        ccf_water -= np.polyval(et.robust_polyfit(dd,ccf_water,1,3)[0],dd)
        ccf_others -= np.polyval(et.robust_polyfit(dd,ccf_others,1,3)[0],dd)

        #ccf_water -= np.median(ccf_water)
        #ccf_others -= np.median(ccf_others)

        # subtract the median of the ccf outside the core of the gaussian. We take this to be the 'external' part of
        # of the scan range
        ccf_water -= np.nanmedian(ccf_water[(np.abs(dd)>ccf_scan_range/2)])
        if ~force_airmass:
            ccf_others -= np.nanmedian(ccf_others[(np.abs(dd)>ccf_scan_range/2)])



        # we measure absorption velocity by fitting a gaussian to the absorption profile. This updates the
        # dv_abso value for the next steps.

        #if np.min(np.isfinite(ccf_water)) == False:
        #    # if CCF is NaN, we have the graceful exit
        #    #return graceful_error_in_tellu(image_e2ds, wave_e2ds,hdr_e2ds)

        if recenter_ccf:
            if ite ==0:
                p0 = [dd[np.argmin(ccf_water)],4, np.nanmin(ccf_water)]
                #plt.plot(dd,ccf_water)
                #plt.show()
                popt, pcov = curve_fit(gauss, dd, ccf_water, p0 = p0)
                dv_water = popt[0]


                if (hdr_e2ds['INSTRUME'] == 'HARPS'):#*(template_flag == False):
                    dv_others = dv_water

                else:
                    p0 = [0,4, np.nanmin(ccf_others)]
                    popt, pcov = curve_fit(gauss, dd, ccf_others, p0 = p0)
                    dv_others = popt[0]
                    #print('dv water = {0} km/s, dv dry = {1} km/s'.format(dv_water, dv_others))
                dv_abso = (dv_others+dv_water)/2.0 # mean of water and others
                #ccf_scan_range/=2


                if ite >2:

                    gfit = gauss(dd,*popt)
                    grad_gfit = np.gradient(gfit)
                    grad_gfit/=np.nansum(grad_gfit**2)


                    ddv = np.nansum(ccf_water*grad_gfit)
                    #print(ddv,dv_abso)
                    dv_abso-=(ddv)


        # get the amplitude of the middle of the CCF
        amp_water = np.nansum(ccf_water[np.abs(dd-dv_abso)<ww])
        if ~force_airmass:
            amp_others = np.nansum(ccf_others[np.abs(dd-dv_abso)<ww])


        # log amplitude of current exponent values
        if ~force_airmass:
            amps_others = np.append(amps_others,amp_others)
            expos_others = np.append(expos_others,expo_others)

        amps_water = np.append(amps_water,amp_water)
        expos_water = np.append(expos_water,expo_water)

        if ite == 0:
            # next exponent to be used
            expo_others = hdr_e2ds['AIRMASS'] # in all cases

            # typical value for water absorption
            expo_water = default_water_abso

        else:
            # we have 2 values or more, if this is iteration<5, we fit a line and find the
            # point where amp would be 0. If we have >5, we get smarter and fit a 2nd order polynomial
            if ite>5:
                if ~force_airmass:
                    ord = np.argsort(np.abs(amps_others))
                    fit_others = np.polyfit(amps_others[ord[0:4]],expos_others[ord[0:4]],1)

                ord = np.argsort(np.abs(amps_water))
                fit_water = np.polyfit(amps_water[ord[0:4]],expos_water[ord[0:4]],1)

            else:
                fit_others = np.polyfit(amps_others,expos_others,1)
                fit_water = np.polyfit(amps_water,expos_water,1)

            # next guess
            if force_airmass:
                expo_others = hdr_e2ds['AIRMASS']  # in all cases
            else:
                expo_others = fit_others[1]

                if expo_others < others_bounds[0]:
                    expo_others = others_bounds[0]
                    flag_error = True

                if expo_others > others_bounds[1]:
                    expo_others = others_bounds[1]
                    flag_error = True

            expo_water = fit_water[1] # find best guess for exponent

            if expo_water<water_bounds[0]:
                expo_water = water_bounds[0]
                flag_error = True # could be used as QC check

            if expo_water > water_bounds[1]:
                expo_water = water_bounds[1]
                flag_error = True # could be used as QC check

            # have we converged yet?
            if force_airmass:
                dexpo = np.abs(expo_water_prev-expo_water)
            else:
                dexpo = np.sqrt((expo_water_prev - expo_water) ** 2 + (expo_others_prev - expo_others) ** 2)

        # keep track of the convergence params
        expo_water_prev = expo_water
        expo_others_prev = expo_others

        if doplot:
            if ite <=1:
                # show CCFs to see if the correlation peaks have been 'killed'
                ax[0].plot(dd,ccf_water,alpha = 0.5)
                ax[0].set(xlabel = 'dv [km/s]', ylabel = 'ccf power', title = 'Water ccf')
                ax[1].plot(dd,ccf_others,alpha=0.5)
                ax[1].set(xlabel = 'dv [km/s]', ylabel = 'ccf power', title = 'Dry ccf')

        #if ite ==0:
        #    dd = np.array(np.arange(-10*ww,10*ww)+dv_abso,dtype = int)
        #    # placeholder for the ccfs
        #    ccf_others = np.zeros_like(dd, dtype=float)
        #    ccf_water = np.zeros_like(dd, dtype=float)


        ite+=1

    if ite==19:
        flag_error = True  # we did not converge

    if doplot:
        # show CCFs to see if the correlation peaks have been 'killed'
        ax[0].plot(dd,ccf_water,alpha = 0.8,color = 'black')
        ax[0].set(xlabel = 'dv [km/s]', ylabel = 'ccf power', title = 'Water ccf')
        ax[1].plot(dd,ccf_others,alpha=0.8,color = 'black')
        ax[1].set(xlabel = 'dv [km/s]', ylabel = 'ccf power', title = 'Dry ccf')


    if doplot:
        # debug plot to show absorption spectrum
        plt.tight_layout()
        plt.show()

        #plt.clf()
        #plt.plot(dd, ccf_water/np.max(np.abs(ccf_water)))
        #plt.plot(dd, grad_gfit/np.max(np.abs(grad_gfit)), alpha=0.5)
        #plt.show()


        med = np.nanmedian(sp)
        sp/=med
        sp_ini/=med

        mask = np.ones_like(trans)
        mask[trans<np.exp(-1)] = np.nan
        mask[~np.isfinite(sp)] = np.nan

        scale = np.nanpercentile(sp/trans*mask,99.5)

        label =  'input - OH'
        if not clean_oh:
            label = 'input'
        plt.plot(wave,sp/scale,color = 'red',label =label)

        if clean_oh:
            plt.plot(wave,sp_ini/(trans*mask)/scale , color='magenta',alpha = .5, label='input')

        label = '(input - OH)/abso'
        if not clean_oh:
            label = 'input/abso'

        plt.plot(wave,sp/trans*mask/scale,color = 'green',label =label )
        plt.plot(wave,trans,color='orange',label = 'abso', alpha = 0.5)

        if template_flag:
            plt.plot(wave,sp/(template/np.nanmedian(template)),color='cyan',alpha = 0.5,label = 'sp/template')

        plt.legend()
        plt.title(hdr_e2ds['OBJECT'])
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Normalized flux \n transmission')
        plt.ylim([0,np.nanmax(sp_ini/(trans*mask)/scale)*1.05])
        plt.tight_layout()
        plt.show()


    wave_trans = np.array(wave_e2ds)
    if berv_offset:
        wave_trans = et.doppler(wave_trans, BERV * 1000)


    # get the final absorption spectrum to be used on the science data. No trimming done on the wave grid
    abso_e2ds = get_abso_sp(wave_trans, expo_others, expo_water, spl_others, spl_water, ww = ww,
           ex_gau = ex_gau, dv_abso = dv_abso, wave0 = wave0, wave1 = wave1)#, dv0 = dv0)

    # all absorption deeper than exp(-1) ~ 30% is considered too deeep to be corrected. We set values there to NaN
    mask = abso_e2ds < np.exp(-1)
    abso_e2ds[mask] = np.nan

    if not clean_oh:
        sky_model = np.zeros_like(image_e2ds_ini)

    corrected_e2ds = (image_e2ds_ini-sky_model)/abso_e2ds
    corrected_e2ds[np.isfinite(image_e2ds_ini) ==  False] = np.nan

    if return_ccf_power == False:
        return corrected_e2ds, mask, abso_e2ds, sky_model,  expo_water, expo_others, dv_water, dv_others

    else:
        print('debug mode')
        # this is a debug mode. We return the CCF power just to better adjust the shape of the LSF
        return np.nansum(np.gradient(ccf_water)**2),np.nansum(np.gradient(ccf_others)**2)



def process_harps_tellu(params, files):

    files = np.array(files)

    # get the tapas splines, no need to re-compute it at each iteration. We could be even smarter and do it above
    # this function
    spl_others, spl_water = get_tapas_spl(params)

    for i in tqdm(range(len(files))):
        file = files[i]

        if 's1d' in file:
            params['berv_offset'] = True
        else:
            params['berv_offset'] = False

        outdir = '/'.join(file.split('/')[0:-1])+'-tc'
        if not os.path.isdir(outdir):
            cmd = 'mkdir '+outdir
            print(cmd)
            os.system(cmd)
        outname = outdir+'/'+file.split('/')[-1]

        if os.path.isfile(outname)*(params['force'] == False):
            continue

        if 'e2ds' in file:
            image_e2ds, hdr_e2ds = fits.getdata(file,header = True)
            hdr_e2ds = et.harps2spirou(hdr_e2ds)
            wave_e2ds = et.fits2wave(hdr_e2ds)

        else:
            image_e2ds, hdr_e2ds = fits.getdata(file, header=True)
            hdr_e2ds = et.harps2spirou(hdr_e2ds)
            wave_e2ds = (np.arange(len(image_e2ds)) * hdr_e2ds['CDELT1'] + hdr_e2ds['CRVAL1']) / 10.0

            image_e2ds = np.reshape(image_e2ds, [1, len(image_e2ds)])
            wave_e2ds = np.reshape(wave_e2ds, [1, len(wave_e2ds)])

        params['BERV'] = hdr_e2ds['BERV'] # only relevant for e2ds in HARPS


        if False:
            params['return_ccf_power'] = True
            params['doplot'] = False


            if True:
                params['ww'] =1.4
                params['ex_gau'] = 2.2

                for ww in np.arange(1.0,2.0,.1):
                    params['ww'] = ww
                    ccf1, ccf2 =  hybrid_fit_tellu(image_e2ds, wave_e2ds, hdr_e2ds, params=params)
                    plt.plot(ww,ccf1,'go')
                    print(ww,ccf1,ccf2)
                plt.show()

                params['ww'] =1.4
                params['ex_gau'] = 2.2

                for ex_gau in np.arange(0.5,3.0,.2):
                    params['ex_gau'] = ex_gau
                    ccf1, ccf2 =  hybrid_fit_tellu(image_e2ds, wave_e2ds, hdr_e2ds, params=params)
                    plt.plot(ex_gau,ccf1,'go')
                    print(ex_gau,ccf1,ccf2)
                plt.show()

        if True:
            corrected_e2ds, mask, abso_e2ds, sky_model,  expo_water, expo_others, dv_water, dv_others = \
                hybrid_fit_tellu(params, image_e2ds, wave_e2ds, hdr_e2ds,spl_others, spl_water)

            hdr_e2ds['TAU_H2O'] = expo_water
            hdr_e2ds['TAU_OTHE'] = expo_others

            if 's1d' in outname:
                corrected_e2ds = corrected_e2ds.ravel()
            fits.writeto(outname, corrected_e2ds, hdr_e2ds, overwrite=True)
        else:
            print('err {}'.format(outname))


def batch(yaml_params='tellu_clean_params.yaml'):
    objs = ['GL699','EPSINDIA','PROXIMA']

    params = et.load_yaml(yaml_params)
    path = params['path_to_lbl_folder']

    params['doplot'] = False
    for obj in objs:
        template = ''
        params['template_file'] = template
        params['force'] = True
        process_harps_tellu(params, glob.glob(path+'s1d/{}/*.fits'.format(obj)))

        for ite in range(1):
            template = path+'templates/Template_{}-tc_HARPS.fits'.format(obj)
            params['template_file'] = template

            mk_harps_template(glob.glob(path+'s1d/{}/*.fits'.format(obj)),template,obj,3500,force = True)
            process_harps_tellu(params,glob.glob(path+'s1d/{}/*.fits'.format(obj)))

        process_harps_tellu(params, glob.glob(path+'e2ds/{}/*.fits'.format(obj)))



