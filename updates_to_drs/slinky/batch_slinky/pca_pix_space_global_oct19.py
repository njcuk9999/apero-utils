import numpy as np
import glob
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
from etienne_tools import hdr2wave, doppler, lowpassfilter, sigma
import getpass
import warnings
from tqdm import tqdm
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from wpca import EMPCA
import os
from scipy.ndimage import shift

def offset_vector(im, dx):
    """
    :param im: Input 2D image with spectra an N offsets to be applied
    :param dx: offsets in pixel space
    :return: 2D image with spectra shifted by dx
    """

    # initialize output
    im2 = np.zeros_like(im)
    # loop over images
    for i in range(im.shape[0]):
        im2[i,:] = shift(im[i,:], dx[i])

    return im2

def replicate_offset(sp, dx):
    """
    :param sp: 1D spectrum to be shifted
    :param dx: offset in pixel space
    :return: matrix of spectra replicated at the different offsets
    """

    # initialize output
    out = np.zeros([len(dx),len(sp)])
    # loop over offsets
    for i in range(len(dx)):
        out[i,:] = shift(sp, dx[i])

    return out

def get_splined_template(sci_file, template_string = '/Volumes/courlan/decorr/templates/Template_s1dv_{}_sc1d_v_file_AB.fits'):
    """
    :param h: header of the template file
    :return: a spline of the template
    """

    # extract the header
    h = fits.getheader(sci_file)
    # get the object name
    objn = h['DRSOBJN']

    # get the wavelength solution
    wave = hdr2wave(h)
    # get the template file
    template_file = template_string.format(objn)

    # get the template
    with warnings.catch_warnings(record=True) as _:
        # we mask the warnings because the template files as sometimes
        # there's more than one extension
        template = Table.read(template_file)
    # retrieve the wavelength and flux
    wave_template = np.array(template['wavelength'])
    flux_tempalte = np.array(template['flux'])
    # valid pixels
    g = np.isfinite(flux_tempalte)
    # spline the template
    wave_doppler = doppler(wave, -h['BERV']*1000)
    spl = ius(wave_template[g], flux_tempalte[g], ext=1, k=3)
    # return the spline on the grid of the science file
    template2  = spl(wave_doppler)

    return template2

def get_params(instrument):
    """
    :param instrument: as an input
    :return: parameters for the file saving and header keywords
    """
    params = dict()

    if instrument.upper() == 'SPIROU':
        instrument = 'SPIRou'
    if instrument.upper() == 'NIRPS':
        instrument = 'NIRPS'

    params['instrument'] = instrument
    # width of the high-pass filter prior to saving the residuals
    params['Npca'] = 15
    # number of principal components to keep
    params['hpwidth'] = 100
    # we bin over 7 days and extend the determination of the sample by +-time_bin
    params['time_bin'] = 7.0
    # we cut out the outliers at 5 sigma (or any other value)
    params['nsig_cut'] = 5.0

    if getpass.getuser() == 'eartigau':
        if instrument == 'SPIRou':
            params['key_mjd'] = 'MJDATE'
            params['key_obj'] = 'DRSOBJN'
            params['template_string'] = '/Volumes/courlan/lbl_SPIROU/templates/Template_s1dv_{}_sc1d_v_file_AB.fits'
            params['residual_path'] = '/Volumes/courlan/decorr/residuals/'
            params['search_tcorr_path'] = '/Volumes/courlan/decorr/tcorr/*tcorr*fits'
        if instrument == 'NIRPS':
            params['key_mjd'] = 'MJD-OBS'
            params['key_obj'] = 'DRSOBJN'
            params['template_string'] = '/Volumes/courlan/lbl_NIRPS_HE/templates/Template_s1dv_{}_sc1d_v_file_A.fits'
            params['residual_path'] = '/Volumes/courlan/nirps_residuals/'
            params['search_tcorr_path'] = '/Volumes/courlan/decorr/tcorr/*tcorr*fits'

    if getpass.getuser() == 'spirou':
        if instrument == 'SPIRou':
            params['key_mjd'] = 'MJDATE'
            params['key_obj'] = 'DRSOBJN'

            params['template_string'] = \
                '/cosmos99/spirou/lbl-data/spirou_07282_online/templates/Template_s1dv_{}_sc1d_v_file_AB.fits'
            params['residual_path'] = '/space/spirou/LBL-PCA/residuals_spirou/'
            params['search_tcorr_path'] = '/space/spirou/LBL-PCA/SPIROU/science/*/*tcorr*'

        if instrument == 'NIRPS':
            params['key_mjd'] = 'MJD-OBS'
            params['key_obj'] = 'DRSOBJN'

            params['template_string'] = \
                '/cosmos99/nirps/lbl-data/nirps_he_online/templates/Template_s1dv_{}_sc1d_v_file_A.fits'
            params['residual_path'] = '/space/spirou/LBL-PCA/residuals_nirps/'
            params['search_tcorr_path'] = '/space/spirou/LBL-PCA/NIRPS_HE/science/*/*tcorr*'

    return params

def create_residuals(params):
    all_files = np.array(glob.glob(params['search_tcorr_path']))
    # check that the files are not PCA files
    keep = ['_PCA' not in all_files[i] for i in range(len(all_files))]
    all_files = all_files[keep]

    all_files = all_files[np.argsort(all_files)]
    mjds = np.zeros(len(all_files))
    objs = np.zeros(len(all_files), dtype = 'U99')
    template_exists = np.zeros(len(all_files), dtype = bool)
    tbl_summary = Table()


    errors = []
    for i in tqdm(range(len(all_files)), desc = 'reading through headers', leave = False):
        try:
            h = fits.getheader(all_files[i])
            mjds[i] = h[params['key_mjd']]
            objs[i] = h[params['key_obj']]
            template_exists[i] = os.path.exists(params['template_string'].format(objs[i]))
        except:
            errors.append('error reading header for {}\n'.format(all_files[i]))
    if len(errors) > 0:
        print('errors reading headers:')
        for ierr in errors:
            print(ierr)

    # only keep objects for which a template exists
    all_files = all_files[template_exists]
    mjds = mjds[template_exists]
    objs = objs[template_exists]
    outname = [params['residual_path']+f.split('/')[-1].replace('_tcorr_','_residual_') for f in all_files]

    tbl_summary['file'] = all_files
    tbl_summary['residual_file'] = outname
    tbl_summary[params['key_mjd']] = mjds
    tbl_summary[params['key_obj']] = objs

    uu = np.arange(0,1.1,0.1)
    u = np.unique(np.round(tbl_summary[params['key_mjd']] % 1, 1))
    midday =  np.nanmedian(uu[[i not in u for i in uu]])

    tbl_summary['day'] = np.round(tbl_summary[params['key_mjd']] - midday).astype(int)
    tbl_summary['unique_obj_day'] = np.zeros(len(tbl_summary), dtype = 'U99')
    for i in range(len(tbl_summary)):
        tbl_summary['unique_obj_day'][i] = tbl_summary[params['key_obj']][i] + '_' + tbl_summary['day'][i].astype(str)


    for i in tqdm(range(len(all_files)), desc = 'files', leave = False):
        if os.path.exists(outname[i]):
            continue
        h = fits.getheader(all_files[i])
        sp = np.array(fits.getdata(all_files[i]), dtype = float)
        template = get_splined_template(all_files[i], template_string = params['template_string'])
        g = np.isfinite(template)
        g &= np.isfinite(sp)
        g &= (template > 0)
        sp[g]/=template[g]
        sp[~g] = np.nan

        for iord in tqdm(range(sp.shape[0]), desc = 'order', leave = False):
            sp[iord,:] /= lowpassfilter(sp[iord,:], params['hpwidth'])

        with warnings.catch_warnings(record=True) as _:
            sp[sp<=0] = np.nan
            residual = np.log(sp)
        fits.writeto(outname[i], residual, h, overwrite = True)

    # now we create a table with the list of files and the corresponding output names
    files = tbl_summary['file']
    outnames = [f.split('/')[-1].replace('_tcorr_', '_PCA{}_'.format(str(params['Npca']).zfill(2))) for f in files]

    for i in tqdm(range(len(files)), desc = 'files', leave = False):
        path = '/'.join(files[i].split('/')[:-1]) + '_PCA{}'.format(str(params['Npca']).zfill(2)) + '/'
        if not os.path.exists(path):
            print('creating directory {}'.format(path))
            os.mkdir(path)
        outnames[i] = path + outnames[i]

    tbl_summary['pca_corrected_file'] = outnames


    # we also add a time bin column
    tbl_summary['TIMEBIN'] = np.array((tbl_summary[params['key_mjd']] - np.min(tbl_summary[params['key_mjd']])) // params['time_bin'], dtype=int)
    tbl_summary = tbl_summary[np.argsort(tbl_summary[params['key_mjd']])]



    tbl_summary.write('{}tbl_summary.csv'.format(params['residual_path']), overwrite = True)

    return tbl_summary

def pca_pix(instrument):
    # Which instrument, for now its only SPIRou or NIRPS
    #instrument = 'SPIRou'
    #instrument = 'NIRPS'

    # retrieve the parameters for the instrument
    params = get_params(instrument)
    # create the table with the list of files and the corresponding output names
    tbl = create_residuals(params)
    # read the first file to get the wavelength solution
    h = fits.getheader(tbl['file'][0])
    wavemap = hdr2wave(h)
    wave_middle = np.nanmedian(wavemap, axis = 1)
    # get the size of the wavelength map and overall size
    sz = wavemap.shape

    for utimebin in np.unique(tbl['TIMEBIN']):
        # loop on the time bins and create the list of files to be included in
        # the PCA
        sub_table = tbl[np.abs(tbl['TIMEBIN']- utimebin) <=1]
        all_files = sub_table['residual_file']

        # list of files that will be corrected
        files_of_interest = sub_table['TIMEBIN'] == utimebin

        # check if all files have been corrected. If so, skip this time bin
        pca_corrected_file = sub_table['pca_corrected_file'][files_of_interest]
        file_exists = np.array([os.path.exists(f) for f in pca_corrected_file])
        if False not in file_exists:
            print('all files exist for utimebin {}'.format(utimebin))
            continue

        Nobjs =  len(np.unique(sub_table['unique_obj_day']))

        # list the number of files that will be corrected
        print_outputs = (utimebin, len( np.unique(tbl['TIMEBIN'])), len(all_files),
              np.sum(files_of_interest),Nobjs)

        print('UtimeBin {}/{}, {} files in PCA, {} files corrected, {} Unique objs'.format(*print_outputs))

        # create the map of the residuals
        sample2 = np.zeros([len(all_files),np.product(sz)])
        # corresponding weights
        sigmap = np.zeros([len(all_files),np.product(sz)])


        # loop on the files and pad the residuals in a ravel() format
        for i in tqdm(range(len(all_files)), desc = 'files', leave = False):
            im = fits.getdata(all_files[i])
            rms = np.zeros(sz)
            with warnings.catch_warnings(record=True) as _:
                n1, med, p1 = np.nanpercentile(im, [16, 50, 84], axis=1)
            for iord in range(sz[0]):
                rms[iord,:] = (p1[iord] - n1[iord])/2
            rms[~np.isfinite(im)] = np.inf

            sigmap[i,:] = rms.ravel()
            sample2[i] = im.ravel()

        # weights inverse to the sigma, *not* the sigma**2
        with warnings.catch_warnings(record=True) as _:
            weights2 = 1/sigmap
            weights2[np.abs(sample2/sigmap)>8] = 0.0

        # remove the bad pixels from sample2 and weights2
        sample3 = np.array(sample2)

        nsig = sample3 / sigma(sample2[weights2 != 0])
        bad = ~np.isfinite(nsig)
        with warnings.catch_warnings(record=True) as _:
            # sigma-clipping at 5 or N sigma
            weights2[np.abs(nsig) > params['nsig_cut']] = 0.0
        # set spurious values to a weight
        weights2[bad] = 0.0
        sample2[bad] = 0.0

        # create the placeholder for the PCA model
        pca_model_all = np.zeros([np.sum(files_of_interest), np.product(sz)], dtype = float)


        sample2_utimebin = sample2[files_of_interest]
        weights2_utimebin = weights2[files_of_interest]

        uobj = np.unique(sub_table['unique_obj_day'])
        sample2_unique = np.zeros([len(uobj), np.product(sz)], dtype = float)
        weights2_unique = np.zeros([len(uobj), np.product(sz)], dtype = float)

        sample2[sample2 == 0] = np.nan

        for iuobj in tqdm(range(len(uobj)), desc = 'unique obj', leave = False):
            g = sub_table['unique_obj_day'] == uobj[iuobj]
            with warnings.catch_warnings(record=True) as _:# ignore the warnings about nanmedian
                sample2_unique[iuobj] = np.nanmedian(sample2[g], axis = 0)
                weights2_unique[iuobj] = np.nanmean(weights2[g], axis = 0)

        for iord in range(sz[0]):
            # loop on the orders and perform the PCA within each order
            # Thes eare files that are used to compute the PCA, utimebin -1 to
            sample2_ord = sample2_unique[:, iord*sz[1]:(iord+1)*sz[1]]
            weights2_ord = weights2_unique[:, iord*sz[1]:(iord+1)*sz[1]]
            # just the spectra that will be corrected. These are utimebin only
            sample2_utimebin_ord = sample2_utimebin[:, iord*sz[1]:(iord+1)*sz[1]]
            weights2_utimebin_ord = weights2_utimebin[:, iord*sz[1]:(iord+1)*sz[1]]

            # just in case somehting went wrong, remove the bad pixels

            with warnings.catch_warnings(record=True) as _:
                good = np.abs(sample2_ord*weights2_ord) < params['nsig_cut']
                good &= np.isfinite(sample2_ord)
                good &= np.isfinite(weights2_ord)
            weights2_ord[~good] = 0.0
            sample2_ord[~good] = 0.0

            # remove the orders that have too many bad pixels, 50% of the pixels
            frac_zero = np.nanmean(weights2_ord == 0,axis=1)
            g = frac_zero < 0.5
            if np.sum(g) <= 2*params['Npca']:
                continue

            sample2_ord = sample2_ord[g]
            weights2_ord = weights2_ord[g]

            g = np.nansum(weights2_ord != 0, axis=0) >= (params['Npca'] + 1)

            pca = EMPCA(n_components=params['Npca'])

            try:
                fit_pca = pca.fit(sample2_ord[:,g], weights=weights2_ord[:, g])
            except:
                print('PCA failed for order {}'.format(iord))
                continue

            # PCA reconstruction of the data for each line
            pca_model = np.zeros_like(sample2_utimebin_ord) + np.nan
            try:
                pca_model2 = fit_pca.reconstruct(sample2_utimebin_ord[:,g], weights=weights2_utimebin_ord[:, g])
            except:
                print('PCA reconstruction failed for order {}'.format(iord))
                continue

            pca_model[:, g] = pca_model2
            pca_model_all[:, iord*sz[1]:(iord+1)*sz[1]] = pca_model

            valid = weights2_utimebin_ord !=0
            valid &= (weights2_utimebin_ord > np.nanmedian(weights2_utimebin_ord[valid]))
            try:
                sig, sig2 = sigma(sample2_utimebin_ord[valid] - pca_model[valid]),sigma(sample2_utimebin_ord[valid])
            except:
                print('Error determination wrong for order {}'.format(iord))
                continue

            out_txt = 'order {:3.0f}/{:3.0f} - sigma before/after PCA = {:.4f} / {:.4f}, ratio {:.3f}, {:.3f}µm'
            print(out_txt.format(iord,sz[0],sig2, sig, sig/sig2,wave_middle[iord]))

            if False:
                w = hdr2wave(h)
                # plot on in 5 orders if doplot is True
                fig, ax = plt.subplots(3, 1, figsize=(10, 10), sharex=True, sharey=True)
                ax[0].imshow(sample2_utimebin_ord, aspect='auto', vmin=-2 * sig, vmax=2 * sig)
                ax[0].set(title='Initial\n{:.1f}nm - order {}'.format(np.nanmedian(w[iord]), iord), xlabel='Pixel',
                          ylabel='File')
                ax[1].imshow(pca_model, aspect='auto', vmin=-2 * sig, vmax=2 * sig)
                ax[1].set(title='PCA model', xlabel='Pixel', ylabel='File')
                ax[2].imshow(sample2_utimebin_ord - pca_model, aspect='auto', vmin=-2 * sig, vmax=2 * sig)
                ax[2].set(title='Corrected', xlabel='Pixel', ylabel='File')
                plt.tight_layout()
                plt.show()


        # we have the PCA model for each files, now we need to save the files by applying the correction
        # onto the original _tcorr_ files
        files = sub_table['file'][files_of_interest]
        outnames = sub_table['pca_corrected_file'][files_of_interest]

        for i in tqdm(range(len(sub_table[files_of_interest]   )), desc='Writing files'):
            # save file after correction with PCA
            corr = pca_model_all[i]
            corr.shape = sz
            with warnings.catch_warnings(record=True) as _:
                corr = np.exp(corr)
                # remove the pixels that too extreme
                corr[corr<0.01] = np.nan

            sp = fits.getdata(files[i])
            with fits.open(files[i]) as hdulist:
                if len(hdulist) == 1:
                    hdulist[0].data = sp / corr
                else:
                    hdulist[1].data = sp / corr

                hdulist.writeto(outnames[i], overwrite = True)

def get_residual_map(files, params):

    iord = params['iord']
    hpwidth = params['hpwidth']
    tplstr = params['template_string']

    residuals = np.zeros([len(files),4088])
    for i in tqdm(range(len(files)),desc='Computing residuals',leave=False):
        sp = fits.getdata(files[i])
        model = get_splined_template(files[i], template_string=tplstr)

        sp = sp[iord]
        model = model[iord]
        ratio = sp/model
        ratio/=lowpassfilter(ratio,hpwidth)
        residuals[i] = np.log(ratio)

    return residuals

def plot_residuals():

    params = get_params('NIRPS')
    params['iord'] = 26

    path1 = '/Volumes/courlan/lbl_SPIROU/science/GJ4274-PCA03/*.fits'
    files1 = np.sort(glob.glob(path1))
    #files1 = files1[:100]

    path2 = '/Volumes/courlan/lbl_SPIROU/science/GJ4274/*.fits'
    files2 = np.sort(glob.glob(path2))

    h = fits.getheader(files2[0])
    w = hdr2wave(h)

    n = np.min([len(files1),len(files2)])
    files1 = files1[:n]
    files2 = files2[:n]

    residuals1 = get_residual_map(files1, params)
    residuals2 = get_residual_map(files2, params)

    n1,p1 = np.nanpercentile(residuals1, [16, 84], axis=1)
    n2,p2 = np.nanpercentile(residuals2, [16, 84], axis=1)
    sig1 = (p1-n1)/2
    sig2 = (p2-n2)/2

    plt.plot(sig1, label = 'PCA')
    plt.plot(sig2, label = 'No PCA')
    plt.legend()
    plt.show()

    print(np.nanmedian(w[params['iord']]))
    fig, ax = plt.subplots(2,1,figsize=(10,5),sharex = True, sharey = True)
    ax[1].imshow(residuals1, aspect = 'auto', vmin = -0.05, vmax = 0.05)
    ax[1].set(title = 'PCA')
    ax[0].imshow(residuals2, aspect = 'auto', vmin = -0.05, vmax = 0.05)
    ax[0].set(title='No PCA')
    plt.show()
