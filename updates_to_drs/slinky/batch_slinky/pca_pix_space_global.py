import getpass
import glob
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline as ius
# do *not* use EMPCA as it has some random numbers inside and does not
# give the same results every time
from wpca import WPCA
from astropy.time import Time
import sys

from etienne_tools import hdr2wave, doppler, lowpassfilter, sigma, save_pickle,\
    read_pickle, printc, snail


def create_symlinks():
    path = '/cosmos99/nirps/apero-data/nirps_he_offline/red/*/*_pp_e2dsff_tcorr_A.fits'
    path_links = '/space/spirou/LBL-PCA/NIRPS_HE/science/{}/'
    files = glob.glob(path)

    bad_objs = ['SKY','MOON']

    for file in files:
        h = fits.getheader(file)
        if h['DRSOBJN'] in bad_objs:
            continue

        if not os.path.exists(path_links.format(h['DRSOBJN'])):
            print('creating directory {}'.format(path_links.format(h['DRSOBJN'])))
            os.mkdir(path_links.format(h['DRSOBJN']))
        cmd = 'ln -s {} {}'.format(file, path_links.format(h['DRSOBJN']))
        print(cmd)
        os.system(cmd)

    cmd = 'ln -s '

def push_to_maestria():
    cmd = 'scp -r /Users/eartigau/pycodes/etienne_hacks/' \
          'pixel_level_decorrelation/pca_pix_space_global.py ' \
          'spirou@maestria:/space/spirou/LBL-PCA'
    print(cmd)
    os.system(cmd)

def get_splined_template(sci_file,
                         template_string='/Volumes/courlan/decorr/templates/'
                                         'Template_s1dv_{}_sc1d_v_file_AB.fits'):
    """
    :param sci_file: science file
    :param template_string: template string
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
    wave_doppler = doppler(wave, -h['BERV'] * 1000)
    spl = ius(wave_template[g], flux_tempalte[g], ext=1, k=3)
    # return the spline on the grid of the science file
    template2 = spl(wave_doppler)

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
    params['Npca'] = 99
    params['Npca_adjust'] = True
    # number of principal components to keep
    params['hpwidth_kms'] = 50  # km/s
    # we bin over 3*timbe_bin days and extend the determination of the sample by +-time_bin
    params['time_bin'] = 15.0
    # we cut out the outliers at 5 sigma (or any other value)
    # params['nsig_cut'] = 5.0
    # do we plot
    params['doplot'] = True
    params['plot_orders'] = [0,5,10,15,20,30,63,66,68]
    params['hotstars'] = ['BET_PIC','ZETVIR','HD195094','HR875','HR4889','HR5671','HR6743','HR7590','HR8709','TAU3ERIDANI']
    params['ratio_rms'] = 1.33

    if getpass.getuser() == 'eartigau':
        if instrument == 'SPIRou':
            params['key_mjd'] = 'MJDATE'
            params['key_obj'] = 'DRSOBJN'
            params['template_string'] = '/Volumes/courlan/lbl_SPIROU/templates/Template_s1dv_{}_sc1d_v_file_AB.fits'
            params['residual_path'] = '/Volumes/courlan/decorr/residuals/'
            params['search_tcorr_path'] = '/Volumes/courlan/lbl_SPIROU/science/*/*tcorr*fits'
            params['pix_scale'] = 2.2  # km/s/pixel
            params['sample_blaze_file'] = '/Users/eartigau/.sample_blaze/sample_blaze_spirou.fits'

        if instrument == 'NIRPS':
            params['key_mjd'] = 'MJD-OBS'
            params['key_obj'] = 'DRSOBJN'
            params['template_string'] = '/Volumes/courlan/lbl_NIRPS_HE/templates/Template_s1dv_{}_sc1d_v_file_A.fits'
            params['residual_path'] = '/Volumes/courlan/decorr/residuals/'
            params['search_tcorr_path'] = '/Volumes/courlan/lbl_NIRPS_HE/science/*/*tcorr*fits'
            params['pix_scale'] = 0.95  # km/s/pixel
            params['sample_blaze_file'] = '/Users/eartigau/.sample_blaze/sample_blaze_nirps.fits'

    if getpass.getuser() == 'spirou':
        if instrument == 'SPIRou':
            params['key_mjd'] = 'MJDATE'
            params['key_obj'] = 'DRSOBJN'

            params['template_string'] = \
                '/cosmos99/spirou/apero-data/spirou_offline/lbl/templates/Template_s1dv_{}_sc1d_v_file_AB.fits'
            params['residual_path'] = '/space/spirou/LBL-PCA/residuals_spirou/'
            params['search_tcorr_path'] = '/space/spirou/LBL-PCA/SPIROU/science/*/*tcorr*'
            params['pix_scale'] = 2.2  # km/s/pixel
            params['sample_blaze_file'] = '/space/spirou/LBL-PCA/sample_blaze/sample_blaze_spirou.fits'

        if instrument == 'NIRPS':
            params['key_mjd'] = 'MJD-OBS'
            params['key_obj'] = 'DRSOBJN'

            params['template_string'] = \
                '/cosmos99/nirps/lbl-data/nirps_he_online/templates/Template_s1dv_{}_sc1d_v_file_A.fits'
            params['residual_path'] = '/space/spirou/LBL-PCA/residuals_nirps/'
            params['search_tcorr_path'] = '/space/spirou/LBL-PCA/NIRPS_HE/science/*/*tcorr*'
            params['pix_scale'] = 0.95  # km/s/pixel
            params['sample_blaze_file'] = '/space/spirou/LBL-PCA/sample_blaze/sample_blaze_nirps.fits'

    # to get an odd number
    params['hpwidth'] = (int(params['hpwidth_kms'] / params['pix_scale']) // 2) * 2 + 1


    tmp = fits.getdata(params['sample_blaze_file'])

    for i in range(tmp.shape[0]):
        tmp[i] /= np.nanmedian(tmp[i])

    params['sample_blaze'] = tmp


    return params


def create_residuals(params):
    all_files = np.array(glob.glob(params['search_tcorr_path']))
    all_files = all_files[np.argsort(np.random.rand(len(all_files)))]

    # check that the files are not PCA files
    keep = ['_PCA' not in all_files[i] for i in range(len(all_files))]
    all_files = all_files[keep]

    # all_files = all_files[np.argsort(all_files)]
    mjds = np.zeros(len(all_files))
    objs = np.zeros(len(all_files), dtype='U99')
    template_exists = np.zeros(len(all_files), dtype=bool)
    tbl_summary = Table()

    errors = []
    for i in snail(range(len(all_files)), desc='reading through headers', leave=False):
        try:
            temp_hdr = params['residual_path'] + all_files[i].split('/')[-1].replace('.fits', '_hdr.pickle')
            if not os.path.exists(temp_hdr):
                h = dict(fits.getheader(all_files[i]))
                save_pickle(temp_hdr, h)
            else:
                h = read_pickle(temp_hdr)

            mjds[i] = h[params['key_mjd']]
            objs[i] = h[params['key_obj']]
            template_exists[i] = os.path.exists(params['template_string'].format(objs[i]))
        except:
            # exit code if this is a keyboard interrupt
            if sys.exc_info()[0] == KeyboardInterrupt:
                printc('We have a keywoard interrupt', 'bad2')
                raise KeyboardInterrupt

            errors.append('error reading header for {}\n'.format(all_files[i]))
    if len(errors) > 0:
        printc('errors reading headers:', 'bad3')
        for ierr in errors:
            printc(ierr,'bad3')

    # only keep objects for which a template exists
    all_files = all_files[template_exists]
    mjds = mjds[template_exists]
    objs = objs[template_exists]
    outname = [params['residual_path'] + f.split('/')[-1].replace('_tcorr_', '_residual_') for f in all_files]
    outname = [outname[i].replace('.fits', '.pick') for i in range(len(outname))]

    tbl_summary['file'] = all_files
    tbl_summary['residual_file'] = outname
    tbl_summary[params['key_mjd']] = mjds
    tbl_summary[params['key_obj']] = objs
    tbl_summary['hotstar'] = [obj in params['hotstars'] for obj in objs]

    uu = np.arange(0, 1.1, 0.1)
    u = np.unique(np.round(tbl_summary[params['key_mjd']] % 1, 1))
    midday = np.nanmedian(uu[[i not in u for i in uu]])

    tbl_summary['day'] = np.round(tbl_summary[params['key_mjd']] - midday).astype(int)
    tbl_summary['SNR'] = np.zeros(len(tbl_summary), dtype=float)
    tbl_summary['SHAPE_DX'] = np.zeros(len(tbl_summary), dtype=float)
    tbl_summary['SHAPE_DY'] = np.zeros(len(tbl_summary), dtype=float)
    tbl_summary['unique_obj_day'] = np.zeros(len(tbl_summary), dtype='U99')

    for i in range(len(tbl_summary)):
        str1 = tbl_summary[params['key_obj']][i]
        str2 = tbl_summary['day'][i].astype(str)
        uobj = '{}_{}'.format(str1, str2)
        tbl_summary['unique_obj_day'][i] = uobj

    for i in snail(range(len(tbl_summary)), desc='files', leave=False):
        if os.path.exists(tbl_summary['residual_file'][i]):
            tmp = read_pickle(tbl_summary['residual_file'][i])
            h = tmp['header']
            tbl_summary['SNR'][i] = h['SNR']
            tbl_summary['SHAPE_DX'][i] = h['SHAPE_DX']
            tbl_summary['SHAPE_DY'][i] = h['SHAPE_DY']
            continue
        h = fits.getheader(tbl_summary['file'][i])
        tbl_summary['SHAPE_DX'][i] = h['SHAPE_DX']
        tbl_summary['SHAPE_DY'][i] = h['SHAPE_DY']

        sp = np.array(fits.getdata(tbl_summary['file'][i]), dtype=float)

        try:
            if tbl_summary['hotstar'][i] == False:
                template = get_splined_template(all_files[i], template_string=params['template_string'])
            else:
                template = np.ones_like(sp)
        except:
            errors.append('error reading header for {}\n'.format(all_files[i]))
            continue

        g = np.isfinite(template)
        g &= np.isfinite(sp)
        g &= (template > 0)
        sp[g] /= template[g]
        sp[~g] = np.nan

        for iord in snail(range(sp.shape[0]), desc='order', leave=False):
            with warnings.catch_warnings(record=True) as _:
                sp[iord, :] /= lowpassfilter(sp[iord, :], params['hpwidth'])

        with warnings.catch_warnings(record=True) as _:
            sp[sp <= 0] = np.nan
            residual = np.log(sp)

        with warnings.catch_warnings(record=True) as _:
            h['SNR'] = 1 / sigma(sp.ravel() - np.roll(sp, 1, axis=1).ravel())
        tbl_summary['SNR'][i] = h['SNR']

        tmp = dict()
        tmp['data'] = residual
        tmp['header'] = dict(h)
        save_pickle(outname[i], tmp)

    files = tbl_summary['file']
    if params['Npca_adjust']:
        # now we create a table with the list of files and the corresponding output names
        outnames = [f.split('/')[-1].replace('_tcorr_', '_PCAx_') for f in files]
    else:
        outnames = [f.split('/')[-1].replace('_tcorr_', '_PCA{}_'.format(str(params['Npca']).zfill(2))) for f in files]

    for i in snail(range(len(files)), desc='files', leave=False):
        if params['Npca_adjust']:
            path = '/'.join(files[i].split('/')[:-1]) + '_PCAx' + '/'
        else:
            path = '/'.join(files[i].split('/')[:-1]) + '_PCA{}'.format(str(params['Npca']).zfill(2)) + '/'
        if not os.path.exists(path):
            printc('creating directory {}'.format(path),'info')
            os.mkdir(path)
        outnames[i] = path + outnames[i]

    tbl_summary['pca_corrected_file'] = outnames

    # we also add a time bin column
    tbl_summary['TIMEBIN'] = np.array(
        (tbl_summary[params['key_mjd']] - np.min(tbl_summary[params['key_mjd']])) // params['time_bin'], dtype=int)
    tbl_summary = tbl_summary[np.argsort(tbl_summary[params['key_mjd']])]

    tbl_summary.write('{}tbl_summary.csv'.format(params['residual_path']), overwrite=True)

    return tbl_summary


def pca_pix(instrument, doplot = False):
    # Which instrument, for now its only SPIRou or NIRPS
    # instrument = 'SPIRou'
    # instrument = 'NIRPS'

    # retrieve the parameters for the instrument
    params = get_params(instrument)

    # to plot or not to plot, that is the question
    params['doplot'] = doplot

    # create the table with the list of files and the corresponding output names
    # this ensures that have files have a residual vector
    tbl = create_residuals(params)

    # read the first file to get the wavelength solution
    h = fits.getheader(tbl['file'][0])
    # this is only for the title of the plots and printout for user,
    # it has no impact on the code
    wavemap = hdr2wave(h)
    wave_middle = np.nanmedian(wavemap, axis=1)
    # get the size of the wavelength map and overall size
    sz = wavemap.shape

    # check that all files have the corresponding residual file
    keep = np.zeros_like(tbl, dtype=bool)
    for i in range(len(tbl)):
        keep[i] = os.path.exists(tbl['residual_file'][i])

    # print the number of files that will be corrected and the number of objects
    # that were passed initially. Normally these two are the same.
    printc('Keeping {} files out of {}'.format(np.sum(keep), len(tbl)), 'info')
    tbl = tbl[keep]

    # loop on the time bins and create the list of files to be included in
    # the PCA. The PCA is computed on files +-1 time bin from the time bin
    # considered at each moment
    utimebins = np.unique(tbl['TIMEBIN'])
    for itimebin in range(len(utimebins)):
        #
        printc('~')

        utimebin = utimebins[itimebin]
        # loop on the time bins and create the list of files to be included in
        # the PCA
        sub_table = tbl[np.abs(tbl['TIMEBIN'] - utimebin) <= 1]
        all_files = sub_table['residual_file']

        if params['Npca_adjust']:
            params['Npca'] = len(np.unique(sub_table['day'])) // 2
            if params['Npca'] < 2:
                params['Npca'] = 2

        printc('Npca = {}'.format(params['Npca']),'info')

        # list of files that will be corrected
        files_of_interest = sub_table['TIMEBIN'] == utimebin

        # check if all files have been corrected. If so, skip this time bin
        pca_corrected_file = sub_table['pca_corrected_file'][files_of_interest]
        file_exists = np.array([os.path.exists(f) for f in pca_corrected_file])

        if False not in file_exists:
            # all files exist, skip this time bin
            msg = (utimebin, itimebin, len(utimebins))
            printc('all files exist for utimebin {} [{}/{}]'.format(*msg), 'bad1')
            continue

        # Number of unique combination of objects * dates
        Nobjs = len(np.unique(sub_table['unique_obj_day']))

        # list the number of files that will be corrected
        print_outputs = (itimebin, len(utimebins), len(all_files),
                         np.sum(files_of_interest), Nobjs, len(np.unique(sub_table['day'])))

        printc('UtimeBin {}/{}, {} files in PCA, {} files corrected, {} Unique objs, {} Unique nights'.format(
            *print_outputs),'number')

        # for the printout, list the first and last day of the files that will
        # be corrected and the first and last day of the files that will be
        # used to compute the PCA
        uday1 = np.unique(sub_table['day'])
        d1 = Time(uday1[0],format='mjd').iso.split(' ')[0]
        d2 = Time(uday1[-1],format='mjd').iso.split(' ')[0]
        printc('First day {} -> last day {} of PCA files'.format(d1,d2),'number')
        uday2 = np.unique(sub_table['day'][files_of_interest])
        d1 = Time(uday2[0],format='mjd').iso.split(' ')[0]
        d2 = Time(uday2[-1],format='mjd').iso.split(' ')[0]
        printc('First day {} -> last day {} of corrected files'.format(d1,d2),'number')

        # create the map of the residuals. This has the size of the (order * 4088) x the number of files
        # in the time bin. Sample2 is a ravel() version of the spectrum.
        sample2 = np.zeros([len(all_files), np.product(sz)])
        # corresponding weights
        sigmap = np.zeros([len(all_files), np.product(sz)])

        # loop on the files and pad the residuals in a ravel() format
        for i in snail(range(len(all_files)), desc='files', leave=False):
            tmp = read_pickle(all_files[i])
            im = tmp['data']*params['sample_blaze']
            rms = np.zeros(sz)
            with warnings.catch_warnings(record=True) as _:
                n1, med, p1 = np.nanpercentile(im - np.roll(im,1,axis=1), [16, 50, 84], axis=1)
            for iord in range(sz[0]):
                rms[iord, :] = (p1[iord] - n1[iord]) / 2
            rms[~np.isfinite(im)] = np.inf

            sigmap[i, :] = rms.ravel()
            sample2[i] = im.ravel()

        # weights inverse to the sigma, *not* the sigma**2
        with warnings.catch_warnings(record=True) as _:
            weights2 = 1 / sigmap

        # remove the bad pixels from sample2 and weights2
        sample3 = np.array(sample2)

        with warnings.catch_warnings(record=True) as _:
            # any point with a residual value of >0.5 (remembre, this is in log space)
            # is really bad. This corresponds to an error of ~30%
            bad = ~((np.abs(sample3) < 0.5) *(np.isfinite(weights2)))


        weights2[bad] = 0.0
        sample2[bad] = 0.0

        # create the placeholder for the PCA model
        pca_model_all = np.zeros([np.sum(files_of_interest), np.product(sz)], dtype=float)

        sample2_utimebin = sample2[files_of_interest]
        weights2_utimebin = weights2[files_of_interest]

        uobj = np.unique(sub_table['unique_obj_day'])
        sample2_unique = np.zeros([len(uobj), np.product(sz)], dtype=float)

        sig_unique = np.zeros(len(uobj))
        for iuobj in snail(range(len(uobj)), desc='unique obj', leave=False):
            g = sub_table['unique_obj_day'] == uobj[iuobj]
            with warnings.catch_warnings(record=True) as _:  # ignore the warnings about nanmedian
                sample2_unique[iuobj] = np.nanmedian(sample2[g], axis=0)

        #for iord in range(sz[0]):
        params['Npca'] = (len(np.unique(sub_table['day'])) // 2)
        for iord in range(sz[0]):
            # loop on the orders and perform the PCA within each order
            # Thes eare files that are used to compute the PCA, utimebin -1 to
            sample2_ord = sample2_unique[:, iord * sz[1]:(iord + 1) * sz[1]]
            weights2_ord = np.zeros_like(sample2_ord)

            for i in snail(range(sample2_ord.shape[0]), desc='unique obj', leave=False):
                diff = sample2_ord[i] - np.roll(sample2_ord[i], 1)
                bad = sample2_ord[i] == 0
                diff[bad] = np.nan
                diff = diff[(diff != 0) * (np.isfinite(diff))]
                if len(diff) < sample2_ord.shape[1] // 4:
                    weights2_ord[i] = 0.0
                    continue
                sig = np.sqrt( sigma(diff)**2 + 0.005**2 )
                #plt.plot(i,sig,'go')

                weights2_ord[i] = 1 / sig

            bad = ~np.isfinite(sample2_ord+weights2_ord)
            weights2_ord[bad] = 0.0
            sample2_ord[bad] = 0.0

            # just the spectra that will be corrected. These are utimebin only
            sample2_utimebin_ord = sample2_utimebin[:, iord * sz[1]:(iord + 1) * sz[1]]
            weights2_utimebin_ord = weights2_utimebin[:, iord * sz[1]:(iord + 1) * sz[1]]

            # just in case somehting went wrong, remove the bad pixels


            """
            for i in range(sample2_ord.shape[0]):
                v = sample2_ord[i]*weights2_ord[i]
                v = v[(v!=0)*(np.isfinite(v))]
                try:
                    n2,p2 = np.nanpercentile(v,[16,85])
                    sig2 = (p2-n2)/2
                    invalid = np.abs(sample2_ord[i]*weights2_ord[i]) > 8*sig2
                    v = sample2_ord[i]
                    v[invalid] = np.nan
                    sample2_ord[i] = v
                except:
                    pass
            """

            with warnings.catch_warnings(record=True) as _:
                # good = np.abs(sample2_ord*weights2_ord) < params['nsig_cut']
                good = np.isfinite(sample2_ord)
                good &= np.isfinite(weights2_ord)

            weights2_ord[~good] = 0.0
            sample2_ord[~good] = 0.0
            with warnings.catch_warnings(record=True) as _:
                weights2_ord/=np.nanmean(weights2_ord[weights2_ord!=0])

            # remove the files that have too many bad pixels, 50% of the pixels
            # must be valid
            frac_zero = np.nansum(weights2_ord == 0, axis=0)
            g = frac_zero < 2 * params['Npca']

            with warnings.catch_warnings(record=True) as _:
                pca = WPCA(n_components=(weights2_ord.shape[0]//2-1))

            try:
                with warnings.catch_warnings(record=True) as _:
                    #fit_pca = pca.fit(sample2_ord[:, g], weights=weights2_ord[:, g])
                    #npca = np.sum(fit_pca.explained_variance_ratio_>params[
                    #    'ratio_rms']*fit_pca.explained_variance_ratio_[-1])

                    #if npca> nmaxpca:
                    #    npca =nmaxpca
                    #printc('order {:3.0f}/{:3.0f} - npca = {}'.format(iord, sz[0], npca),'number')
                    pca = WPCA(n_components=params['Npca'])
                    fit_pca = pca.fit(sample2_ord[:, g], weights=weights2_ord[:, g])


            except:
                # exit code if this is a keyboard interrupt
                if sys.exc_info()[0] == KeyboardInterrupt:
                    printc('We have a keywoard interrupt', 'bad2')

                    #return

                printc('PCA failed for order {}'.format(iord),'bad2')

                continue

            # PCA reconstruction of the data for each line
            pca_model = np.zeros_like(sample2_utimebin_ord)  # + np.nan
            weights_tmp = weights2_utimebin_ord[:, g]
            try:
                with warnings.catch_warnings(record=True) as _:
                    pca_model2 = fit_pca.reconstruct(sample2_utimebin_ord[:, g], weights=weights_tmp)

                residual = sample2_utimebin_ord[:, g] - pca_model2
                n1, med, p1 = np.nanpercentile(residual, [16, 50, 84], axis=1)
                sig = (p1 - n1) / 2

                for i in range(residual.shape[0]):
                    nsig = ((residual[i] - med[i]) / sig[i])
                    with warnings.catch_warnings(record=True) as _:
                        bad = np.abs(nsig) > 8
                    weights_tmp[i][bad] = 0.0
                with warnings.catch_warnings(record=True) as _:
                    pca_model2 = fit_pca.reconstruct(sample2_utimebin_ord[:, g], weights=weights_tmp)

            except:
                # exit code if this is a keyboard interrupt
                if sys.exc_info()[0] == KeyboardInterrupt:
                    printc('We have a keywoard interrupt', 'bad2')

                    #return

                printc('PCA reconstruction failed for order {}'.format(iord),'bad2')
                continue

            pca_model2[weights2_utimebin_ord[:, g] == 0] = 0.0
            pca_model[:, g] = pca_model2
            pca_model_all[:, iord * sz[1]:(iord + 1) * sz[1]] = pca_model

            valid = weights2_utimebin_ord != 0
            valid &= sample2_utimebin_ord != 0
            # mask &= (weights2_utimebin_ord > np.nanmedian(weights2_utimebin_ord[mask]))
            nanmask = np.zeros_like(valid, dtype=float)
            nanmask[~valid] = np.nan

            n1, p1 = np.nanpercentile(sample2_utimebin_ord + nanmask, [16, 84], axis=1)
            sig1 = (p1 - n1) / 2

            n1_after, p1_after = np.nanpercentile(sample2_utimebin_ord - pca_model + nanmask, [16, 84], axis=1)
            sig1_after = (p1_after - n1_after) / 2

            gain = np.nanmedian(sig1_after / sig1)

            out_txt = 'order {:3.0f}/{:3.0f} - sigma before/after PCA = {:.4f} / {:.4f}, med(gain) {:.3f}, ' \
                      ' NPCA= {} {:.3f}µm'
            printc(out_txt.format(iord, sz[0], np.mean(sig1), np.mean(sig1_after), gain, params["Npca"],
                                  wave_middle[ iord]),'number')


            nanmask = np.zeros_like(sample2_utimebin_ord, dtype=float)
            nanmask[~np.isfinite(sample2_utimebin_ord)] = np.nan
            nanmask[pca_model == 0] = np.nan
            nanmask[sample2_utimebin_ord == 0] = np.nan
            #sig = np.nanpercentile(sig1, 75)

            if params['doplot'] and (iord in params['plot_orders']):

                valid = np.nanmean(np.isfinite(nanmask), axis=0) > 0.5
                sig = 0.05
                rms_map = np.ones_like(sample2_utimebin_ord)
                for i in range(rms_map.shape[0]):
                    rms_map[i] *= sig1[i]

                xmin = np.nanmin(np.where(valid)[0])
                xmax = np.nanmax(np.where(valid)[0])


                w = hdr2wave(h)
                # plot on in 5 orders if doplot is True
                fig, ax = plt.subplots(3, 1, figsize=(10, 10), sharex='all', sharey='all')
                ax[0].imshow(sample2_utimebin_ord, aspect='auto', vmin=-3 * sig, vmax=3 * sig)
                ax[0].set(title='Initial\n{:.1f}nm - order {}'.format(np.nanmedian(w[iord]), iord), xlabel='Pixel',
                          ylabel='File')
                ax[1].imshow(pca_model + nanmask, aspect='auto', vmin=-3 * sig, vmax=3 * sig)
                ax[1].set(title='PCA model', xlabel='Pixel', ylabel='File')
                ax[2].imshow((sample2_utimebin_ord-pca_model), aspect='auto', vmin=-3 * sig, \
                    vmax=3 * sig)
                ax[2].set(title='Corrected', xlabel='Pixel', ylabel='File')

                ax[2].set(xlim=[xmin, xmax])

                plt.tight_layout()
                plt.show()

        # we have the PCA model for each file, now we need to save the files by applying the correction
        # onto the original _tcorr_ files
        files = sub_table['file'][files_of_interest]
        outnames = sub_table['pca_corrected_file'][files_of_interest]

        printc('Writing files [n={}]'.format(len(sub_table)), 'info')
        for i in snail(range(len(sub_table[files_of_interest])), desc='Writing files', leave=False):
            # save file after correction with PCA
            corr = pca_model_all[i]/params['sample_blaze'].ravel()
            corr.shape = sz
            with warnings.catch_warnings(record=True) as _:
                corr = np.exp(corr)
                # remove the pixels that too extreme
                corr[corr < 0.01] = np.nan

            sp = fits.getdata(files[i])
            try:
                with fits.open(files[i]) as hdulist:
                    if len(hdulist) == 1:
                        hdulist[0].data = sp / corr
                    else:
                        hdulist[1].data = sp / corr

                    hdulist.writeto(outnames[i], overwrite=True)
            except:
                printc('Error writing {}'.format(outnames[i]),'bad2')
                continue

def get_residual_map(files, params):
    iord = params['iord']
    hpwidth = params['hpwidth']
    tplstr = params['template_string']

    residuals = np.zeros([len(files), 4088])
    for i in snail(range(len(files)), desc='Computing residuals', leave=False):
        sp = fits.getdata(files[i])
        model = get_splined_template(files[i], template_string=tplstr)

        sp = sp[iord]
        model = model[iord]
        ratio = sp / model
        ratio /= lowpassfilter(ratio, hpwidth)
        residuals[i] = np.log(ratio)

    return residuals


def plot_residuals():
    params = get_params('NIRPS')
    params['iord'] = 67

    obj = 'TOI1078'

    path1 = '/Volumes/courlan/lbl_NIRPS_HE/science/{}_PCAx/*.fits'.format(obj)
    files1 = np.sort(glob.glob(path1))
    # files1 = files1[:100]

    path2 = '/Volumes/courlan/lbl_NIRPS_HE/science/{}/*.fits'.format(obj)
    files2 = np.sort(glob.glob(path2))

    h = fits.getheader(files2[0])
    w = hdr2wave(h)

    n = np.min([len(files1), len(files2)])
    files1 = files1[:n]
    files2 = files2[:n]
    date = [fits.getheader(f)['DATE'].split('T')[0] for f in files1]
    residuals1 = get_residual_map(files1, params)
    residuals2 = get_residual_map(files2, params)

    n1, p1 = np.nanpercentile(residuals1, [16, 84], axis=1)
    n2, p2 = np.nanpercentile(residuals2, [16, 84], axis=1)
    sig1 = (p1 - n1) / 2
    sig2 = (p2 - n2) / 2

    plt.plot(sig1, label='PCA')
    plt.plot(sig2, label='No PCA')
    plt.legend()
    plt.show()

    #printc(np.nanmedian(w[params['iord']]),'info')
    fig, ax = plt.subplots(2, 1, figsize=(10, 5), sharex='all', sharey='all')
    ax[1].imshow(residuals1, aspect='auto', vmin=-0.05, vmax=0.05)
    ax[1].set(title='PCA')
    ax[0].imshow(residuals2, aspect='auto', vmin=-0.05, vmax=0.05)
    ax[0].set(title='No PCA')

    for i in range(len(date)-1):
        if date[i] == date[i+1]:
            continue

        ax[0].hlines(i+.5,0, 4088, color='w', alpha=0.5)
        ax[1].hlines(i+.5,0, 4088, color='w', alpha=0.5)

    plt.show()



#pca_pix('NIRPS', doplot=False)