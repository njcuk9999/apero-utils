import getpass
import glob
import os
import sys
import warnings
from tqdm import tqdm

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from sklearn.cluster import Birch
# do *not* use EMPCA as it has some random numbers inside and does not
# give the same results every time
from wpca import WPCA, EMPCA

from astromdatabase import fetch_astrom
from pixtools import doppler, lowpassfilter, sigma, save_pickle, \
    read_pickle, printc, snail, dict2mef, mef2dict, load_yaml, read_t, write_t

import datetime

from concurrent.futures import ThreadPoolExecutor, as_completed
import concurrent.futures

def sigma(tmp):
    return np.nanpercentile(np.abs(tmp - np.nanmedian(tmp)), 68) / 0.6745

def get_med_residual(res):
    for ite in range(3):
        med = np.nanmedian(res, axis=0)
        for i in range(res.shape[0]):
            g = np.isfinite(res[i]*med)
            amps = np.nansum(res[i][g]*med[g])/np.nansum(med[g]**2)
            res[i]/=amps

    med = np.nanmedian(res, axis=0)
    return med

def get_splined_template(sci_file, 
                         template_string = None, params = None):
    """
    :param sci_file: science file
    :param template_string: template string
    :return: a spline of the template
    """
    if 'dict_splines' not in params:
        params['dict_splines'] = dict()
        # This dictionnary has one entry per object and then one spline per order

    if params is None:
        raise ValueError('params is None')

    if template_string is None:
        template_string = params['template_string']

    # extract the data
    data_dict = read_t(sci_file)

    fiber_setup = params['fibresetup']

    # get the object name
    objn = data_dict[f'Flux{fiber_setup}_header']['DRSOBJN']

    # get the wavelength solution
    wave = data_dict['Wave'+fiber_setup]
    sp = data_dict['FluxA']

    if objn not in params['dict_splines']:
        printc(f'We add {objn} to the spline dict','info')
        params['dict_splines'][objn] = dict()

        # get the master wavelength solution
        wave_template = fits.getdata(params['waveref_file_name'])
        if data_dict['FluxA_header']['DRSOBJN'] not in params['hotstars']:

            spectrum_template = fits.getdata(template_string.format(objn))
            
            for iord in range(sp.shape[0]):
                wave_ord = wave_template[iord]
                params['dict_splines'][objn][iord] = dict()
                spectrum_template_ord = spectrum_template[iord]
                g = np.isfinite(spectrum_template_ord)
                if np.sum(g)>10:
                    params['dict_splines'][objn][iord] = ius(wave_ord[g],spectrum_template_ord[g], ext=1, k=2)
                else:
                    params['dict_splines'][objn][iord] = ius(wave_ord,np.zeros_like(wave_ord), ext=1, k=2)

        else:
            for iord in range(sp.shape[0]):
                wave_ord = wave_template[iord]
                params['dict_splines'][objn][iord] = dict()
                params['dict_splines'][objn][iord] = ius(wave_ord,np.ones_like(wave_ord), ext=1, k=1)



    template2 = np.zeros_like(wave)

    # spline the template
    # return the spline on the grid of the science file
    wave_doppler = doppler(wave, -data_dict['FluxA_header']['BERV'] * 1000)
    for iord in range(template2.shape[0]):
        template2[iord] = params['dict_splines'][objn][iord](wave_doppler[iord])

    template2[template2 == 0] = np.nan


    # Apply masks to the template and science data
    g = np.isfinite(template2)
    g &= np.isfinite(sp)
    g &= (template2 > 0)

    bad = ~g
    sp[bad] = np.nan


    # run the per-order processing in parallel
    def process_order(iord):
        ratio = lowpassfilter(sp[iord] / template2[iord], params['hpwidth'])
        template2[iord] *= ratio

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        executor.map(process_order, snail(range(sp.shape[0]), desc='order', leave=False))


    return template2, params


def get_params(yaml_file=None):
    """
    :param instrument: as an input
    :return: parameters for the file saving and header keywords
    """

    if getpass.getuser() == 'eartigau':
        yaml_path = '/Users/eartigau/pycodes/pixeldecorr/yamls/'
    if getpass.getuser() == 'spirou':
        yaml_path = '/space/spirou/LBL-PCA/yamls/'

    if yaml_file is None:
        yaml_files = glob.glob(yaml_path + '*.yaml')
        printc('Choose a yaml file:', 'number')
        for i, yaml_file in enumerate(yaml_files):
            printc('[{:3.0f}] {}'.format(i, yaml_file), 'number')
        i = int(input('Enter the number of the yaml file: '))
        yaml_file = yaml_files[i]

    params = load_yaml(yaml_file)

    if params['instrument'].upper() == 'SPIROU':
        params['instrument'] = 'SPIRou'
    if params['instrument'].upper() == 'NIRPS':
        params['instrument'] = 'NIRPS'

    params['err_file'] = yaml_file.replace('.yaml', '.err')

    if 'Npca' not in params.keys():
        params['Npca'] = 99
    if 'Npca_adjust' not in params.keys():
        params['Npca_adjust'] = True
    if 'npca-max' not in params.keys():
        # number of principal components to keep
        params['npca-max'] = 15

    if 'hpwidth_kms' not in params.keys():
        # width of the high-pass filter prior to saving the residuals
        params['hpwidth_kms'] = 50  # km/s
    if 'time_bin' not in params.keys():
        # we bin over 3*timbe_bin days and extend the determination of the sample by +-time_bin
        params['time_bin'] = 15.0

    if 'doplot' not in params.keys():
        params['doplot'] = False

    if 'plot_orders' not in params.keys():
        params['plot_orders'] = [0, 5, 10, 15, 20, 30, 63, 66, 68]

    astrometric_db = fetch_astrom()
    params['hotstars'] = np.array(astrometric_db['OBJNAME'][astrometric_db['TELLURIC'] == True])

    if 'object_of_interest' not in params.keys():
        params['object_of_interest'] = ['GJ707', 'GJ3737', 'GL406', 'TOI406']

    if 'ratio_rms' not in params.keys():
        params['ratio_rms'] = 1.33

    if getpass.getuser() == 'eartigau':
        if params['instrument'] == 'SPIRou':
            params['key_mjd'] = 'MJDATE'
            params['key_obj'] = 'DRSOBJN'
            params['template_string'] = '/Volumes/courlan/lbl_SPIROU/templates/Template_s1dv_{}_sc1d_v_file_AB.fits'
            params['residual_path'] = '/Volumes/courlan/decorr/residuals/'
            #params['search_t_slinky_path'] = '/Volumes/courlan/SLINKY/data_SPIROU_output/*_slinky/*t_slinky.fits'
            params['pix_scale'] = 2.2  # km/s/pixel
            params['pca_mef_dir'] = '/Volumes/courlan/decorr/pca_mef/'
            params['yaml_folder'] = '/Users/eartigau/pycodes/pixeldecorr/yamls/'

        if params['instrument'] == 'NIRPS_HE':
            params['key_mjd'] = 'MJD-OBS'
            params['key_obj'] = 'DRSOBJN'
            params['template_string'] = '/Volumes/courlan/lbl_NIRPS_HE/templates/Template_{}_sc1d_v_file_A.fits'
            params['residual_path'] = '/Volumes/courlan/decorr/residuals/'
            #params['search_t_slinky_path'] = '/Volumes/courlan/lbl_NIRPS_HE/science/*/*t.fits'
            params['pix_scale'] = 0.95  # km/s/pixel
            params['pca_mef_dir'] = '/Volumes/courlan/decorr/pca_mef/'

    if getpass.getuser() == 'spirou':
        if params['instrument'] == 'SPIRou':
            params['key_mjd'] = 'MJDATE'
            params['key_obj'] = 'DRSOBJN'

            params['template_string'] = \
                '/cosmos99/spirou/apero-data/spirou_offline/red/other/Template_{}_tellu_obj_AB.fits'
            params['residual_path'] = '/space/spirou/LBL-PCA/residuals_spirou/'
            #params['search_t_slinky_path'] = ' /space/spirou/SLINKY/data_SPIROU_output/*_slinky/*t_slinky.fits'
            params['pix_scale'] = 2.2  # km/s/pixel
            params['pca_mef_dir'] = '/space/spirou/LBL-PCA/residuals_spirou/pca_mef/'
            params['patched_wavesol'] = '/space/spirou/SLINKY/calib_SPIROU_updatedwavesol/'
            params['calib_dir'] = '/space/spirou/SLINKY/calib_SPIROU/'
            params['fibresetup'] = 'AB'
            params['sample_order'] = 35

        if params['instrument'] == 'NIRPS_HE':
            params['key_mjd'] = 'MJD-OBS'
            params['key_obj'] = 'DRSOBJN'

            params['template_string'] = \
                '/cosmos99/nirps/apero-data/nirps_he_online/red/other/Template_{}_tellu_obj_A.fits'
            params['residual_path'] = '/space/spirou/LBL-PCA/residuals_nirps/'
            #
            params['pix_scale'] = 0.95  # km/s/pixel
            params['pca_mef_dir'] = '/space/spirou/LBL-PCA/residuals_nirps/pca_mef/'
            params['patched_wavesol'] = '/space/spirou/SLINKY/calib_NIRPS_HE_updatedwavesol/'
            params['calib_dir'] = '/space/spirou/SLINKY/calib_NIRPS_HE/'
            params['fibresetup'] = 'A'
            params['sample_order'] = 65


            params['waveref_file_name'] = glob.glob('/cosmos99/nirps/apero-data/nirps_he_online/calib/*_pp_e2dsff_A_wavesol_ref_A.fits')[0]



    if 'MAX_CLUSTERS' not in params.keys():
        params['MAX_CLUSTERS'] = 30  # max number of clusters in the Birch algorithm with SHAPE_DX and SHAPE_DY

    if 'CUT_MAX_SPAN_CLUSTER' not in params.keys():
        params['CUT_MAX_SPAN_CLUSTER'] = 0.2  # pixels

    if 'hpwidth' not in params.keys():
        # to get an odd number
        params['hpwidth'] = (int(params['hpwidth_kms'] / params['pix_scale']) // 2) * 2 + 1

    if not params['pca_mef_dir']:
        cmd = 'mkdir ' + params['pca_mef_dir']
        printc(cmd, 'info')
        os.system(cmd)

    if 'method' not in params.keys():
        # Binning method (either TIMEBINS or CLUSTER)
        params['method'] = 'CLUSTER'

    if 'pca-sample' not in params.keys():
        # Which stars
        params['pca-sample'] = 'TELLURICS'  # 'ALL' or 'TELLURICS'


    if 'waveref_file_name' not in params.keys():
        params['waveref_file_name'] = '/'.join(params['template_string'].split('/')[:-1])+'/waveref.fits'

    params['whoami'] = getpass.getuser()


    if 'NIRPS_HE' in params['instrument'].upper():
        params['path_to_red'] = '/cosmos99/nirps/apero-data/nirps_he_online/objects/*/*t.fits'
        params['WAVEFILE_KEY'] = 'WAVEFILE'
        params['ext_wavefile'] = 1
        params['fiber'] = 'A'
    elif 'NIRPS_HA' in params['instrument'].upper():
        params['path_to_red'] = '/cosmos99/nirps/apero-data/nirps_ha_online/objects/*/*t.fits'
        params['WAVEFILE_KEY'] = 'WAVEFILE'
        params['ext_wavefile'] = 1
        params['fiber'] = 'A'
    elif 'SPIROU' in params['instrument'].upper():
        params['path_to_red'] = '/cosmos99/spirou/apero-data/spirou_offline/objects/*/*t.fits'
        params['WAVEFILE_KEY'] = 'WAVEFILE'
        params['ext_wavefile'] = 1
        params['fiber'] = 'AB'
    else:
        raise ValueError(f'Instrument {params["instrument"]} not recognized')


    if 'doplot' not in params.keys():
        params['doplot'] = False
    if 'force' not in params.keys():
        params['force'] = False




    if params['instrument'].upper() == 'SPIROU':
        params['INSTRUMENT'] = 'SPIROU'
        params['DATA_DIR'] = '/space/spirou/LBL-PCA/SPIROU'
        inst_folder = 'spirou'
        rsync_machine = 'spirou-client@maestria'
        inst_short = 'spirou'

    if 'NIRPS' in params['instrument'].upper():
        params['INSTRUMENT'] = 'NIRPS_HE'
        params['DATA_DIR'] = '/space/spirou/LBL-PCA/NIRPS_HE'
        inst_folder = 'nirps_he'
        rsync_machine = 'nirps-client@maestria'
        inst_short = 'nirps'


    if params['whoami'] == 'eartigau':
        params['path0'] = '/Users/eartigau/glitch_fp/'+params['instrument']
        if not os.path.exists(params['path0']):
            os.system('mkdir -p '+params['path0'])
        params['output_slinky'] = os.path.join('/Users/eartigau/glitch_fp/', params['instrument'], 'data_'+params['instrument']+'_output/')
        if not os.path.exists(params['output_slinky']):
            os.system('mkdir -p '+params['output_slinky'])

    if params['whoami'] == 'spirou':
        params['path0'] = '/space/spirou/SLINKY/'+params['instrument']
        if not os.path.exists(params['path0']):
            os.system('mkdir -p '+params['path0'])
        params['output_slinky'] = os.path.join('/space/spirou/SLINKY/', params['instrument'], 'data_'+params['instrument']+'_output/')
        if not os.path.exists(params['output_slinky']):
            os.system('mkdir -p '+params['output_slinky'])


    params['search_t_slinky_path'] = os.path.join(params['output_slinky'], '*_slinky/*t_slinky.fits')


    now = datetime.datetime.now()
    now_str = (now.isoformat().replace('-','').replace(':','').replace('T','_')).split('.')[0]

    path_to_summary = '/space/spirou/LBL-PCA/wraps/batch_summary/'+params['INSTRUMENT']
    path_to_summary = path_to_summary+'/'+now_str
    os.system('mkdir '+path_to_summary)

    params['path_to_summary'] = path_to_summary
    params['now_str'] = now_str

    all_files = glob.glob(params['search_t_slinky_path'])
    # loop over all files and check that they are fine
    for f in all_files:
        try:
            os.path.getsize(f)
        except Exception as e:
            print(f'File {f} is corrupted: {e}')
            os.remove(f)

    params['plot_dir'] = os.path.join(params['path_to_summary'], 'plots')
    if not os.path.exists(params['plot_dir']):
        printc(f'Creating directory {params["plot_dir"]}', 'info')
        os.makedirs(params['plot_dir'], exist_ok=True)

    path_key = ['calib_dir','residual_path','patched_wavesol']
    for key in path_key:
        if key in params.keys():
            if not os.path.exists(params[key]):
                printc(f'Creating directory {params[key]}', 'info')
                os.makedirs(params[key], exist_ok=True)

    return params

def dump_error(params,text_comment = ''):
    with open(params['err_file'], 'a') as f:
        # timestamp
        time_iso = datetime.datetime.now().isoformat()
        f.write('\n\t' + time_iso + '\n')
        f.write(text_comment + '\n')
        f.write('\n')

def create_residuals(params):
    # Get all files matching the search path pattern
    all_files = np.array(glob.glob(params['search_t_slinky_path']))

    # Combine objects of interest and hot stars into one array
    good_objs = np.append(params['object_of_interest'], params['hotstars'])

    # Initialize a boolean array to keep track of files to keep
    keep = np.zeros(len(all_files), dtype=bool)
    for ifile, file in enumerate(all_files):
        for obj in good_objs:
            if '/'+obj+'_' in file:
                keep[ifile] = True

    # Filter the files to keep only those of interest
    all_files = all_files[keep]

    # Initialize arrays for MJD, object names, and template existence
    mjds = np.zeros(len(all_files))
    objs = np.zeros(len(all_files), dtype='U99')
    template_exists = np.zeros(len(all_files), dtype=bool)
    tbl_summary = Table()

    # Initialize a list to collect errors
    errors = []
    keep = np.zeros(len(all_files), dtype=bool)

    fiber_setup = params['fibresetup']

    def process_file(i):
        try:
            name_hdr = params['residual_path'] + all_files[i].split('/')[-1].replace('.fits', '_hdr.pick')
            if not os.path.exists(name_hdr):
                data_dict = read_t(all_files[i])
                h = dict(data_dict[f'Flux{fiber_setup}_header'])
                save_pickle(name_hdr, h)
            else:
                h = read_pickle(name_hdr)

            mjds[i] = h[params['key_mjd']]
            objs[i] = h[params['key_obj']]
            template_exists[i] = os.path.exists(params['template_string'].format(objs[i]))

            nn = str(len(all_files))
            ii = str(i).zfill(len(nn))

            printc(f'hdr -> pickle [{ii}/{nn}] {all_files[i]}', 'info')
            if (h[params['key_obj']] in good_objs) and template_exists[i]:
                keep[i] = True
            if h[params['key_obj']] in  params['hotstars']:
                keep[i] = True

        except KeyboardInterrupt:
            printc('We have a keyboard interrupt', 'bad2')
            raise KeyboardInterrupt
        except Exception as e:
            errors.append(f'error reading header for {all_files[i]}\n')
            dump_error(params, f'error reading header for {all_files[i]}\n')

    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = [executor.submit(process_file, i) for i in range(len(all_files))]
        for future in as_completed(futures):
            future.result()


    # Print errors if any
    if len(errors) > 0:
        printc('errors reading headers:', 'bad3')
        for ierr in errors:
            printc(ierr, 'bad3')

    # Filter the files, MJDs, and object names to keep only those with templates
    all_files = all_files[keep]
    mjds = mjds[keep]
    objs = objs[keep]

    # Generate output file names for residuals
    outname = [params['residual_path'] + f.split('/')[-1].replace('.fits', '_residual.fits') for f in all_files]
    outname = [outname[i].replace('.fits', '.pick') for i in range(len(outname))]

    # Populate the summary table with file information
    tbl_summary['file'] = all_files
    tbl_summary['residual_file'] = outname
    tbl_summary[params['key_mjd']] = mjds
    tbl_summary[params['key_obj']] = objs
    tbl_summary['hotstar'] = [obj in params['hotstars'] for obj in objs]

    # Calculate the median day offset
    uu = np.arange(0, 1.1, 0.1)
    u = np.unique(np.round(tbl_summary[params['key_mjd']] % 1, 1))
    midday = np.nanmedian(uu[[i not in u for i in uu]])

    # Add day, SNR, SHAPE_DX, SHAPE_DY, and unique_obj_day columns to the summary table
    tbl_summary['day'] = np.round(tbl_summary[params['key_mjd']] - midday).astype(int)
    tbl_summary['SNR'] = np.zeros(len(tbl_summary), dtype=float)
    tbl_summary['SHAPE_DX'] = np.zeros(len(tbl_summary), dtype=float)
    tbl_summary['SHAPE_DY'] = np.zeros(len(tbl_summary), dtype=float)
    tbl_summary['unique_obj_day'] = np.zeros(len(tbl_summary), dtype='U99')

    # Populate the unique_obj_day column
    for i in range(len(tbl_summary)):
        str1 = tbl_summary[params['key_obj']][i]
        str2 = tbl_summary['day'][i].astype(str)
        uobj = '{}_{}'.format(str1, str2)
        tbl_summary['unique_obj_day'][i] = uobj

    # randomize the order of the files
    tbl_summary = tbl_summary[np.random.permutation(len(tbl_summary))]
    

    # Process each file to compute residuals
    for i in range(len(tbl_summary)):
        file = tbl_summary['residual_file'][i]

        if os.path.exists(tbl_summary['residual_file'][i]):
            # Read the residual file if it already exists
            tmp = read_pickle(tbl_summary['residual_file'][i])
            h = tmp['header']
            tbl_summary['SNR'][i] = h['SNR']
            tbl_summary['SHAPE_DX'][i] = h['SHAPE_DX']
            tbl_summary['SHAPE_DY'][i] = h['SHAPE_DY']

            printc(f'[{i+1}/{len(tbl_summary)}] {file} exists', 'info')

            continue

        # Read the data and header from the file
        data_dict = read_t(tbl_summary['file'][i])
        h = data_dict['FluxA_header']
        tbl_summary['SHAPE_DX'][i] = h['SHAPE_DX']
        tbl_summary['SHAPE_DY'][i] = h['SHAPE_DY']

        data_dict = read_t(tbl_summary['file'][i])
        sp = np.array(data_dict['FluxA'], dtype=float)

        # Get the splined template for non-hot stars
        template,params = get_splined_template(tbl_summary['file'][i],  params = params)

        with warnings.catch_warnings(record=True) as _:
            residual_multiplicative = np.log(sp /template)

        # Apply a low-pass filter to each order
        for iord in snail(range(sp.shape[0]), desc='order', leave=False):
            with warnings.catch_warnings(record=True) as _:
                residual_multiplicative[iord, :] -= lowpassfilter(residual_multiplicative[iord, :], params['hpwidth'])

        # Compute the SNR
        with warnings.catch_warnings(record=True) as _:
            h['SNR'] = 1 / sigma(sp.ravel() - np.roll(sp, 1, axis=1).ravel())
        tbl_summary['SNR'][i] = h['SNR']

        # Save the residuals and header to a pickle file
        tmp = dict()
        tmp['residual_multiplicative'] = residual_multiplicative
        tmp['header'] = dict(h)
        printc(f'[{i+1}/{len(tbl_summary)}] writing {file}', 'info')
        save_pickle(tbl_summary['residual_file'][i], tmp)

    # we loop on objects and make residual plots in the plots directory
    all_objs = np.unique(tbl_summary[params['key_obj']])
    for obj in all_objs:
        printc(f'Plotting residuals for {obj}', 'info')
        tbl_obj = tbl_summary[tbl_summary[params['key_obj']] == obj]
        tbl_obj = tbl_obj[np.argsort(tbl_obj[params['key_mjd']])]
        wave_tmp = fits.getdata(tbl_obj['file'][0],2)[params['sample_order']]
        wplot = 1200
        map_res = np.zeros([len(tbl_obj['residual_file']),wplot], dtype=float)
        for i,f in snail(enumerate(tbl_obj['residual_file']), desc='files', leave=False):
            res = read_pickle(f)['residual_multiplicative'][params['sample_order']]
            map_res[i] = res[len(res)//2-wplot//2:len(res)//2+wplot//2]

        plt.imshow(map_res,vmin = -0.02, vmax = 0.02, aspect='auto', origin='lower',
                   extent=[wave_tmp[len(wave_tmp)//2-wplot//2],wave_tmp[len(wave_tmp)//2+wplot//2],0,len(tbl_obj['residual_file'])], cmap='RdBu_r')
        plt.colorbar()
        plt.title(f'Residuals for {obj}')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Nth File')
        plt.savefig(os.path.join(params['plot_dir'], f'residuals_{obj}.png'))
        plt.close()

    # Generate output file names for PCA-corrected files
    files = tbl_summary['file']
    if params['Npca_adjust']:
        outnames = [f.split('/')[-1].replace('.fits', '_' + params['batchname'] + '.fits') for f in files]
    else:
        outnames = [f.split('/')[-1].replace('.fits', '_' + params['batchname'] + '.fits') for f in files]

    # Create directories for PCA-corrected files if they do not exist
    for i in snail(range(len(files)), desc='files', leave=False):
        path = '/'.join(files[i].split('/')[:-1]) + '_{}/'.format(params['batchname'])
        if not os.path.exists(path):
            printc('creating directory {}'.format(path), 'info')
            os.mkdir(path)
        outnames[i] = path + outnames[i]

    tbl_summary['pca_corrected_file'] = outnames

    # Add a time bin column to the summary table
    tbl_summary['TIMEBIN'] = np.array(
        (tbl_summary[params['key_mjd']] - np.min(tbl_summary[params['key_mjd']])) // params['time_bin'], dtype=int)
    tbl_summary = tbl_summary[np.argsort(tbl_summary[params['key_mjd']])]

    # Fetch the astrometric database and identify telluric stars
    astrom_db = fetch_astrom()
    tellurics = astrom_db['OBJNAME'][astrom_db['TELLURIC'] == True]
    tbl_summary['TELLURIC'] = [obj in tellurics for obj in tbl_summary[params['key_obj']]]

    # Find the optimal number of clusters for the Birch algorithm
    X = np.array([tbl_summary['SHAPE_DX'], tbl_summary['SHAPE_DY']]).T
    for nclusters in range(2, params['MAX_CLUSTERS']):
        model = Birch(threshold=0.01, n_clusters=nclusters)
        model.fit(X)
        yhat = model.predict(X)
        clusters = np.unique(yhat)

        max_span = []
        for cluster in clusters:
            row_ix = np.where(yhat == cluster)
            max_span.append(np.max(np.max(X[row_ix], axis=0) - np.min(X[row_ix], axis=0)))
            if params['doplot']:
                plt.scatter(X[row_ix, 0], X[row_ix, 1])
        if params['doplot']:
            plt.show()

        opt = nclusters, np.max(max_span), params['CUT_MAX_SPAN_CLUSTER']
        printc('nclusters = {}, max_span = {:.2f} pix, max = {:.2f}'.format(*opt), 'info')
        if np.max(max_span) < params['CUT_MAX_SPAN_CLUSTER']:
            break

    tbl_summary['cluster'] = yhat

    # Add a column to indicate objects of interest
    if 'object_of_interest' in params.keys():
        tbl_summary['object_of_interest'] = [obj in params['object_of_interest'] for obj in
                                             tbl_summary[params['key_obj']]]

    # Save the summary table to a CSV file
    tbl_summary.write('{}tbl_summary.csv'.format(params['residual_path']), overwrite=True)

    return tbl_summary


def apply_pca_pix(tbl, params, force=False, doplot=False):
    # to plot or not to plot, that is the question
    params['doplot'] = doplot

    fiber_setup = params['fibresetup']

    keep = tbl['object_of_interest']

    if True not in keep:
        printc('No object of interest in the table', 'bad2')
        return

    tbl0 = Table(tbl[keep])

    # read the first file to get the wavelength solution
    data_dict = read_t(tbl0['file'][0])
    # this is only for the title of the plots and printout for user,
    # it has no impact on the code
    wavemap = data_dict['Wave'+fiber_setup]
    sz = wavemap.shape

    pca_bins = np.unique(tbl0['cluster']).astype(int)

    
    for pca_bin in pca_bins:
        tbl1 = tbl0[tbl0['cluster'] == pca_bin]

        all_pcas = dict()
        for iord in range(sz[0]):
            pca_name = 'IORD{}CLUSTER{}-{}'.format(iord, pca_bin, params['batchname'])
            pca_outname = '{}{}_{}_{}_{}.fits'.format(params['pca_mef_dir'],
                                                      params['instrument'],
                                                      pca_bin, iord, params['batchname'])

            if os.path.exists(pca_outname):
                all_pcas[pca_name] = mef2dict(pca_outname)
            else:
                all_pcas[pca_name] = None

        pca = EMPCA(n_components=params['Npca'])

        for i in range(len(tbl1)):
            file = tbl1['file'][i]
            outname = tbl1['pca_corrected_file'][i]

            if not force and os.path.exists(outname):
                printc('file {} exists, skipping'.format(outname), 'bad1')
                continue

            data = read_pickle(tbl1['residual_file'][i])['residual_multiplicative']
            w = np.ones_like(data)
            bad = ~np.isfinite(data)
            w[bad] = 1e-6*sigma(data[~bad])
            data[bad] = 0

            recon = np.zeros_like(data)

            for iord in snail(range(sz[0]), desc='orders', leave=False):
                pca_name = 'IORD{}CLUSTER{}-{}'.format(iord, pca_bin, params['batchname'])

                if all_pcas[pca_name] is None:
                    # printc('PCA for {} does not exist, skipping'.format(pca_name), 'bad1')
                    continue
                else:
                    dict0 = all_pcas[pca_name]
                    dict1 = dict(dict0)
                    for key in dict0.keys():
                        dict1[key.lower()] = dict0[key]
                    dict1['mean_'] = np.zeros(sz[1])
                    pca.__dict__ = dict1
                    try:
                        recon[iord] = pca.reconstruct(data[iord].reshape(-1, 1).T, w[iord].reshape(-1, 1).T)[0]
                    except:
                        # if a keyword interupt then exit
                        if sys.exc_info()[0] == KeyboardInterrupt:
                            printc('We have a keywoard interrupt', 'bad2')
                            raise KeyboardInterrupt

                        dump_error(params)
                        # printc('Error in PCA reconstruction', 'bad1')
                        continue

            with warnings.catch_warnings(record=True) as _:
                corr = np.exp(recon)
            corr[corr < 0.01] = np.nan

            try:  # TODO change to try
                data_dict = read_t(file)
                data_dict['Flux'+fiber_setup] /= corr

                #write_t(data_dict, outname)

                fibresetup = params['fibresetup']

                # which instrument
                if params['instrument'].upper() == 'SPIROU':
                    wavefile = data_dict[f'Flux{fibresetup}_header']['WAVEFILE']

                if 'NIRPS' in params['instrument'].upper():
                    wavefile = data_dict[f'Flux{fibresetup}_header']['WAVEFILE']

                patched_wavefile = params['patched_wavesol'] + wavefile.split('/')[-1]

                if os.path.exists(patched_wavefile):
                    printc('Reading {}'.format(patched_wavefile), 'info')
                    printc('Updating the wavelength solution in the data_dict', 'info')
                    wave = fits.getdata(patched_wavefile)
                    data_dict[f'Wave{fibresetup}'] = wave

                    #outname2 = outname.replace('TELLU', 'SLINKY')
                    path2 = '/'.join(outname.split('/')[:-1])
                    if not os.path.exists(path2):
                        cmd = 'mkdir {}'.format(path2)
                        printc(cmd, 'info')
                        os.system(cmd)
                    printc('Writing {}'.format(outname), 'info')

                    write_t(data_dict, outname)

                else:
                    dump_error(params, 'problem with the wavefile {}'.format(patched_wavefile))
                    printc('The file {} does not exist'.format(patched_wavefile), 'bad2')

                #printc('Writing {}'.format(outname), 'info')
            except:
                printc('Error writing {}'.format(outname), 'bad2')
                # We print the error message :
                printc(sys.exc_info()[0])
                # if a keyword interupt then exit
                if sys.exc_info()[0] == KeyboardInterrupt:
                    printc('We have a keywoard interrupt', 'bad2')
                    raise KeyboardInterrupt

                continue


def create_symlinks(params):
    instrument = params['instrument']
    printc(instrument)

    if instrument.upper() == 'SPIROU':
        inst = 'spirou'
        inst_short = 'spirou'
        onoff_line = 'offline'
    if 'NIRPS' in instrument.upper():
        inst = 'nirps_he'
        inst_short = 'nirps'
        onoff_line = 'online'

    path_links = f'/space/spirou/LBL-PCA/' + inst.upper() + '/science/{}/'

    obj_folders = np.append(params['object_of_interest'], params['hotstars'])

    files = []
    for obj in obj_folders:
        path = f'/cosmos99/{inst_short}/apero-data/{inst}_{onoff_line}/objects/{obj}/*t.fits'

        files_tmp =  np.array(glob.glob(path))
        files = np.append(files,files_tmp)


    valid_files = []
    for file in files:
        skip = True
        for obj in obj_folders:
            if obj in file:
                skip = False
                folder_obj = obj
                break
        if skip:
            printc('Skip')
            continue

        if not os.path.exists(path_links.format(folder_obj)):
            printc('creating directory {}'.format(path_links.format(folder_obj)), 'info')
            os.mkdir(path_links.format(folder_obj))

        try:
            os.path.getsize(file)

            link_name = path_links.format(folder_obj)+file.split('/')[-1]
            if not os.path.exists(link_name):
                cmd = 'ln -s {} {}'.format(file, link_name)
                printc(cmd, 'info')
                # create the symlink
                os.system(cmd)
            else:
                printc('Link {} already exists'.format(link_name), 'bad1')
            valid_files.append(file)

        except:
            printc(f'File {file} is corrupted', 'bad2')
            os.remove(file)

    templates = glob.glob(f'/cosmos99/{inst_short}/apero-data/{inst}_{onoff_line}/red/other/Template_*_tellu_obj_*.fits')
    for i in range(len(templates)):
        template0 = templates[i]
        # check if the any objs in template0
        if True not in ([obj in template0 for obj in obj_folders]):
            continue

        obj = obj_folders[ np.where([obj in template0 for obj in obj_folders])[0] ][0]

        template_path = params['template_string'].format(obj)

        cmd = f'cp {template0} {template_path}'
        printc(cmd, 'info')
        os.system(cmd)

    waveref_files = glob.glob(f'/cosmos99/{inst_short}/apero-data/{inst}_{onoff_line}/calib/*_pp_e2dsff_A_wavesol_ref_*.fits')
    waveref_file_name = params['waveref_file_name']
    for i in range(len(waveref_files)):
        waveref_file = waveref_files[i]
        cmd = f'cp {waveref_file} {waveref_file_name}'
        printc(cmd, 'info')
        os.system(cmd)



def wrap(params):
    # create the table with the list of files and the corresponding output names
    tbl = create_residuals(params)
    
    # create the PCA models for the residuals
    compute_pca_pix(tbl, params)

    # we add preclean to the batch name to indicate that the files have been cleaned with PCA
    apply_pca_pix(tbl, params)




def plot_wrap_obj():
    params = get_params()
    # create a folder for the PCA files
    if getpass.getuser() == 'eartigau':
        path = '/Users/eartigau/Dropbox/test_pixpca/{}/'.format(params['instrument'])
    if getpass.getuser() == 'spirou':
        path = '/space/spirou/LBL-PCA/wrap_summary/{}/'.format(params['instrument'])

    objs = params['object_of_interest']
    # pick one object and get user input
    printc('Choose an object:', 'number')
    for i, obj in enumerate(objs):
        printc('[{:3.0f}] {}'.format(i, obj), 'number')
    i = int(input('Enter the number of the object: '))
    obj = objs[i]

    tbl_quality = Table()
    tbl_quality['obj'] = objs
    tbl_quality['Median error'] = np.zeros(len(objs))
    tbl_quality['RMS'] = np.zeros(len(objs))

    files = glob.glob(path + 'lbl_*{}*.rdb'.format(obj))

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    for file in files:
        tbl = Table.read(file, format='ascii')
        tbl['vrad'] -= np.nanmedian(tbl['vrad'])
        plt.plot_date(tbl['plot_date'], tbl['vrad'], '.')
        plt.errorbar(tbl['plot_date'], tbl['vrad'], yerr=tbl['svrad'], fmt='.', label=file.split('/')[-1])
    plt.legend()
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.xlabel('Date')
    plt.ylabel('RV [m/s]')
    plt.title(obj)
    plt.tight_layout()
    plt.savefig(path + 'RV_{}.pdf'.format(obj))
    plt.show()


def compute_pca_pix(tbl, params, force=False):
    keep = tbl['TELLURIC']
    tbl0 = Table(tbl[keep])

    if len(tbl0) <=2*params['Npca']:
        printc('Not enough files for PCA', 'bad2')
        return

    fibre_setup = params['fibresetup']

    print(len(tbl0))
    print(tbl0)

    print(tbl0['file'][0])

    # read the first file to get the wavelength solution
    wavemap = read_t(tbl0['file'][0])['Wave'+fibre_setup]

    params['sample_blaze'] = read_t(tbl0['file'][0])['Blaze'+fibre_setup]
    for iord in range(wavemap.shape[0]):
        params['sample_blaze'][iord] /= np.nanmedian(params['sample_blaze'][iord])

    # this is only for the title of the plots and printout for user,
    # it has no impact on the code
    wave_middle = np.nanmedian(wavemap, axis=1)
    # get the size of the wavelength map and overall size
    sz = wavemap.shape

    # check that all files have the corresponding residual file
    keep = np.zeros_like(tbl0, dtype=bool)
    for i in range(len(tbl0)):
        keep[i] = os.path.exists(tbl0['residual_file'][i])

    # print the number of files that will be corrected and the number of objects
    # that were passed initially. Normally these two are the same.
    printc('Keeping {} files out of {}'.format(np.sum(keep), len(tbl0)), 'info')
    tbl0 = tbl0[keep]

    # loop on the time bins and create the list of files to be included in
    # the PCA. The PCA is computed on files +-1 time bin from the time bin
    # considered at each moment

    if params['method'] == 'TIMEBINS':
        tbl0['PCABIN'] = tbl0['TIMEBIN']
        tbl['PCABIN'] = tbl['TIMEBIN']
    elif params['method'] == 'CLUSTER':
        tbl0['PCABIN'] = tbl0['cluster']
        tbl['PCABIN'] = tbl['cluster']

    pca_bins = np.unique(tbl0['PCABIN'].data)

    for iord in range(sz[0]):
        params['sample_blaze'][iord] /= np.nanmedian(params['sample_blaze'][iord])

    outdir =  params['pca_mef_dir']
    instrument = params['instrument']

    for ipca_bin in range(len(pca_bins)):
        # We leave open the possibility that bins are either in 'clusters'
        # or in 'timebins'.
        printc('~')

        pca_bin = pca_bins[ipca_bin]

        if params['method'] == 'TIMEBINS':
            # loop on the time bins and create the list of files to be included in
            # the PCA
            sub_table = tbl0[np.abs(tbl0['PCABIN'] - pca_bin) <= 1]
        elif params['method'] == 'CLUSTER':
            sub_table = tbl0[tbl0['PCABIN'] == pca_bin]

        if params['Npca_adjust']:
            params['Npca'] = len(np.unique(sub_table['day'])) // 2
            if params['Npca'] < 2:
                params['Npca'] = 2

            if params['Npca'] > params['npca-max']:
                params['Npca'] = params['npca-max']

        printc('Npca = {}'.format(params['Npca']), 'info')

        # files are both tellurics and in pca bin
        pca_files = (tbl0['PCABIN'] == pca_bin) & tbl0['TELLURIC']
        pca_table = tbl0[pca_files]

        # faint-target files
        faint_files = (tbl['PCABIN'] == pca_bin) & ~tbl['TELLURIC']
        faint_table = tbl[faint_files]



        # Number of unique combination of objects * dates
        unique_obj_day_pca = np.unique(pca_table['unique_obj_day'])
        Nobjs_pca = len(unique_obj_day_pca)

        unique_obj_day_faint = np.unique(faint_table['unique_obj_day'])
        Nobjs_faint = len(unique_obj_day_faint)

        # list the number of files that will be corrected
        print_outputs = (ipca_bin, len(pca_bins), len(pca_table), Nobjs_pca, Nobjs_faint)

        printc('UPcaBins {}/{}, {} files in PCA, {} Unique objs+nights, {} unique faint'.format(
            *print_outputs), 'number')

        if Nobjs_pca < 2 * params['Npca']:
            printc('Not enough unique objects+nights for PCA', 'bad2')
            continue

        if not os.path.isdir(params['pca_mef_dir']):
            # create the directory if it does not exist
            cmd = 'mkdir ' + params['pca_mef_dir']
            printc(cmd, 'info')
            os.system(cmd)

        instrument = params['instrument']
        for iord in range(sz[0]):
            pca_outname = '{}{}_{}_{}_{}.fits'.format(params['pca_mef_dir'],
                                                      instrument,
                                                      pca_bin, iord, params['batchname'])

            if os.path.exists(pca_outname) and not force:
                printc('file {} exists, skipping'.format(pca_outname), 'bad1')
                continue

            # create the map of the residuals. This has the size of the (pca_table * 4088)
            # in the time bin. Sample2 is a ravel() version of the spectrum.
            sample2 = np.zeros([Nobjs_pca, sz[1]])
            # corresponding weights
            sigmap = np.zeros([Nobjs_pca, sz[1]])

            n_per_day = np.zeros(Nobjs_pca)
            # loop on the files and pad the residuals in a ravel() format
            for i in snail(range(len(pca_table)), desc='files, iord {}/{}'.format(iord, sz[0]), leave=False):
                tmp = read_pickle(pca_table['residual_file'][i])['residual_multiplicative'][iord]
                im = tmp * params['sample_blaze'][iord]

                ipca = np.where(unique_obj_day_pca == pca_table['unique_obj_day'][i])[0][0]
                # we add the residual on the corresponding line for the unique object+night
                sample2[ipca] += im
                n_per_day[ipca] += 1

            for i in range(Nobjs_pca):
                sample2[i] /= n_per_day[i]

            for i in range(Nobjs_pca):
                sigmap[i] = sigma(sample2[i] - np.roll(sample2[i], 1))

            # weights inverse to the sigma, *not* the sigma**2
            with warnings.catch_warnings(record=True) as _:
                weights2 = 1 / sigmap

            # remove the bad pixels from sample2 and weights2
            sample3 = np.array(sample2)

            with warnings.catch_warnings(record=True) as _:
                # any point with a residual value of >0.5 (remembre, this is in log space)
                # is really bad. This corresponds to an error of ~30%
                cond1 = np.abs(sample3) < 0.5
                cond2 = np.isfinite(weights2)

                # is not meeting one or the other
                bad = ~cond1 | ~cond2

            weights2[bad] = np.nanstd(weights2[~bad])*1e-6
            sample2[bad] = 0

            try:
                pca = EMPCA(n_components=params['Npca'])

                with warnings.catch_warnings(record=True) as _:
                    fit_pca = pca.fit(sample2, weights=weights2)
                    dict2mef(pca_outname, fit_pca.__dict__)
                printc('order {:3.0f}/{:3.0f} - npca = {}'.format(iord, sz[0], params['Npca']), 'number')
                printc('We write {}'.format(pca_outname), 'info')

            except:
                # exit code if this is a keyboard interrupt
                if sys.exc_info()[0] == KeyboardInterrupt:
                    printc('We have a keywoard interrupt', 'bad2')
                    return
                printc('PCA failed for order {}'.format(iord), 'bad2')


def get_residual_map(files, params, residual_type = 'multiplicative'):

    if residual_type not in ['multiplicative', 'additive']:
        printc('Invalid residual type', 'bad2')
        raise ValueError

    iord = params['iord']
    hpwidth = params['hpwidth']

    fiber_setup = params['fibresetup']

    residuals = np.zeros([len(files), 4088])
    for i in snail(range(len(files)), desc='Computing residuals', leave=False):
        data_dict = read_t(files[i])

        sp = data_dict['Flux'+fiber_setup]
        model,params = get_splined_template(files[i], params=params)

        sp = sp[iord]
        model = model[iord]

        if type == 'additive':
            ratio = lowpassfilter(sp / model, hpwidth)
            model*=ratio
            sp -= model
            residuals[i] = sp
        else:
            ratio = np.log(sp / model)
            ratio -= lowpassfilter(ratio, hpwidth)
            residuals[i] = np.log(ratio)

    return residuals


def plot_residuals():
    params = get_params('/Users/eartigau/pycodes/pixeldecorr/yamls/params_tellu15B.yaml')
    params['iord'] = 50

    fiber_setup = params['fibresetup']

    obj = 'GJ707'

    path1 = '/Volumes/courlan/lbl_NIRPS_HE/science/{}_PCAx/*{}*.fits'.format(obj, params['batchname'])
    files1 = np.sort(glob.glob(path1))
    # files1 = files1[:100]

    path2 = '/Volumes/courlan/lbl_NIRPS_HE/science/{}/*.fits'.format(obj)
    files2 = np.sort(glob.glob(path2))

    w = read_t(files1[0])['Wave'+fiber_setup]

    n = np.min([len(files1), len(files2)])
    files1 = files1[:n]
    files2 = files2[:n]
    date = [read_t(f)[f'Flux{fiber_setup}_header']['DATE'].split('T')[0] for f in files1]
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

    # printc(np.nanmedian(w[params['iord']]),'info')
    fig, ax = plt.subplots(2, 1, figsize=(10, 5), sharex='all', sharey='all')
    ax[1].imshow(residuals1, aspect='auto', vmin=-0.05, vmax=0.05)
    ax[1].set(title='PCA')
    ax[0].imshow(residuals2, aspect='auto', vmin=-0.05, vmax=0.05)
    ax[0].set(title='No PCA')

    for i in range(len(date) - 1):
        if date[i] == date[i + 1]:
            continue

        ax[0].hlines(i + .5, 0, 4088, color='w', alpha=0.5)
        ax[1].hlines(i + .5, 0, 4088, color='w', alpha=0.5)

    plt.show()

# pca_pix('NIRPS', doplot=False)
