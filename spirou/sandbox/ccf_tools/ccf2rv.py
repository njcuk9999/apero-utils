import numpy as np
import os
import glob
from astropy.io import fits
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.table import Table
from bisector import *
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from fits2wave import *
import warnings
import sys
#
# we get the user name. as we are just a few people in the team, we could all
# have our paths hard coded here.
#
name = os.popen('whoami').read().split('\n')[0]

if name == 'eartigau':
    PATH = '/Users/eartigau/wrap_drs_rv'

elif name == 'clairemoutou':
    PATH = '/Users/clairemoutou/apero/data'

elif name=='andres':
    PATH = '/some/folder/somewhere'

elif name == 'pascal':
    PATH = '/some/folder/somewhere'

else:
    default: PATH = ''

def gauss(v,v0,ew,zp,amp):
    # gaussian with a constant offset. As we know that the ccfs are negative structures, amp will be negative
    return zp+amp*np.exp( -0.5*(v-v0)**2/ew**2)

def dispatch_object(object, all_ccf_dir = 'all_ccfs'):
    # find all ccf files and copy those that match your object into a folder with the object's name
    #
    # set object to "all" to dispatch all files in folders
    #
    all_ccf_files = glob.glob(PATH+'/'+all_ccf_dir+'/*ccf*AB.fits')

    ngood_AB = 0
    ngood_C = 0

    print('we are dispatching CCF files into the proper directory for {0}'.format(object))
    for file in tqdm(all_ccf_files):
        # read ccf file
        im, h = fits.getdata(file, header=True)
        outname = PATH+'/'+h['DRSOBJN']+'/'+file.split('/')[-1]


        if (h['DRSOBJN'] == object) or (object == 'all'):
            # if the directory with the object name does not exist, we create it
            if os.path.isdir(PATH+'/'+h['DRSOBJN']) == False:
                os.system('mkdir {0}/{1}'.format(PATH,h['DRSOBJN']))
            if os.path.isfile(outname) == False:
                os.system('cp {0} {1}'.format(file,outname))

            if '_AB.fits' in file:
                ngood_AB+=1
            if '_C.fits' in file:
                ngood_C+=1

    print('We found {0} files for that object. It includes {1} AB files and {2} C drift files'.format(ngood_AB+ngood_C,ngood_AB, ngood_C))

def get_object_rv(object,
                  mask = 'sept18_andres_trans50',
                  method = 'template',
                  exclude_orders = [-1],
                  weight_table = '',
                  force = True,
                  snr_min = 0.0,
                  weight_type = '',
                  sanitize = False,
                  bandpass = 'YJHK',
                  velocity_window = 10,
                  dvmax_per_order = 1.0,
                  doplot = True,
                  do_blacklist = False):
    #
    # parameters :
    #
    # object -> name of the object to be analyzed, linked to the folder where the data should be. You need
    #           to use dispatch_object first to put the data where it belongs
    #
    #
    # mask = X -> name of the mask used for the CCF. You can have a number of different masks for the same object
    #           and these should be processed individually
    #
    # if you want to systematically exclude orders, then set exclude_orders = [a,b,c,d]. Leave to -1 if you want to
    #          keep all valid orders
    #
    # method = X -> Method to measured velocity. Add you own! For now we have---
    #
    #          bisector_N_M returns the velocity between N and Mth percentile of line depth.
    #
    #         gaussfit --> fits a gaussian to the mean CCF
    #
    # if you want to force a weight onto orders rather than let the program decide automatically, then provide a
    # value for weight_table (e.g., weigth_table = "my_weights.csv". It must have 49 rows, one column name
    # "WEIGHT" (all capitals) and be read as a proper numpy.table file
    #

    if sanitize:
        sp_type = 'sani'
    else:
        sp_type = 'tcorr'


    ccf_files = np.array(glob.glob('{0}/{1}/*{3}*{2}*.fits'.format(PATH,object,mask,sp_type)))

    if do_blacklist:
        ccf_files = check_blacklist(ccf_files)

    # definition of photometric bandpasses
    #http: // svo2.cab.inta - csic.es / svo / theory / fps3 / index.php?mode = browse & gname = CFHT & gname2 = Wircam & asttype =
    # see values here
    # Y = 938.600	1113.400
    # J = 1153.586	1354.422
    # H = 1462.897	1808.544
    # k = 1957.792	2400 # modified to get to the edge of the SPIRou domain

    # get typical wavelength solution from first file and central wavelength grid per order
    wave_middle = np.nanmean(fits2wave(fits.getheader(ccf_files[0], ext=1)), axis=1)

    keep_orders = np.zeros(49)
    if 'Y' in bandpass:
        keep_orders[(wave_middle>938)*(wave_middle<1113)] = True
    if 'J' in bandpass:
        keep_orders[(wave_middle>1153)*(wave_middle<1354)] = True
    if 'H' in bandpass:
        keep_orders[(wave_middle>1462)*(wave_middle<1808)] = True
    if 'K' in bandpass:
        keep_orders[(wave_middle>1957)*(wave_middle<2400)] = True

    for i in range(49):
        if i in exclude_orders:
            keep_orders[i] = False

    # update orders rejected because of domain as well as input requirements.
    exclude_orders = np.where(keep_orders == False)[0]

    for i in range(49):
        if i in exclude_orders:
            keep_orders[i] = False

    # form a unique batch name with mask, object and method
    batch_name = '{0}/{1}_mask_{2}_{3}'.format(PATH,object,mask,method)

    if force == False:
        if os.path.isfile('{0}.csv'.format(batch_name)):
            return Table.read('{0}.csv'.format(batch_name))

    tbl = Table()# output table to be saved as CSV file with RV measurements
    tbl['FILES'] = ccf_files
    tbl['RV'] = np.zeros_like(ccf_files,dtype = float) # measured RV
    tbl['ERROR_RV'] = np.zeros_like(ccf_files,dtype = float) # measured RV error
    tbl['CCF_RESIDUAL_RMS']= np.zeros_like(ccf_files,dtype = float) # RMS of CCF - median(CCF)

    # keywords from file headers to be added to the CSV table.
    keywords = ['MJDATE','BERV', 'RV_DRIFT','EXTSN035','AIRMASS','TLPEH2O','TLPEOTR','RV_WAVFP','RV_SIMFP','DATE',
                'MJDMID','DATE-OBS','EXPTIME']


    # add method-specific keywords
    if 'bisector' in method:
        tbl['BIS_SLOPE'] = np.zeros_like(ccf_files,dtype = float)  # bisector slope
        tbl['BIS_WIDTH'] = np.zeros_like(ccf_files,dtype = float)  # bisector slope

    # add method-specific keywords
    if 'gaussian' in method:
        tbl['GAUSS_WIDTH'] = np.zeros_like(ccf_files,dtype = float)  # bisector slope
        tbl['GAUSS_AMP'] = np.zeros_like(ccf_files,dtype = float)  # bisector slope
        tbl['GAUSS_ZP'] = np.zeros_like(ccf_files,dtype = float)  # bisector slope

    #... feel free to add columns here


    npy_file = '{0}/{0}_{1}_{2}_ccf_cube.npy'.format(object,sp_type,mask)

    if not os.path.isfile(npy_file):
        print('we load all CCFs into one big cube')
        for i in (range(len(ccf_files))):
            # we loop through all files
            ccf_tbl = fits.getdata(ccf_files[i])

            ccf_RV = ccf_tbl['RV'] # ccf velocity offset, not to be confused with measured RV

            # We must absolutely have always the same RV grid for the CCF. We must check that consecutive ccf_RV are i
            # identical
            if i != 0:
                if np.sum(ccf_RV != ccf_RV_previous):
                    print('We have a big problem! The RV vector of CCF files are not all the same')
                    print('Files {0} and {1} have different RV vectors.'.format(ccf_files[i-1],ccf_files[i]))
                    sys.exit()
            ccf_RV_previous = np.array(ccf_RV)

            print('V[min/max] {0:.1f} / {1:.1f} km/s, file {2}'.format(np.min(ccf_RV),np.max(ccf_RV),ccf_files[i]))
            # if this is the first file, we create a cube that contains all CCFs for all orders for all files
            if i ==0:
                ccf_cube = np.zeros([49,len(ccf_tbl),len(ccf_files)])+np.nan

            # we input the CCFs in the CCF cube
            for j in range(49):
                # if we need to exlude orders, we do it here.
                if j in exclude_orders:
                    continue

                tmp =  ccf_tbl['ORDER'+str(j).zfill(2)]

                if False not in np.isfinite(tmp):
                    # we normalize to a continuum of 1
                    tmp /= np.polyval(np.polyfit(ccf_RV, tmp, 1), ccf_RV)
                    ccf_cube[j,:,i] = tmp

        print('We save {0}, this will speed things up next time you run this code'.format(npy_file))
        np.save(npy_file,ccf_cube)

    else:
        print('We load {0}, this is speedier'.format(npy_file))
        ccf_cube = np.load(npy_file)

        # we need to load the first file just to get the velocity grid
        ccf_tbl = fits.getdata(ccf_files[0])
        ccf_RV = ccf_tbl['RV']

    print('We load values from headers, slower on first run, faster later')
    for i in (range(len(ccf_files))):
        hdr = fits.getheader(ccf_files[i],ext = 1)
        if i ==0:
            # now that we have a first header, we add the relevant columns to the CSV table
            for key in keywords:
                if key in hdr:

                    key_type = type(hdr[key])
                    # if we have a string, we set the table to accept long values (up to 99 characters)
                    if key_type == str:
                        key_type = '<U99'
                else:
                    # keyword not in header, we need to assume something. The safest is string
                    key_type = str

                # add the column to the CSV file
                tbl[key] = np.zeros_like(ccf_files,dtype = key_type)

        for key in keywords:
            if key in hdr:
                tbl[key][i] = hdr[key]


    # we apply the SNR threshold
    keep = tbl['EXTSN035']>snr_min
    tbl = tbl[keep]
    ccf_cube = ccf_cube[:,:,keep]
    ccf_files = ccf_files[keep]

    with warnings.catch_warnings(record=True) as _:
        # some slices in the sum are NaNs, that's OK
        # this is the median CCF  all epochs for the 49 orders
        med_ccf =  np.nanmedian(ccf_cube,axis = 2)

    for iord in range(49):
        if iord not in exclude_orders:
            if not np.isfinite(np.mean(med_ccf[iord,:])):
                print('Order {0} has a CCF full of NaNs. Added to the rejected orders '.format(iord))
                exclude_orders = np.append(exclude_orders,iord)


    # find minimum for CCF. This is used to fit a gaussian to each order and force
    # velocity to zero
    id_min = np.nanargmin( np.nanmedian(med_ccf,axis=0))

    # find if a given CCF is off by more than the pre-defined threshold
    dv_CCF_min = (ccf_RV[np.argmin(med_ccf,axis = 1)] - ccf_RV[id_min])
    bad_orders = dvmax_per_order<np.abs(dv_CCF_min )

    for iord in range(49):
        if iord not in exclude_orders:
            if bad_orders[iord]:
                print('The CCF of order {0} has its minima {1:.2f} km/s from median CCF, above threshold of +-{2:.2f} km/s'.
                      format(iord,dv_CCF_min[iord],dvmax_per_order))
                exclude_orders = np.append(exclude_orders,iord)

    ccf_cube[exclude_orders, :, :] = np.nan


    # find valid pixels to measure CCF properties
    g = np.abs(ccf_RV - ccf_RV[id_min]) < velocity_window

    with warnings.catch_warnings(record=True) as _:
        # some slices in the sum are NaNs, that's OK
        ccf_Q = np.nansum(np.gradient(med_ccf[:,g],axis=1)**2, axis = 1)

    ccf_Q[ccf_Q == 0] = np.nan
    ccf_depth = 1-med_ccf[:,id_min]



    if weight_type == 'DVRMS_CC':
        weights = 1/np.nanmedian(DVRMS_CC,axis=1)**2
        weights[np.isfinite(weights) == False] = 0

        for iord in range(49):
            if iord in exclude_orders:
                weights[iord] = 0

        weights = weights/np.sum(weights)
    else:
        if weight_table == "":
            # now we find the RMS of the Nth spectrum relative to the median
            rms = np.zeros([len(ccf_files),49])
            for i in range(len(ccf_files)):
                with warnings.catch_warnings(record=True) as _:
                # some slices in the median are NaNs, that's OK
                    rms[i,:] = np.nanmedian(np.abs(ccf_cube[:,:,i]-med_ccf),axis=1)
                rms[i, :] /= np.nanmedian(rms[i, :])

            rms[:,exclude_orders] = np.nan

            if doplot:
                plt.imshow(rms)
                plt.xlabel('Nth order')
                plt.ylabel('Nth frame')
                plt.title('RMS of CCF relative to median')
                plt.show()

            with warnings.catch_warnings(record=True) as _:
                # some slices in the sum are NaNs, that's OK
                # this is the typical noise from the ccf dispersion
                ccf_rms = np.nanmedian(rms,axis=0)

            # set to NaN values that are invalid
            ccf_rms[ccf_rms == 0] = np.nan


            # assuming that the CCF has the same depth everywhere, this is the correct weighting of orders
            weights = ccf_Q/ccf_rms**2
            weights[weights == 0] = np.nan
            weights[exclude_orders] = np.nan
            # we normalize the sum of the weights to one
            weights /= np.nansum(weights)

            if doplot:
                fig,ax = plt.subplots(nrows = 3, ncols=1,sharex = True)
                ax[0].plot(weights,'go')
                ax[0].set(title = '{0}, mask {1}'.format(object,mask),xlabel = 'Nth order', ylabel = 'Relative order weight')

                ax[1].plot(ccf_Q,'go')
                ax[1].set(xlabel = 'Nth order', ylabel = 'ccf Q')

                ax[2].plot(1/ccf_rms**2,'go')
                ax[2].set(xlabel = 'Nth order', ylabel = '1/$\sigma_{CCF}^2$')
                plt.tight_layout()
                plt.savefig('{0}_weights.pdf'.format(batch_name))
                plt.show()

            tbl_weights = Table()
            tbl_weights['order'] = np.arange(49)
            tbl_weights['weights'] = weights
            tbl_weights['ccf_depth'] = ccf_depth
            tbl_weights['ccf_Q'] = ccf_Q
            tbl_weights.write('{0}_weights.csv'.format(batch_name),overwrite = True)

        else:
            print('You provided a weight file, we load it and apply weights accordingly')
            tbl_weights = Table.read(weight_table)
            weights = np.array(tbl_weights['weights'],dtype = float)
            weights /= np.nansum(weights)

    if doplot:
        plt.imshow(med_ccf,aspect = 'auto',vmin = 0.8,vmax= 1.05,extent = [np.min(ccf_RV),np.max(ccf_RV),49,0])
        plt.xlabel('Velocity bin [km/s] ')
        plt.ylabel('Nth order')
        plt.title('Median CCF')
        plt.savefig('{0}_medianccf.pdf'.format(batch_name))
        plt.show()

        plt.imshow(ccf_cube[:,:,0]-med_ccf,aspect = 'auto',vmin = -0.1,vmax= 0.1,extent = [np.min(ccf_RV),np.max(ccf_RV),49,0])
        plt.xlabel('Velocity bin [km/s]')
        plt.ylabel('Nth order')
        plt.title('Sample residual CCF map')
        plt.savefig('{0}_residualccf.pdf'.format(batch_name))
        plt.show()

        plt.plot(tbl['MJDATE'], tbl['EXTSN035'], 'g.')
        plt.xlabel('MJDATE')
        plt.ylabel('SNR for order 35\n(around 1.6 $\mu$m)')
        plt.title('Signal-to-noise ratio')
        plt.savefig('{0}_snr.pdf'.format(batch_name))
        plt.show()

    ccf_cube_norm = np.zeros_like(ccf_cube)
    for i in range(49):
        if np.isfinite(weights[i]):
            ccf_cube_norm[i, :, :] = (ccf_cube[i,:,:] * weights[i])


    # get a per-file weighted mean CCF
    mean_ccf = np.nansum(ccf_cube_norm,axis=0)

    if doplot:
        fig,ax = plt.subplots(nrows = 2, ncols = 1)
        for i in range(len(ccf_files)):
            color = [i/len(ccf_files),1-i/len(ccf_files),1-i/len(ccf_files)]
            ax[0].plot(ccf_RV,mean_ccf[:,i],color = color,alpha = 0.2)
            ax[1].plot(ccf_RV, mean_ccf[:, i], color=color, alpha=0.2)

        ax[0].set(xlabel = 'Velocity [km/s]',ylabel = 'CCF depth', title = 'Mean CCFs')
        ax[1].set(xlabel = 'Velocity [km/s]',ylabel = 'CCF depth', title = 'Mean CCFs',xlim = [ccf_RV[id_min]-10,
                                                                                               ccf_RV[id_min]+10])
        plt.tight_layout()
        plt.savefig('{0}_CCFs.pdf'.format(batch_name))
        plt.show()

    #if 'template' in method:
    #print('')

    g = np.abs(ccf_RV - ccf_RV[id_min]) < velocity_window

    rv_prev = np.array(tbl['RV'])
    corr_ccf = np.array(mean_ccf)

    nite_max = 20
    ite = 0
    rms_rv_ite = np.inf
    # we iterate until we have an rms from iteration to iteration of <10 cm/s or we reached a max of 20 iterations
    print('\n')
    while (rms_rv_ite>1e-4) and (ite<nite_max):
        med_corr_ccf = np.nanmedian(corr_ccf,axis=1)

        # normalize continuum to 1
        continuum = np.abs(ccf_RV-ccf_RV[id_min])>velocity_window
        med_corr_ccf/=np.nanmedian(med_corr_ccf[continuum])

        fit = np.polyfit(ccf_RV[continuum], med_corr_ccf[continuum], 2)
        corr = np.polyval(fit, ccf_RV)
        corr -= np.mean(corr)
        med_corr_ccf -= corr


        for i in range(len(ccf_files)):
            spline = ius(ccf_RV,mean_ccf[:,i],ext=3,k=5)
            corr_ccf[:,i] = spline(ccf_RV+tbl['RV'][i])

            # correcting median of CCF
            med =  np.nanmedian(corr_ccf[:,i] - med_corr_ccf)
            mean_ccf[:, i] -= med

            # correcting depth of CCF
            amp = np.nansum( (corr_ccf[:,i] - np.mean(corr_ccf[:,i]))*(med_corr_ccf - np.mean(med_corr_ccf)) )/np.nansum((med_corr_ccf - np.mean(med_corr_ccf))**2)
            mean_ccf[:, i] = (mean_ccf[:,i] - np.mean(mean_ccf[:,i]))/np.sqrt(amp)+np.mean(mean_ccf[:,i])

            # correcting 2rd order polynomial structures in continuum
            fit = np.polyfit(ccf_RV,med_corr_ccf-corr_ccf[:,i],2)
            corr = np.polyval(fit, ccf_RV)
            mean_ccf[:, i] += corr/2

        deriv = np.gradient(med_corr_ccf) / np.gradient(ccf_RV)
        deriv = deriv[g]
        deriv = deriv / np.nansum(deriv ** 2)

        for i in range(len(ccf_files)):
            tbl['RV'][i]-=np.nansum( (corr_ccf[:,i] - med_corr_ccf)[g]*deriv)

        tbl['RV'] -= np.nanmean(tbl['RV'])
        #plt.plot( tbl['RV'],'.')
        rms_rv_ite = np.nanstd(rv_prev - tbl['RV'])
        print('Template CCF iteration number {0:3}, rms RV change {1:3.4f} km/s for this step'.format(ite+1,rms_rv_ite))
        rv_prev = np.array(tbl['RV'])
        ite+=1

    tbl['RV']+=ccf_RV[id_min]

    if doplot:
        fig,ax = plt.subplots(nrows = 2, ncols = 1)
        for i in range(len(ccf_files)):
            color = [i/len(ccf_files),1-i/len(ccf_files),1-i/len(ccf_files)]
            ax[0].plot(ccf_RV,corr_ccf[:,i],color = color,alpha = 0.2)
            ax[1].plot(ccf_RV, corr_ccf[:, i], color=color, alpha=0.2)

        ax[0].set(xlabel = 'Velocity [km/s]',ylabel = 'CCF depth', title = 'Mean CCFs')
        ax[1].set(xlabel = 'Velocity [km/s]',ylabel = 'CCF depth', title = 'Mean CCFs',xlim = [ccf_RV[id_min]-10,
                                                                                               ccf_RV[id_min]+10])
        plt.tight_layout()
        plt.savefig('{0}_template.pdf'.format(batch_name))
        plt.show()

    # we use the bisector method, the name should be something like
    # method = 'bisector_30_70' to get the mean bisector between the 30th and 70th percentile
    if 'bisector' in method:
        # we find the min of CCF and will only compute bisector of +-50 km/s to avoid problems at ccf edges
        imin = np.argmin(np.nanmedian(mean_ccf, axis=1))
        good_RV = np.abs(ccf_RV - ccf_RV[imin]) < 50

        bis_min, bis_max = np.array(method.split('_')[1:3], dtype=float) / 100.

        for i in range(len(ccf_files)):


            try:
                depth, bis, width = bisector(ccf_RV[good_RV], mean_ccf[good_RV,i], low_high_cut=0.2)
                fit = np.polyfit(depth[(depth > bis_min) & (depth < bis_max)] - (bis_min + bis_max) / 2,
                           bis[(depth > bis_min) & (depth < bis_max)], 1)

                tbl['RV'][i] =  fit[1]
                tbl['BIS_SLOPE'][i] =  fit[0]
                tbl['BIS_WIDTH'][i] = np.mean(width[(depth > bis_min) & (depth < bis_max)])
            except:
                print('We had an error with file {0} computing the bisector'.format(ccf_files[i]))
                print('Values will be reported as NaN')
                tbl['RV'][i] =  np.nan
                tbl['BIS_SLOPE'][i] =  np.nan
                tbl['BIS_WIDTH'][i] = np.nan

        fig,ax = plt.subplots(nrows = 2, ncols=1,sharex = True)
        ax[0].plot(tbl['MJDATE'], tbl['RV'], 'g.')
        ax[0].set(title='Velocity',xlabel = 'MJDATE',ylabel = 'RV [km/s]')
        ax[1].plot(tbl['MJDATE'], tbl['BIS_SLOPE'], 'g.')
        ax[1].set(title='Bisector slope',xlabel = 'MJDATE',ylabel = 'slope [km/s/fract. depth]')
        plt.tight_layout()
        plt.savefig('{0}_RV.pdf'.format(batch_name))
        plt.show()

    if 'gaussian' in method:
        imin = np.argmin(np.nanmedian(mean_ccf, axis=1))

        for i in range(len(ccf_files)):
            p0 = [ccf_RV[imin],1,1,-0.1]
            fit, pcov = curve_fit(gauss,ccf_RV,mean_ccf[:,i],p0 = p0)

            tbl['RV'][i] = fit[0]
            tbl['GAUSS_WIDTH'][i] = fit[1]
            tbl['GAUSS_AMP'][i] = fit[3]
            tbl['GAUSS_ZP'][i] = fit[2]

        fig,ax = plt.subplots(nrows = 2, ncols=1,sharex = True)
        ax[0].plot(tbl['MJDATE'], tbl['RV'], 'g.')
        ax[0].set(title='Velocity',xlabel = 'MJDATE',ylabel = 'RV [km/s]')
        ax[1].plot(tbl['MJDATE'], tbl['GAUSS_WIDTH']*2.354, 'g.')
        ax[1].set(title='Gaussian width',xlabel = 'MJDATE',ylabel = 'Gaussian FWHM [km/s]')
        plt.tight_layout()
        plt.savefig('{0}_RV.pdf'.format(batch_name))
        plt.show()



    # we add a measurement of the STDDEV of each mean CCF relative to the median CCF after correcting for the measured
    # velocity. If you are going to add 'methods', add them before this line
    corr_ccf = np.array(mean_ccf)
    for i in range(len(ccf_files)):
        spline = ius(ccf_RV,mean_ccf[:,i],ext=3)
        corr_ccf[:,i] = spline(ccf_RV+tbl['RV'][i]-np.mean(tbl['RV']))

    med_corr_ccf = np.nanmedian(corr_ccf, axis=1)
    g = np.abs(ccf_RV - ccf_RV[id_min]) < 10

    if doplot:
        plt.plot(ccf_RV, med_corr_ccf, color='black', alpha=0.4,label = 'median CCF', linewidth=2)


    # pix scale expressed in CCF pixels
    # SPIRou pixels are about 2.3 km/s
    pix_scale = 2.3/np.nanmedian( np.gradient(ccf_RV))
    for i in range(len(ccf_files)):
        residual = corr_ccf[:,i] - med_corr_ccf
        tbl['CCF_RESIDUAL_RMS'][i] = np.std(residual[g])

        # 1/dvrms -avoids division by zero
        inv_dvrms = (np.gradient(med_corr_ccf) / np.gradient(ccf_RV))/((np.nanstd(residual) * np.sqrt(pix_scale)) )
        tbl['ERROR_RV'][i] = 1 / np.sqrt(np.nansum(inv_dvrms** 2))

        color = [i/len(ccf_files),1-i/len(ccf_files),1-i/len(ccf_files)]
        if doplot:
            plt.plot(ccf_RV,residual+1,color = color,alpha = 0.2)

    if doplot:
        plt.title('Residual CCFs')
        plt.xlabel('velocity [km/s]')
        plt.ylabel('CCF depth')
        plt.legend()
        plt.savefig('{0}_residual_CCF.pdf'.format(batch_name))

        plt.show()

    # output to csv file
    tbl.write('{0}.csv'.format(batch_name),overwrite = True)

    return tbl
