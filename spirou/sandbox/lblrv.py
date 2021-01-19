import numpy as np
import glob
from astropy.io import fits
import os
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from astropy.table import Table
from tqdm import tqdm
from scipy import constants
from datetime import datetime
import etienne_tools as et

def lblrv(obj_sci,obj_template,doplot_ccf = False,doplot_debug = False, force = False,
          lblrv_path = 'lblrv/',mask_path = 'masks/',template_path = 'templates/',
          science_path = 'tellurics/',ref_blaze_file = '2498F798T802f_pp_blaze_AB.fits', check_fp = False):
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Inputs to the code
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    #obj_sci,obj_template = 'TOI-1452', 'GL699'
    #doplot_ccf,doplot_debug,force  = False, False, False

    #lblrv_path = 'lblrv/'  # must exist, path to catalog of lines
    #mask_path = 'masks/' # path to mask files, must have _pos files in there
    #template_path = 'templates/' # stellar templates
    #science_path = 'tellurics/' # direcory with per-object telluric-corect (tcorr) files are hidden

    template = fits.getdata(template_path+'Template_s1d_'+obj_template+'_sc1d_v_file_AB.fits')

    if 'FP' not in obj_sci:
        scifiles = glob.glob(science_path+obj_sci+'/2*tcorr*AB.fits')
    else:
        scifiles = glob.glob(science_path+obj_sci+'/*C.fits')

        if check_fp:
            print('Checking if this is OBJ_FP')
            for i in tqdm(range(len(scifiles))):
                hdr = fits.getheader(scifiles[i])
                if hdr['DPRTYPE'] != 'OBJ_FP':
                    os.system('mv '+scifiles[i]+' '+science_path+obj_sci+'/others/')
                    scifiles[i] = ''
            scifiles = scifiles[scifiles != '']

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    End of inputs to the code
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Start of checks
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # reference table that will always be used afterward as a starting point
    ref_table_name = 'ref_table_{0}.csv'.format(obj_template)
    maskfile = mask_path+obj_template+'_pos.csv'


    if not os.path.isfile(ref_table_name):
        print('File {0} will be created, that''s fine but you have been warned!'.format(ref_table_name))
    else:
        print('File {0} exists, it will be read, if you want to generate it again, delete the file'.format(ref_table_name))


    if not os.path.isdir(lblrv_path):
        raise ValueError('Directory {0} does not exist, that''s really bad, we exit'.format(lblrv_path))
    else:
        print('Direcotry {0} exists, we are happy'.format(lblrv_path))

    if not os.path.isdir(mask_path):
        raise ValueError('Directory {0} does not exist, that''s really bad, we exit'.format(mask_path))
    else:
        print('Direcotry {0} exists, we are happy'.format(mask_path))

    if not os.path.isdir(template_path):
        raise ValueError('Directory {0} does not exist, that''s really bad, we exit'.format(template_path))
    else:
        print('Direcotry {0} exists, we are happy'.format(template_path))

    if not os.path.isdir(science_path):
        raise ValueError('Directory {0} does not exist, that''s really bad, we exit'.format(science_path))
    else:
        print('Direcotry {0} exists, we are happy'.format(science_path))

    if not os.path.isdir(science_path+'/'+obj_sci):
        raise ValueError('Directory {0} does not exist, that''s really bad, we exit'.format(science_path+'/'+obj_sci))
    else:
        print('Direcotry {0} exists, we are happy'.format(science_path+'/'+obj_sci))


    if not os.path.isfile(maskfile):
        raise ValueError('File {0} does not exist, that''s really bad, we exit'.format(maskfile))
    else:
        print('File {0} exists, we are happy!'.format(maskfile))

    if not os.path.isfile(ref_blaze_file):
        raise ValueError('File {0} does not exist, that''s really bad, we exit'.format(ref_blaze_file))
    else:
        print('File {0} exists, we are happy!'.format(ref_blaze_file))


    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    End of checks
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # sort files by name so that they are consecutive in time
    scifiles.sort()

    if os.path.isfile(ref_table_name) == False:
        # we create an empty table that serves as the baseline for all future observations
        # using that template as a starting point. This table is read if it exists.
        wavegrid = et.fits2wave(scifiles[0])
        tbl = Table.read(maskfile, format='ascii')

        wave_start = []
        wave_end = []
        order = []
        xpix = []
        weight_line = [] # used only to get systemic velocity as a starting point
        for i in range(wavegrid.shape[0]):
            gg = (tbl['ll_mask_s']>np.min(wavegrid[i]))*(tbl['ll_mask_s']<np.max(wavegrid[i]))
            if np.sum(gg) !=0:
                order = np.append(order,np.ones(np.sum(gg)-1)*i)
                wave_start = np.append(wave_start, np.array(tbl['ll_mask_s'][gg][:-1]))
                wave_end = np.append(wave_end, np.array(tbl['ll_mask_s'][gg][1:]))
                weight_line = np.append(weight_line,np.array(tbl['w_mask'][gg][:-1]) )
                xpix = np.append(xpix, ius(wavegrid[i],np.arange(len(wavegrid[i])))(np.array(tbl['ll_mask_s'][gg][:-1])))


        tbl = Table()
        tbl['ORDER'] = np.array(order,dtype = int)
        tbl['WAVE_START'] = wave_start #np.array(wave_start*(1 + objrv * 1000 / constants.c))
        tbl['WAVE_END'] = wave_end#np.array(wave_end* (1 + objrv * 1000 / constants.c))
        tbl['WEIGHT_LINE'] = weight_line # from mask
        tbl['XPIX'] = xpix #pixel position along array

        print('We write {0}'.format(ref_table_name))
        tbl.write(ref_table_name)
    else:
        print('We read {0}'.format(ref_table_name))
        tbl = Table.read(ref_table_name)
        order = np.array(tbl['ORDER'] )
        wave_start = np.array(tbl['WAVE_START'])
        wave_end = np.array(tbl['WAVE_END'])
        weight_line = np.array(tbl['WEIGHT_LINE'] )

    # we create the spline of the template to be used everywhere further down
    valid = np.isfinite(template['flux'])
    bl = fits.getdata(ref_blaze_file)

    if 'FP' not in obj_sci:
        # we get an estimate of the model velocity
        systemic_vel = et.get_rough_ccf_rv(template['wavelength'],template['flux'],wave_start,np.ones_like(weight_line),
                                           doplot = doplot_ccf)
    else:
        systemic_vel = 0

    # template removed from its systemic velocity, the spline is redefined for good
    spline = ius(et.doppler(template['wavelength'][valid],-systemic_vel),template['flux'][valid],k=3,ext=1)

    # we create a mask to know if the splined point  is valid
    spline_mask = ius(et.doppler(template['wavelength'],-systemic_vel),np.isfinite(template['flux']),k=1,ext=1)


    current_epoch = 1e9 # dummy epoch so we know that we did not previously have a good RV
    for ifile in range(len(scifiles)):

        # output name
        outname = lblrv_path + scifiles[ifile].split('.fits')[0].split('/')[-1]+'_'+obj_sci+'_'+obj_template+'.lblrv'


        if os.path.isfile(outname) and not force:
            print('file {0} exists and force = {1}'.format(outname, force))
            continue

        # get the science file info
        sp, hdr = fits.getdata(scifiles[ifile], header=True)
        # get the wave grid from science file
        wave = et.fits2wave(hdr)

        # if this is the first file for the day, we get a rough velocity
        if current_epoch != np.floor(hdr['MJDATE']):
            try:
                # get CCF rv as a starting point
                if 'FP' not in obj_sci:
                    rv = et.get_rough_ccf_rv(wave,sp,wave_start,np.ones_like(weight_line),doplot = doplot_ccf )
                else:
                    rv = 0

                print('For epoch {0}, the system velocity is {1:.2f} km/s'.format(np.floor(hdr['MJDATE']),rv/1000 ))
                current_epoch =  np.floor(hdr['MJDATE'])
            except:
                # shit happens
                print('we have a problem with the fitting, we just take the previous epoch value')

        # loop on the model spline velocity
        rv_zp = 0
        ite_rv = 0
        rv_mean = 1e9

        keep = np.ones_like(tbl,dtype = bool)

        bulk_error = 1e9 # will be updated with the estimated error
        while (ite_rv < 10) and (np.abs(rv_mean)>bulk_error/10.0):
            now = datetime.now()

            sp, hdr = fits.getdata(scifiles[ifile], header=True)
            wave = et.fits2wave(hdr)

            if (ite_rv == 0):
                # spline the model onto the wavelenght grid at the specified velocity
                model = np.zeros_like(sp)
                model1ms = np.zeros_like(sp)
                rms = np.zeros_like(sp)
                lowf = np.zeros_like(sp)
                model_mask = np.zeros_like(sp)


            for ii in range(sp.shape[0]):
                # RV shift the spline and give it the shape of the model
                model[ii] = spline(et.doppler(wave[ii],-rv) ) * bl[ii]
                # offset by 1 m/s to get numerical derivative
                model1ms[ii] = spline(et.doppler(wave[ii],-rv-1)) * bl[ii]

                # if on first iteration, get the low-frequency component out
                if (ite_rv == 0):
                    model_mask[ii] = spline_mask(et.doppler(wave[ii], -rv))
                    tmp = et.lowpassfilter(model[ii]/sp[ii])
                    tmp[tmp==0] = np.nan
                    lowf[ii] = tmp



                # match low frequency of model to that of spectrum
                model[ii] /= lowf[ii]
                model1ms[ii] /= lowf[ii]

                #plt.plot(model[ii])
                #plt.plot(sp[ii])
                #plt.show()

                # get residuals and running RMS
                tmp = sp[ii]-model[ii]
                # approximately the running RMS
                rms[ii] = et.lowpassfilter(np.abs(tmp))/0.664

            # keep only legit splined points
            model[model_mask < .99] = np.nan

            # create a dummy array that will contain velocities and corresponding errors
            dv = np.zeros(len(tbl['WAVE_START']))+np.nan
            dvrms = np.zeros(len(tbl['WAVE_START']))+np.nan

            # keep track of which order we are looking at
            iord_prev = -1

            # order of each line in the table
            iord = np.array(tbl['ORDER'])

            # noinspection PyInterpreter
            for i in tqdm(range(0,len(tbl)), leave = False):

                # if the line has been flagged as bad, skip right away
                if (ite_rv !=0) and (keep[i] == False):
                    continue

                # if not the same as previous order, then get residuals for that order
                if iord[i] != iord_prev:
                    iord_prev = iord[i]

                    ww_ord = et.doppler(wave[iord[i]],-rv)
                    sp_ord = sp[iord[i]]
                    rms_ord = rms[iord[i]]
                    model_ord = model[iord[i]]
                    model1ms_ord = model1ms[iord[i]]
                    dd = model_ord - model1ms_ord
                    diff = sp_ord - model_ord


                    if doplot_debug:
                        if (iord[i] == 35)*(ite_rv == 0):
                            plt.plot(ww_ord,sp_ord, color = 'grey', linewidth=6, alpha = 0.3)

                # valid points for the dot-product
                g = (ww_ord>(wave_start[i]))*(ww_ord<(wave_end[i]))*(np.isfinite(dd))*(np.isfinite(diff))

                # you need at least 5 points to do something meaningful
                if (np.sum(g)<5):
                    keep[i] = False
                    continue

                # edge points
                imin = np.min( np.min(np.where(g)))-1 # we offset by 1 to get the edge pixels
                imax = np.min( np.max(np.where(g)))+1  # we offset by 1 to get the edge pixels

                if imin<0:
                    continue
                if imax>( ww_ord.shape[0]-1):
                    continue

                # get weights at the edge of the domain. Pixels inside have a weight of 1, at the edge, it's proportional
                # to the overlap
                weight_mask = np.ones(imax-imin+1)
                weight_mask[0] = (ww_ord[imin+1]-wave_start[i])/(ww_ord[imin+1]-ww_ord[imin])
                weight_mask[-1] = (wave_end[i] - ww_ord[imax-1])/(ww_ord[imax-1]-ww_ord[imax-2])

                # maybe some plots
                if doplot_debug:
                    if (iord[i] == 35)*(ite_rv == 0):
                        color = (['red','green','blue'])[i % 3]
                        plt.plot(ww_ord[imin:imax+1],sp_ord[imin:imax+1],color = color)
                        for ii in range(imax-imin+1):
                            plt.plot(ww_ord[imin+ii],sp_ord[imin+ii] ,alpha = weight_mask[ii],color = color,marker = '.')

                # derivative of the segment
                dd_segment = dd[imin:imax+1]*weight_mask
                # residuals of the segment
                diff_segement = diff[imin:imax+1]*weight_mask


                # From Bouchy 2001 equation, RV error for each pixel
                dvrms_pix = np.sum(rms_ord[imin:imax+1]*weight_mask)/np.sum(weight_mask) / dd_segment

                # RV error for the line
                dvrms[i] = 1/np.sqrt(np.sum((1/dvrms_pix**2)))

                # feed the monster
                dv[i] = np.nansum(diff_segement*dd_segment)/np.nansum(dd_segment**2)

            if doplot_debug:
                if ite_rv ==0:
                    plt.show()

            # get the best estimate of the velocity and update spline
            rv_mean, bulk_error = et.odd_ratio_mean(dv,dvrms)
            rv -= rv_mean

            # some printouts
            print('\tdv = {0:.4f} m/s, ite_rv {1}, bulk error {2:.4f} m/s '.format(rv_mean,ite_rv,bulk_error))
            print('\tDuration of loop iteration : ', datetime.now()-now)
            ite_rv +=1

        # update the table
        if 'FP' not in obj_sci:
            tbl['RV'] = -(dv+rv-hdr['BERV']*1000) # express to have sign fine relative to convention
        else:
            tbl['RV'] = -(dv+rv) # express to have sign fine relative to convention

        tbl['DVRMS'] = dvrms

        # write the table
        print('We write file {0} {1}/{2}\n'.format(outname,ifile, len(scifiles)))

        hdu1 = fits.PrimaryHDU()

        hdr['ITE_RV'] = ite_rv,'Number of iteration to reach <10 cm/s numerical accuracy'
        hdr['SYSTVELO'] = systemic_vel, 'Systemic velocity in m/s'
        hdu1.header = hdr

        hdu2 = fits.BinTableHDU(tbl)
        new_hdul = fits.HDUList([hdu1, hdu2])
        new_hdul.writeto(outname, overwrite=True)

