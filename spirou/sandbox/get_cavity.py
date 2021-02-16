import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import etienne_tools as et
from scipy import constants
import os

def get_wave_sol(hc_lines,fp_lines, wavesol_order = 5, cavity_order = 9, nsig_cut = 5,
                 cavity_table_name = None, no_achromatic = False):

    # only keep pixels that have finite positions
    # it is fine to have orders with no valid lines
    hc_lines = hc_lines[np.isfinite(hc_lines['PIXEL_MEAS'])]
    fp_lines = fp_lines[np.isfinite(fp_lines['PIXEL_MEAS'])]

    # to get started, we assume that we do not know the wavelength of FP lines
    fp_lines['WAVE_MEAS'] = np.nan

    for iord in np.unique(fp_lines['ORDER']):
        # find hc and fp lines for the current order
        g_fp = fp_lines['ORDER'] == iord
        g_hc = hc_lines['ORDER'] == iord

        # get fp lines for this order
        fp_lines2 = fp_lines[g_fp]

        # express step in pixels as a polynomial fit. This is used to count
        # fp peaks afterward
        fit_step = \
        et.robust_polyfit(fp_lines2['PIXEL_MEAS'][1:], fp_lines2['PIXEL_MEAS'][1:] - fp_lines2['PIXEL_MEAS'][:-1], wavesol_order, nsig_cut)[0]

        # counting steps backward
        # maybe first step is wrong, we'll see later by x-matching with HC lines
        # after this step, we know that lines within the order have *relative* numbers
        # that are ok
        for i in range(1, len(fp_lines2)):
            # find expected step between previous and current FP peak
            # We start numbering at 1 as the 0th serves as a relativenstarting point
            dn = ((fp_lines2['PIXEL_MEAS'][i] - fp_lines2['PIXEL_MEAS'][i - 1]) / np.polyval(fit_step,
                                                                                             fp_lines2['PIXEL_MEAS'][
                                                                                                 i - 1]))
            # dn is always very close to an integer value, we round it
            # we subtract the steps, FP peaks go in decreasing number
            fp_lines2['PEAK_NUMBER'][i] = fp_lines2['PEAK_NUMBER'][i - 1] - np.round(dn)

        # if we have more than 5 hc lines, we try to find the absolute numbering
        if np.sum(g_hc)>5:
            hc_lines2 = hc_lines[g_hc]
            # we fit an approximate wavelength solution
            wave_fit = et.robust_polyfit(hc_lines2['PIXEL_MEAS'],hc_lines2['WAVE_REF'],wavesol_order,nsig_cut)[0]

            # we find the steps in FP lines at the position of all HC lines
            step_hc = np.polyval(fit_step,hc_lines2['PIXEL_MEAS'])
            step_wave = np.polyval(np.polyder(wave_fit),hc_lines2['PIXEL_MEAS'])*step_hc

            # convert step in cavity through the order. We assume a constant cavity
            # through the order
            cavity_per_order = np.nanmedian(hc_lines2['WAVE_REF'] ** 2 / step_wave)
            cavity_per_order0 = np.array(cavity_per_order)

            # we explore integer offset in FP peak counting and find the offset that
            # gives the defines the wavelength solution leading to the smallest RMS
            # between predicted and catalog HC peak positions
            offset = np.arange(-5,6)
            sig = np.zeros_like(offset,dtype = float)
            for i in range(len(offset)):
                cavity_per_order = np.array(cavity_per_order0)
                PEAK_NUMBER_OFFSET = fp_lines2['PEAK_NUMBER']+offset[i]

                for ite in range(3):
                    wave_tmp = cavity_per_order/PEAK_NUMBER_OFFSET
                    wave_fit = np.polyfit(fp_lines2['PIXEL_MEAS'],wave_tmp,wavesol_order)
                    med = np.nanmedian(1-hc_lines2['WAVE_REF']/np.polyval(wave_fit,hc_lines2['PIXEL_MEAS']))
                    cavity_per_order*=(1-med)

                sig[i] = et.sigma(1-hc_lines2['WAVE_REF']/np.polyval(wave_fit,hc_lines2['PIXEL_MEAS']))*constants.c
            # we apply the offset that leads to the smallest HC (o-c) RMS
            fp_lines2['PEAK_NUMBER']=   fp_lines2['PEAK_NUMBER']+ offset[np.argmin(sig)]

            # we find the best cavity length estimate
            for ite in range(3):
                # we force the cavity length to lead to a median HC peak position error of zero.
                # we could use a better sigma-clipping, but this is hard with a small number of lines
                wave_tmp = cavity_per_order / fp_lines2['PEAK_NUMBER']
                wave_fit = np.polyfit(fp_lines2['PIXEL_MEAS'], wave_tmp, wavesol_order)
                med = np.nanmedian(1 - hc_lines2['WAVE_REF'] / np.polyval(wave_fit, hc_lines2['PIXEL_MEAS']))
                cavity_per_order *= (1 - med)

            # we now have a best-guess of the wavelength solution, we update
            #  the WAVE_MEAS in the FP line list. This will be used to constrain the
            # cavity length below
            wave_fit = np.polyfit(fp_lines2['PIXEL_MEAS'], wave_tmp, wavesol_order)
            fp_lines2['WAVE_MEAS'] = np.polyval(wave_fit,fp_lines2['PIXEL_MEAS'])

        # put into the table. If we had enough HC lines, the WAVE_MEAS has been updated
        # if not, at least the FP peak counting is valid.
        fp_lines[g_fp] = fp_lines2

    # plot cavity length prior to finding offsets
    fig, ax = plt.subplots(nrows = 2, ncols = 1,sharex = True)
    ax[0].plot(fp_lines['WAVE_MEAS'],fp_lines['WAVE_MEAS']*fp_lines['PEAK_NUMBER'],'r.',alpha = 0.4,
               label = 'cavity prior to glitch corr.')


    for ite in range(2): # do twice just to be sure it's all fine
        for iord in np.unique(fp_lines['ORDER'])[1:]:
            # skip first order, and check order to order if the cavity
            # is consistent with previous. Order=1 is compared to Order=0, then
            # Order=2 to Order=1...

            g_prev = fp_lines['ORDER'] == (iord-1)
            g = fp_lines['ORDER'] == iord

            # current peak numbers
            current_numbering = fp_lines['PEAK_NUMBER'][g]
            extrapolated_numbering = np.nanmedian(fp_lines['WAVE_MEAS'][g_prev]*fp_lines['PEAK_NUMBER'][g_prev])/\
                                     fp_lines['WAVE_MEAS'][g]

            # current peak numbers if you take the previous order cavity length and
            # assume it's the same of current order while using the current order wavelength
            # solution.
            # The median between the two should be close to zero (extrapolated_numbering can
            # be float, not int as current_numbering). If the is a mis-counting, then you get
            # a median very close to an integer value. This tells you the offset to be applied
            # to have a consistent cavity lenth
            offset = np.nanmedian(current_numbering-extrapolated_numbering)
            if  np.isfinite(offset):
                fp_lines['PEAK_NUMBER'][g] -= int(np.round(offset))

        # now that we have patched the order-to-order glitches, we look for a bulk offset
        # in the counts that minimizes the dispersion in cavity length.
        global_offset = np.arange(-10,11)
        sig = np.zeros_like(global_offset,dtype = float)
        for i in range(len(global_offset)):
            cavity = fp_lines['WAVE_MEAS']*(fp_lines['PEAK_NUMBER']+global_offset[i])
            sig[i] = et.sigma(cavity)
        fp_lines['PEAK_NUMBER'] = fp_lines['PEAK_NUMBER'] +global_offset[np.argmin(sig)]

    # now we have valid numbering and best-guess WAVE_MEAS, we find the cavity length
    cavity = et.robust_polyfit(fp_lines['WAVE_MEAS'], fp_lines['WAVE_MEAS']*fp_lines['PEAK_NUMBER'],cavity_order,nsig_cut)[0]

    ax[0].plot(fp_lines['WAVE_MEAS'],fp_lines['WAVE_MEAS']*fp_lines['PEAK_NUMBER'],'g.',alpha = 0.4,
               label = 'cavity after glitch corr.')
    ax[1].plot(fp_lines['WAVE_MEAS'],fp_lines['WAVE_MEAS']*fp_lines['PEAK_NUMBER']-np.polyval(cavity,fp_lines['WAVE_REF']),'k.',alpha= 0.5)
    ax[0].set(xlabel = 'Wavelength [nm]',ylabel = 'Cavity [nm]')
    ax[1].set(xlabel = 'Wavelength [nm]',ylabel = 'Cavity [nm]')
    ax[0].legend()
    plt.tight_layout()

    plt.show()

    # if we gave a table name and the file is absent, we write. If the file is
    # present, then we'll only correct the achromatic term (no_achromatic == False)

    do_write = False
    if cavity_table_name != None:
        if os.path.isfile(cavity_table_name):
            cavity = np.array(Table.read(cavity_table_name)['coeffs'])
        else:
            do_write = True
    else:
        # if no file and no name, we'll save
        cavity_table_name = 'cavity.csv'
        do_write = True


    cavity0 = np.array(cavity)
    mean2err=np.inf

    # we change the achromatic cavity length term to force HC peaks to have a zero velocity error.

    # we will loop until the cavity change is below 1/100th of the error on the cavity change.
    while mean2err>1e-2:
        # get the proper cavity length from the cavity polynomial
        for ite in range(3):
            fp_lines['WAVE_REF'] = np.polyval(cavity,fp_lines['WAVE_REF'])/fp_lines['PEAK_NUMBER']

        # get the wavelength solution for the order and the HC line position that it implies.
        # the diff between the HC position found here and the catalog one is used to
        # change the cavity length

        for iord in np.unique(fp_lines['ORDER']):
            g_fp = fp_lines['ORDER'] == iord
            wave_fit = et.robust_polyfit(fp_lines['PIXEL_MEAS'][g_fp],fp_lines['WAVE_REF'][g_fp],wavesol_order,nsig_cut)[0]
            fp_lines['WAVE_MEAS'][g_fp] = np.polyval(wave_fit,fp_lines['PIXEL_MEAS'][g_fp])

            g_hc = hc_lines['ORDER'] == iord
            if np.sum(g_hc) !=0:
                hc_lines['WAVE_MEAS'][g_hc] = np.polyval(wave_fit, hc_lines['PIXEL_MEAS'][g_hc])

        # in velocity, diff between measured and catalog HC line positions
        diff_hc = (1-hc_lines['WAVE_MEAS']/hc_lines['WAVE_REF'])*constants.c

        # model of the errors in the HC line positions. We assume that
        # they decrease as 1/NSIG
        sig = et.sigma(diff_hc*(hc_lines['NSIG']))/(hc_lines['NSIG'])

        # get smart mean of the velocity error
        mean_hc_vel,err = et.odd_ratio_mean(diff_hc,sig,odd_ratio=1e-2)

        # if we are allowed to change the achromatic cavity length, then we do it, else
        # we just keep track of how much we would have changed it.
        if no_achromatic == False:
            mean2err = np.abs(mean_hc_vel/err)
            cavity[-1]*=(1+mean_hc_vel/constants.c)
        else:
            mean2err = 0

        # mean error in position of HC lines
        print('Mean HC position {0:6.2f}+-{1:.2f} m/s'.format(mean_hc_vel,err))

    # final cavity change
    print('change in cavity length {0:6.2f} nm'.format(cavity[-1]-cavity0[-1]))

    #
    fig,ax = plt.subplots(nrows = 2, ncols = 1)
    ax[0].hist(diff_hc[np.abs(diff_hc)<500],bins = 100)
    ax[0].set(xlabel = 'Velocity [m/s]',title = 'HC diff')
    nsig = diff_hc/sig
    ax[1].hist(nsig[np.abs(nsig)<5],bins = 100)
    ax[1].set(xlabel = 'sigma',title = 'HC diff')
    plt.tight_layout()
    plt.show()

    # we update WAVE_REF
    for ite in range(3):
        fp_lines['WAVE_REF'] = np.polyval(cavity, fp_lines['WAVE_REF']) / fp_lines['PEAK_NUMBER']

    # we construct the output matrix to keep track of wavelength coefficients
    wave_fit = np.zeros([len(np.unique(fp_lines['ORDER'])),wavesol_order+1])
    for iord in np.unique(fp_lines['ORDER']):
        g_fp = fp_lines['ORDER'] == iord
        wave_fit[iord] = et.robust_polyfit(fp_lines['PIXEL_MEAS'][g_fp], fp_lines['WAVE_REF'][g_fp], wavesol_order, 5)[0]

    # if required, we save the cavity file
    if do_write:
        print('We write {}'.format(cavity_table_name))
        tbl  =Table()
        tbl['coeffs'] = cavity
        tbl.write(cavity_table_name,overwrite = True)

    return wave_fit

if os.path.isfile('sample_cavity.csv'):
    os.system('rm sample_cavity.csv')

# we create the cavity file
fp_lines = Table.read('EA_23994F38T42a_pp_e2dsff_AB_wavem_fplines_AB.fits')
hc_lines = Table.read('EA_239944F5T6c_pp_e2dsff_AB_wavem_hclines_AB.fits')

wave_fit = get_wave_sol(hc_lines,fp_lines,cavity_table_name = 'sample_cavity.csv')
print()

# we read the cavity file and perform an achromatic cavity change
fp_lines = Table.read('EA_242680F5T9a_pp_e2dsff_AB_wavem_fplines_AB.fits')
hc_lines = Table.read('EA_241318F5T6c_pp_e2dsff_AB_wavem_hclines_AB.fits')

wave_fit = get_wave_sol(hc_lines,fp_lines,cavity_table_name = 'sample_cavity.csv')
print()

# we read the cavity file but perform do *NOT* perform the achromatic cavity change
fp_lines = Table.read('EA_24004F09T13a_pp_e2dsff_AB_wavem_fplines_AB.fits')
hc_lines = Table.read('EA_240041F6T7c_pp_e2dsff_AB_wavem_hclines_AB.fits')
wave_fit = get_wave_sol(hc_lines,fp_lines,cavity_table_name = 'sample_cavity.csv',
                        no_achromatic=True)
print()
