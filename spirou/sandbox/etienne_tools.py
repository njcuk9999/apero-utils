import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.optimize import curve_fit

def doppler(wave,v):
    # velocity expressed in m/s
    # relativistic calculation
    return wave*np.sqrt(  (1-v/constants.c)/(1+v/constants.c))

def weighted_median(values, weights):
    keep = np.isfinite(values)*np.isfinite(weights)
    values1 = np.array(values[keep], dtype = float)
    weights1 = np.array(weights[keep], dtype = float)
    weights1 /= np.nansum(weights1)

    ord = np.argsort(values1)
    values1 = values1[ord]
    weights1 = weights1[ord]

    cumsum = np.cumsum(weights1)
    imed = np.min(np.where(cumsum>0.5))

    return values1[imed]

def get_rough_ccf_rv(wave,sp,wave_mask,weight_line, doplot = False):

    if len(wave.shape) == 2: # we have the e2ds file, we reshape it
        mask = np.ones_like(wave,dtype = bool)

        for iord in range(1,49):
            mask[iord]*=(wave[iord-1,::-1]<wave[iord])

        for iord in range(0,48):
            mask[iord]*=(wave[iord]<wave[iord+1,::-1])

        mask*=np.isfinite(sp)

        sp2 = sp[mask]
        wave2 = wave[mask]
    else:# we have the s1d, we just need to get rid of NaNs
        sp2 = sp[np.isfinite(sp)]
        wave2 = wave[np.isfinite(sp)]

    spline_sp = ius(wave2,sp2, k=1, ext=1)

    dvs = np.arange(-1.5e5, 1.5e5, 500)
    ccf = np.zeros_like(dvs)
    print('computing CCF')
    for i in range(len(dvs)):
        ccf[i] = np.nansum(weight_line*spline_sp(doppler(wave_mask, dvs[i])))

    imax = np.argmax(ccf)

    guess = [dvs[imax],2000,ccf[imax]-np.nanmedian(ccf),np.nanmedian(ccf),0]
    fit, pcov = curve_fit(gauss, dvs,ccf,p0 = guess)

    #fit = np.polyfit(dvs[imax - 1:imax + 2], ccf[imax - 1:imax + 2], 2)
    systemic_vel = fit[0]#-.5 * fit[1] / fit[0]
    print('CCF Velocity : {0:.2f} m/s'.format(-systemic_vel))


    if doplot:
        plt.plot(-dvs/1000, ccf )
        plt.plot(-dvs/1000,   gauss(dvs,*fit))
        plt.xlabel('RV [km/s]')
        plt.ylabel('Normalized CCF')
        plt.show()

    return systemic_vel

def odd_ratio_mean(value,err, odd_ratio = 1e-4, nmax = 10):
    #
    # Provide values and corresponding errors and compute a
    # weighted mean
    #
    #
    # odd_bad -> probability that the point is bad
    #
    # nmax -> number of iterations

    guess = np.nanmedian(value)

    nite = 0
    while (nite < nmax):
        nsig = (value-guess)/err
        gg = np.exp(-0.5*nsig**2)
        odd_bad = odd_ratio/(gg+odd_ratio)
        odd_good = 1-odd_bad

        w = odd_good/err**2

        guess = np.nansum(value*w)/np.nansum(w)
        nite+=1

    bulk_error =  np.sqrt(1/np.nansum(odd_good/err**2))

    return guess,bulk_error

def sigma(tmp):
    # return a robust estimate of 1 sigma
    sig1 = 0.682689492137086
    p1 = (1-(1-sig1)/2)*100
    return (np.nanpercentile(tmp,p1) -np.nanpercentile(tmp,100-p1))/2.0

def gauss(x,cen, ew, amp, zp, slope):
    return np.exp(-0.5*(x-cen)**2/ew**2)*amp+zp+(x-cen)*slope
