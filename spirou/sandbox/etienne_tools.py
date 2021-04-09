import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from scipy.optimize import curve_fit
from astropy.io import fits
from astropy.table import Table
import warnings
#import numba
from numba import jit
from scipy.interpolate import InterpolatedUnivariateSpline as ius
import requests
import os

def mk_hash_name(hash_name,path,suffix = '_pp.fits',check_exist = True):
    # create a file list from Neil's hash
    hash_name = hash_name.split('_')[0]
    file_type = hash_name[-1]
    hash_name=hash_name[:-1]


    if 'T' in hash_name:
        prefix = hash_name.split('F')[0]
        range = hash_name.split('F')[1].split('T')
        index = np.arange(int(range[0]), int(range[1]) + 1)
        files = []
        for i in index:
            files = np.append(files, path +prefix+ str(i).zfill(len(range[1])) + suffix)
    else:
        files = np.array([path+hash_name+file_type+suffix]) # we just have one file here

    if check_exist:
        for i in np.arange(len(files)):
            if os.path.isfile(files[i]) == False:
                files[i] = ''
        if  np.sum(files != '') == 0:
            return None

        else:
            return files[files != '']

    else:
        return files

def get_xycen(im,x0,y0,w=0):

    x = np.array(x0+0.5,dtype = int)
    y = np.array(y0+0.5,dtype = int)
    w = int(w)

    if w !=0:
        col1,col2,col3 = np.zeros_like(x0,dtype = float),np.zeros_like(x0,dtype = float),np.zeros_like(x0,dtype = float)
        for ww in range(-w,w+1):
            col1 += im[x+ww,y-1]
            col2 += im[x+ww,y]
            col3 += im[x+ww,y+1]
    else:
        col1 = im[x,y-1]
        col2 = im[x,y]
        col3 = im[x,y+1]

    ycen =  .5 * (col1 - col3) / (col1 + col3 - 2 * col2) + y

    if w !=0:
        col1,col2,col3 = np.zeros_like(x0,dtype = float),np.zeros_like(x0,dtype = float),np.zeros_like(x0,dtype = float)
        for ww in range(-w,w+1):
            col1 += im[x - 1, y+ww]
            col2 += im[x, y+ww]
            col3 += im[x + 1, y+ww]

    else:
        col1 = im[x - 1, y]
        col2 = im[x, y]
        col3 = im[x + 1, y]

    xcen = .5 * (col1 - col3) / (col1 + col3 - 2 * col2) + x

    return xcen,ycen

def harps2spirou(hdr):
    # input the header of a HARPS image and add the corresponding SPIRou keywords.
    # the code leaves in place the HARPS keywords.

    if 'COMMENT' in hdr:
        del hdr['COMMENT']

    badkey = ['ESO INS TEMP1 RMS','ESO INS TEMP40 RMS']

    for key in badkey:
        if key in hdr:
            del hdr[key]


    hdr2 = dict(hdr)
    for key in hdr.keys():
        # val =  hdr[key]
        # del hdr[key]
        val = hdr2[key]
        if type(val) not in [int, bool, float, str]:
            hdr[key] = str(val)
            print(et.color('Value {0} had type {1}, now forced to str'.format(val, type(val)), 'red'))
        else:
            hdr[key] = val


    keys = [
        ['MJDATE','HIERARCH ESO DRS BJD'],
        ['MJDMID','HIERARCH ESO DRS BJD'],
        ['EXPTIME','HIERARCH ESO DET WIN1 DIT1'],
        ['AIRMASS','HIERARCH ESO TEL AIRM START'],
        ['FILENAME',''],
        ['DATE-OBS','HIERARCH ESO INS DATE'],
        ['BERV','HIERARCH ESO DRS BERV'],
        ['TAU_H2O',''],
        ['TAU_OTHE',''],
        ['ITE_RV',''],
        ['SYSTVELO',''],
        ['WAVETIME',''],
        ['WAVEFILE',''],
        ['TLPDVH2O',''],
        ['TLPDVOTR',''],
        ['CDBWAVE',''],
        ['OBJECT','HIERARCH ESO OBS TARG NAME'],
        ['SBRHB1_P',''],
        ['SBRHB2_P',''],
        ['SBCDEN_P',''],
        ['SNRGOAL',''],
        ['EXTSN035',''],
        ['BJD','HIERARCH ESO DRS BJD'],
        ['SHAPE_DX', ''],
        ['SHAPE_DY',''],
        ['SHAPE_A',''],
        ['SHAPE_B',''],
        ['SHAPE_C',''],
        ['SHAPE_D','']]
    keys = np.array(keys)

    for i in range(keys.shape[0]):
        if keys[i,1]!='':
            hdr[keys[i, 0]] = hdr[keys[i,1]]
        else:
            hdr[keys[i, 0]] = 0.0

    # we get the wavelength solution in the HARPS header
    for i in range(0,999):
        key = 'HIERARCH ESO DRS CAL TH COEFF LL'+str(i)
        if key in hdr:
            hdr['WAVE'+str(i).zfill(4)] = hdr[key]/10.0 # expressed in nm, not Angstrom for SPIRou

    if 'NAXIS2' in hdr:
        nord = hdr['NAXIS2']
        # that's to match spirou's representation
        hdr['WAVEORDN'] = nord

    hdr['EXTSN035'] = 100 # we fake a good observation as there is no equivalent in HARPS headers

    for key in hdr.keys():
        try:
            _= hdr[key]
        except:
            hdr[key] = ''

    return hdr

# magic grid is a standard (well, my standard) way of representing a wavelength vector
# it is set so that each element is exactly dv_grid step in velocity. If you shift
# your velocity, then you have a simple translation of this vector.
def get_magic_grid(wave0=1500,wave1=1800,dv_grid=0.5):
    # default for the function is 500 m/s
    # the arithmetic is a but confusing here, you first find how many
    # elements you have on your grid, then pass it to an exponential
    # the first element is exactely wave0, the last element is NOT
    # exactly wave1, but is very close and is set to get your exact
    # step in velocity
    len_magic = int(np.ceil(np.log(wave1/wave0)*np.array(constants.c)/dv_grid))
    magic_grid =  np.exp(np.arange(len_magic)/len_magic*np.log(wave1/wave0))*wave0
    return magic_grid


def get_bad_odo():

    URL_BASE = ('https://docs.google.com/spreadsheets/d/'
                '{}/gviz/tq?tqx=out:csv&sheet={}')
    SHEET_ID = '1gvMp1nHmEcKCUpxsTxkx-5m115mLuQIGHhxJCyVoZCM'
    WORKSHEET = 0
    BAD_ODO_URL = URL_BASE.format(SHEET_ID, WORKSHEET)


    # fetch data
    data = requests.get(BAD_ODO_URL)
    tbl = Table.read(data.text, format='ascii')
    # Convert types
    tbl['ODOMETER'] = tbl['ODOMETER'].astype(str)
    tbl['PP'] = tbl['PP'] == 'TRUE'
    tbl['RV'] = tbl['RV'] == 'TRUE'

    return np.array(tbl['ODOMETER'])


def nanpercentile(v,p,axis = None):
    if axis == None:
        return jit_nanpercentile(v, p)
    else:
        return np.nanpercentile(v,p,axis = axis)

@jit(nopython=True)
def jit_nanpercentile(v, p):
    return np.nanpercentile(v,p)

@jit(nopython=True)
def nanstd(v):
    return np.nanstd(v)

@jit(nopython=True)
def nansum(v):
    return np.nansum(v)

@jit(nopython=True)
def sum(v):
    return np.sum(v)

@jit(nopython=True)
def mean(v):
    return np.mean(v)

@jit(nopython=True)
def std(v):
    return np.std(v)

def nanmean(v):
    with warnings.catch_warnings(record=True) as _:
        mean = jitnanmean(v)
    return mean

@jit(nopython=True)
def jitnanmean(v):
    mean = np.nanmean(v)
    return mean

#@jit(nopython=True)
def nanmedian(v):
    med = np.nanmedian(v)
    #g = np.isfinite(v)
    return med

@jit(nopython=True)
def median(v):
    return np.median(v)

@jit(nopython=True)
def exp(v):
    return np.exp(v)


def sigma(im):
    return (nanpercentile(im,85) - nanpercentile(im,15))/2

def running_sigma(v,w):
    # provide a vector and return a running robust dispersion
    # the width (w) is the box size to measure the dispersion
    # pick a w that is wide enough to be representative of the
    # noise but small enough to avoid variations
    ll = len(v)

    sigma = np.zeros(ll)
    for i in range(ll):
        i0 = ll-w//2
        i1 = ll+w//2
        if i0 < 0:
            i0 = 0
        if i1 > ll:
            i1 = ll
        sigma[i] = sigma(v[i0:i1])

    return sigma

def color(message,color):
    COLOURS = dict()
    COLOURS['BLACK'] = '\033[90;1m'
    COLOURS['RED'] = '\033[1;91;1m'
    COLOURS['GREEN'] = '\033[92;1m'
    COLOURS['YELLOW'] = '\033[1;93;1m'
    COLOURS['BLUE'] = '\033[94;1m'
    COLOURS['MAGENTA'] = '\033[1;95;1m'
    COLOURS['ORANGE']  = '\033[1;99;1m'
    COLOURS['CYAN'] = '\033[1;96;1m'
    COLOURS['WHITE'] = '\033[97;1m'
    COLOURS['ENDC'] = '\033[0;0m'

    return COLOURS[color.upper()] + message + COLOURS['ENDC']


def get_ratio(sp1,sp2):
    # get the scaling ratio that minimizes least-square between two spectra with a

    #iord = 40
    #sp1 = np.array(sp)
    #sp2 = np.array(model)

    g = np.isfinite(sp1) * np.isfinite(sp2) * (np.abs(sp1)<5*sigma(sp1)) * (np.abs(sp2)<5*sigma(sp2)) #* (np.abs(diff / et.sigma(diff)) < 10)

    amp = np.sqrt(np.nansum(sp1[g]**2)/np.nansum(sp2[g]**2))

    #sp1 -= np.nanmedian(sp1[g])
    #sp2 -= np.nanmedian(sp2[g])
    for ite in range(5):
        diff = sp1-amp*sp2

        g = np.isfinite(sp1) * np.isfinite(sp2) * (np.abs(diff/sigma(diff))<3)

        v=np.nansum(diff[g]*sp2[g])/np.sqrt(np.nansum(sp2[g]**2)*np.nansum(sp1[g]**2))

        amp/=(1-v)
        #print(amp,v)

    #sp1[~g] = np.nan
    #sp2[~g] = np.nan
    #plt.plot(sp1)
    #plt.plot(sp2*amp,alpha = 0.5)
    #plt.plot(sp1-sp2*amp,color = 'red',alpha = 0.5)
    #plt.show()

    return amp

def doppler_shift(wave,sp,v,k=3):

    # apply doppler shift to data using spline and Doppler function
    # you may adjust k to avoid avoid spikes at the edge of the valid
    # domain. k=1 gives a less accurate spline but not spike.

    wave2 = doppler(wave,v)

    sp2 = np.zeros_like(sp)+np.nan

    for iord in np.arange(sp.shape[0]):
        g = np.isfinite(sp[iord])
        if np.sum(g)<5:
            continue
        spl = ius(wave[iord][g],sp[iord][g],k=k,ext=1)
        spl_mask = ius(wave[iord],np.array(g,dtype = float),k=1,ext=1)

        # keep only points having a contribution from valid pixels>0.99
        valid = spl_mask(wave2[iord]) > 0.99
        sp2[iord][valid] = spl(wave2[iord])[valid]

    return sp2


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

def wave2wave(e2ds_data_input, wave1, wave2):
    # transform e2ds data from one wavelength grid to another.

    e2ds_data = np.array(e2ds_data_input)

    for iord in range(49):
        keep = np.isfinite(e2ds_data[iord])
        spl = ius(wave1[iord][keep],e2ds_data[iord][keep],k=3, ext=1)
        e2ds_data[iord][keep] = spl(wave2[iord][keep])

    return e2ds_data

def fit_gauss(x,y,p0):
    fit, pcov = curve_fit(gauss, x,y,p0 = p0)
    return fit


def fit_super_gauss(x,y,p0):
    # values:
    # cen, ew, amp,expo, zp, slope
    fit, pcov = curve_fit(super_gauss, x,y,p0 = p0)
    return fit


def get_rough_ccf_rv(wave,sp,wave_mask,weight_line, doplot = False):

    if len(wave.shape) == 2: # we have the e2ds file, we reshape it
        mask = np.ones_like(wave,dtype = bool)


        for iord in range(1,wave.shape[0]):
            mask[iord]*=(wave[iord-1,::-1]<wave[iord])

        for iord in range(0,wave.shape[0]-1):
            mask[iord]*=(wave[iord]<wave[iord+1,::-1])

        mask*=np.isfinite(sp)

        sp2 = sp[mask]
        wave2 = wave[mask]
    else:# we have the s1d, we just need to get rid of NaNs
        sp2 = sp[np.isfinite(sp)]
        wave2 = wave[np.isfinite(sp)]

    spline_sp = ius(wave2,sp2, k=1, ext=1)

    dvs = np.arange(-3e5, 3e5, 500)
    ccf = np.zeros_like(dvs)
    print('computing CCF')
    for i in range(len(dvs)):
        ccf[i] = np.nansum(weight_line*spline_sp(doppler(wave_mask, dvs[i])))

    imax = np.argmax(ccf)

    guess = [dvs[imax],2000,ccf[imax]-nanmedian(ccf),nanmedian(ccf),0]
    fit, pcov = curve_fit(gauss, dvs,ccf,p0 = guess)

    systemic_vel = fit[0]
    print('CCF Velocity : {0:.2f} m/s'.format(-systemic_vel))

    if doplot:
        plt.plot(-dvs/1000, ccf )
        plt.plot(-dvs/1000,   gauss(dvs,*fit))
        plt.xlabel('RV [km/s]')
        plt.ylabel('Normalized CCF')
        plt.show()

    # returns the systemic velocity and e width of the CCF
    return systemic_vel, fit[1]

def robust_polyfit(x, y, degree, nsigcut):
    keep = np.isfinite(y)
    # set the nsigmax to infinite
    nsigmax = np.inf
    # set the fit as unset at first
    fit = None
    # while sigma is greater than sigma cut keep fitting
    while nsigmax > nsigcut:
        # calculate the polynomial fit (of the non-NaNs)
        fit = np.polyfit(x[keep], y[keep], degree)
        # calculate the residuals of the polynomial fit
        res = y - np.polyval(fit, x)
        # work out the new sigma values
        sig = nanmedian(np.abs(res))
        if sig == 0:
            nsig = np.zeros_like(res)
            nsig[res != 0] = np.inf
        else:
            nsig = np.abs(res) / sig
        # work out the maximum sigma
        nsigmax = np.max(nsig[keep])
        # re-work out the keep criteria
        keep = nsig < nsigcut
    # return the fit and the mask of good values
    return fit, keep


def sigma(tmp):
    # return a robust estimate of 1 sigma
    sig1 = 0.682689492137086
    p1 = (1-(1-sig1)/2)*100
    return (np.nanpercentile(tmp,p1) -np.nanpercentile(tmp,100-p1))/2.0

def gauss(x,cen, ew, amp, zp, slope):
    return np.exp(-0.5*(x-cen)**2/ew**2)*amp+zp+(x-cen)*slope

def super_gauss(x,cen, ew, amp,expo, zp, slope):
    return np.exp(-0.5*(np.abs(x-cen)/ew)**np.abs(expo))*amp+zp+(x-cen)*slope


def lowpassfilter(input_vect,width = 101):
    # Computes a low-pass filter of an input vector. This is done while properly handling
    # NaN values, but at the same time being reasonably fast.
    # Algorithm:
    #
    # provide an input vector of an arbtrary length and compute a running NaN median over a
    # box of a given length (width value). The running median is NOT computed at every pixel
    # but at steps of 1/4th of the width value. This provides a vector of points where
    # the nan-median has been computed (ymed) and mean position along the input vector (xmed)
    # of valid (non-NaN) pixels. This xmed/ymed combination is then used in a spline to
    # recover a vector for all pixel positions within the input vector.
    #
    # When there are no valid pixel in a 'width' domain, the value is skipped in the creation
    # of xmed and ymed, and the domain is splined over.

    # indices along input vector
    index = np.arange(len(input_vect))

    # placeholders for x and y position along vector
    xmed = []
    ymed = []

    # loop through the lenght of the input vector
    for i in np.arange(-width//2,len(input_vect)+width//2,width//4):

        # if we are at the start or end of vector, we go 'off the edge' and
        # define a box that goes beyond it. It will lead to an effectively
        # smaller 'width' value, but will provide a consistent result at edges.
        low_bound = i
        high_bound = i+int(width)

        if low_bound<0:
            low_bound = 0
        if high_bound>(len(input_vect)-1):
            high_bound = (len(input_vect)-1)

        pixval = index[low_bound:high_bound]

        if len(pixval)<3:
            continue

        # if no finite value, skip
        if np.sum(np.isfinite(input_vect[pixval])) < 3:
            continue

        # mean position along vector and NaN median value of
        # points at those positions
        xmed.append(nanmean(pixval))
        ymed.append(nanmedian(input_vect[pixval]))

    xmed = np.array(xmed,dtype = float)
    ymed = np.array(ymed,dtype = float)

    # we need at least 3 valid points to return a
    # low-passed vector.
    if len(xmed) < 3:
        return np.zeros_like(input_vect)+np.nan

    if len(xmed) != len(np.unique(xmed)):
        xmed2 = np.unique(xmed)
        ymed2 = np.zeros_like(xmed2)
        for i in range(len(xmed2)):
            ymed2[i] = np.mean(ymed[xmed == xmed2[i]])
        xmed = xmed2
        ymed = ymed2

    # splining the vector
    spline = ius(xmed,ymed, k=2, ext=3)
    lowpass = spline(np.arange(len(input_vect)))

    return lowpass

def air_index(wavelength, t=15., p=760.,Unit='nm'):

    #
    # Gracieusté de Romain Allart qui l'a eue d'on ne sait où
    #
    # for harps, transform the wavelength in vacuum on the air
    # wavelength in nm, t=temperature en C, p=pression en millibar
    #

    if Unit=='Ang':
        wavelength=wavelength/10.
    elif Unit=='nm':
        wavelength=wavelength
    else:
        print('Error: wrong units to convert from vaccum to air')
    n = 1e-6 * p * (1 + (1.049-0.0157*t)*1e-6*p) / 720.883 / (1 + 0.003661*t) * \
        (64.328 + 29498.1/(146-(1e3/wavelength)**2) + 255.4/(41-(1e3/wavelength)**2))
    n = n + 1

    return n


def fits2wave(file_or_header):
    info = """
        Provide a fits header or a fits file
        and get the corresponding wavelength
        grid from the header.

        Usage :
          wave = fits2wave(hdr)
                  or
          wave = fits2wave('my_e2ds.fits')

        Output has the same size as the input
        grid. This is derived from NAXIS 
        values in the header
    """

    # check that we have either a fits file or an astropy header
    if type(file_or_header) == str:
        hdr = fits.getheader(file_or_header)
    elif str(type(file_or_header)) == "<class 'astropy.io.fits.header.Header'>":
        hdr = file_or_header
    else:
        print()
        print('~~~~ wrong type of input ~~~~')
        print()

        print(info)
        return []

    if hdr['INSTRUME'] == 'HARPS':
        hdr = harps2spirou(hdr) # make the header SPIRou friendly

    # get the keys with the wavelength polynomials
    wave_hdr = hdr['WAVE0*']
    # concatenate into a numpy array
    wave_poly = np.array([wave_hdr[i] for i in range(len(wave_hdr))])

    # get the number of orders
    nord = hdr['WAVEORDN']

    # get the per-order wavelength solution
    wave_poly = wave_poly.reshape(nord, len(wave_poly) // nord)

    # get the length of each order (normally that's 4088 pix)
    npix = hdr['NAXIS1']

    # project polynomial coefficiels
    wavesol = [np.polyval(wave_poly[i][::-1], np.arange(npix)) for i in range(nord)]

    # return wave grid
    return np.array(wavesol)

def smart_time(time_in_s):
    if time_in_s>3600:
        h = str(int(np.floor(time_in_s/3600)))+'h'
        time_in_s -= (np.floor(time_in_s/3600)*3600)
        flag_h = True
    else:
        h = ''
        flag_h = False

    if time_in_s>60:
        min = str(int(np.floor(time_in_s/60))).zfill(2)+'m'
        time_in_s -= (np.floor(time_in_s/60)*60)
    else:
        min = ''


    if flag_h:
        return h+min
    else:
        sec = str(int(np.round(time_in_s))).zfill(2) + 's'
        return min+sec


def td_convert(instance):
    if isinstance(instance, Table):
        out = dict()
        for col in instance.keys():
            out[col] = np.array(instance[col])
        return out
    if isinstance(instance, dict):
        out = Table()
        for col in instance.keys():
            out[col] = np.array(instance[col])
        return out

def running_rms(sp1):
    sp1b = np.zeros(4096)+np.nan
    sp1b[4:-4] = sp1
    with warnings.catch_warnings(record=True) as _:
        b1 = np.nanpercentile(np.reshape(sp1b, [16, 256]),[16,84], axis=1)
    rms =  (b1[1]-b1[0])/2
    index = np.arange(16)*256+128
    keep = np.isfinite(rms)
    index = index[keep]
    rms = rms[keep]

    return ius(index,rms,k=2,ext=3)(np.arange(len(sp1)))

def sed_ratio(sp1,sp2,doplot = False):

    sp1b = np.zeros(4096)+np.nan
    sp2b = np.zeros(4096)+np.nan
    sp1b[4:-4] = sp1
    sp2b[4:-4] = sp2

    invalid = (np.isfinite(sp1b)*np.isfinite(sp2b)) == False
    sp1b[invalid] = np.nan
    sp2b[invalid] = np.nan

    index = np.arange(128)*32+16
    b1 = np.nansum(np.reshape(sp1b, [128,32]), axis=1)
    b2 = np.nansum(np.reshape(sp2b, [128,32]), axis=1)

    invalid = ( (b1!=0)*(b2!=0) ) == False
    b1[invalid] = np.nan
    b2[invalid] = np.nan

    ratio = b1/b2

    #fit,_ = robust_polyfit(index,ratio,3,3)
    #return np.polyval(fit,np.arange(len(sp1)))

    #plt.plot(index,ratio)
    #plt.plot(index,np.polyval(fit,index))
    #plt.show()

    ratio2 = np.zeros_like(ratio)+np.nan
    for i in range(len(ratio)):
        if np.isfinite(ratio[i]):
            i1 = i-3
            i2 = i+4
            ratio2[i] = nanmedian(ratio[i1:i2])

    keep = np.isfinite(ratio2)
    index = index[keep]
    ratio2 = ratio2[keep]

    if len(ratio2)<4:
        return np.zeros_like(sp1)+nanmedian(ratio2)
    else:
        return ius(index,ratio2,k=2,ext=3)(np.arange(len(sp1)))



@jit(nopython=True)
def odd_ratio_mean(value,err, odd_ratio = 1e-4, nmax = 10):
    #
    # Provide values and corresponding errors and compute a
    # weighted mean
    #
    #
    # odd_bad -> probability that the point is bad
    #
    # nmax -> number of iterations
    keep = np.isfinite(value)*np.isfinite(err)

    if np.sum(keep) == 0:
        return np.nan,np.nan

    value = value[keep]
    err = err[keep]

    guess = np.nanmedian(value)

    nite = 0
    while (nite < nmax):
        nsig = (value-guess)/err
        gg = np.exp(-0.5*nsig**2)
        odd_bad = odd_ratio/(gg+odd_ratio)
        odd_good = 1-odd_bad

        w = odd_good/err**2

        guess = nansum(value*w)/np.nansum(w)
        nite+=1

    bulk_error =  np.sqrt(1/np.nansum(odd_good/err**2))

    return guess,bulk_error


def lin_mini(vector, sample):
    # wrapper function that sets everything for the @jit later
    # In particular, we avoid the np.zeros that are not handled
    # by numba

    # size of input vectors and sample to be adjusted
    sz_sample = sample.shape  # 1d vector of length N
    sz_vector = vector.shape  # 2d matrix that is N x M or M x N

    # define which way the sample is flipped relative to the input vector
    if sz_vector[0] == sz_sample[0]:
        case = 2
    elif sz_vector[0] == sz_sample[1]:
        case = 1
    else:
        emsg = ('Neither vector[0]==sample[0] nor vector[0]==sample[1] '
                '(function = {0})')
        print(emsg)
        raise ValueError(emsg.format(func_name))

    # we check if there are NaNs in the vector or the sample
    # if there are NaNs, we'll fit the rest of the domain
    isnan = (np.sum(np.isnan(vector)) != 0) or (np.sum(np.isnan(sample)) != 0)

    if case == 1:

        if isnan:
            # we create a mask of non-NaN
            keep = np.isfinite(vector) * np.isfinite(np.sum(sample, axis=0))
            # redefine the input vector to avoid NaNs
            vector = vector[keep]
            sample = sample[:, keep]

            sz_sample = sample.shape
            sz_vector = vector.shape

        # matrix of covariances
        mm = np.zeros([sz_sample[0], sz_sample[0]])
        # cross-terms of vector and columns of sample
        v = np.zeros(sz_sample[0])
        # reconstructed amplitudes
        amps = np.zeros(sz_sample[0])
        # reconstruted fit
        recon = np.zeros(sz_sample[1])

    if case == 2:
        # same as for case 1, but with axis flipped
        if isnan:
            # we create a mask of non-NaN
            keep = np.isfinite(vector) * np.isfinite(np.sum(sample, axis=1))
            vector = vector[keep]
            sample = sample[keep, :]

            sz_sample = sample.shape
            sz_vector = vector.shape

        mm = np.zeros([sz_sample[1], sz_sample[1]])
        v = np.zeros(sz_sample[1])
        amps = np.zeros(sz_sample[1])
        recon = np.zeros(sz_sample[0])

    # pass all variables and pre-formatted vectors to the @jit part of the code
    amp_out, recon_out = linear_minimization(vector, sample, mm, v, sz_sample, case,
                                             recon, amps)

    # if we had NaNs in the first place, we create a reconstructed vector
    # that has the same size as the input vector, but pad with NaNs values
    # for which we cannot derive a value
    if isnan:
        recon_out2 = np.zeros_like(keep) + np.nan
        recon_out2[keep] = recon_out
        recon_out = recon_out2

    return amp_out, recon_out


# @jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @nji
def linear_minimization(vector, sample, mm, v, sz_sample, case, recon, amps):
    # raise ValueError(emsg.format(func_name))
    # ​
    # vector of N elements
    # sample: matrix N * M each M column is adjusted in amplitude to minimize
    # the chi2 according to the input vector
    # output: vector of length M gives the amplitude of each column
    #
    if case == 1:
        # fill-in the co-variance matrix
        for i in range(sz_sample[0]):
            for j in range(i, sz_sample[0]):
                mm[i, j] = np.sum(sample[i, :] * sample[j, :])
                # we know the matrix is symetric, we fill the other half
                # of the diagonal directly
                mm[j, i] = mm[i, j]
            # dot-product of vector with sample columns
            v[i] = np.sum(vector * sample[i, :])
        # if the matrix cannot we inverted because the determinant is zero,
        # then we return a NaN for all outputs
        if np.linalg.det(mm) == 0:
            amps = np.zeros(sz_sample[0]) + np.nan
            recon = np.zeros_like(v)
            return amps, recon

        # invert coveriance matrix
        inv = np.linalg.inv(mm)
        # retrieve amplitudes
        for i in range(len(v)):
            for j in range(len(v)):
                amps[i] += inv[i, j] * v[j]

        # reconstruction of the best-fit from the input sample and derived
        # amplitudes
        for i in range(sz_sample[0]):
            recon += amps[i] * sample[i, :]
        return amps, recon

    if case == 2:
        # same as for case 1 but with axis flipped
        for i in range(sz_sample[1]):
            for j in range(i, sz_sample[1]):
                mm[i, j] = np.sum(sample[:, i] * sample[:, j])
                mm[j, i] = mm[i, j]
            v[i] = np.sum(vector * sample[:, i])

        if np.linalg.det(mm) == 0:
            return amps, recon

        inv = np.linalg.inv(mm)
        for i in range(len(v)):
            for j in range(len(v)):
                amps[i] += inv[i, j] * v[j]

        for i in range(sz_sample[1]):
            recon += amps[i] * sample[:, i]
        return amps, recon

