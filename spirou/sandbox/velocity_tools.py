import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as ius
import etienne_tools as et
from scipy import constants
import os
from tqdm import tqdm
from scipy.interpolate import interp1d
from PyAstronomy.pyasl import rotBroad
import glob

def bisector(rv,ccf,doplot=False,low_high_cut=0.1,figure_title='',
        ccf_plot_file='',showplots=True):

    # get minima
    imin = int(np.argmin(ccf))
    # print(imin,type(imin))

    # get point where the derivative changes sign at the edge of the line
    # the bisector is ambiguous passed this poind
    width_blue = imin - np.max(np.where(np.gradient(ccf[0:imin]) > 0))
    width_red = np.min(np.where(np.gradient(ccf[imin:]) < 0))

    # get the width from the side of the center that reaches
    # that point first
    width = int(np.min([width_blue, width_red]))

    # set depth to zero
    ccf -= np.min(ccf)

    # set continuum to one
    ccf /= np.min(ccf[[imin - width, imin + width]])

    # interpolate each side of the ccf slope at a range of depths
    depth = np.arange(low_high_cut, 1-low_high_cut, 0.001)

    # blue and red side of line
    g1 = ((ccf[imin:imin - width:-1] > low_high_cut)
          & (ccf[imin:imin - width:-1] < (1-low_high_cut)))
    spline1 = ius(
            ccf[imin:imin - width:-1][g1],
            rv[imin:imin - width:-1][g1],
            k=2
            )

    g2 = ((ccf[imin:imin + width] > low_high_cut)
          & (ccf[imin:imin + width] < (1-low_high_cut)))
    spline2 = ius(
            ccf[imin:imin + width][g2],
            rv[imin:imin + width][g2],
            k=2)

    # get midpoint
    bisector_position = (spline2(depth)+spline1(depth))/2

    # get bisector widht
    width_ccf = (spline2(depth)-spline1(depth))

    if doplot:
        # some nice plots
        fig = plt.figure()
        plt.plot(
                rv[imin - width:imin + width],
                ccf[imin - width:imin + width],
                label='ccf')
        plt.plot(bisector_position, depth, label='bisector')
        plt.plot(
                ((bisector_position-np.mean(bisector_position))*100
                 + np.mean(bisector_position)),
                depth,
                label='bisector * 100',
                 )
        plt.legend()
        plt.title(figure_title)
        plt.xlabel('Velocity (km/s)')
        plt.ylabel('Depth')
        if ccf_plot_file != '':
            plt.savefig(ccf_plot_file)
        if showplots:
            plt.show()
        plt.close(fig)

    # define depth in the same way as Perryman, 0 is top, 1 is bottom
    return 1-depth, bisector_position, width_ccf

def systemic_velo(file,teff = 3600, logg = 5.50, smart = False, doplot = True,
                  model_path = '/Volumes/courlan/HiResFITS/'):

    """
    Pass an s1d or template file and get the FWHM
    """

    allowed_logg = [4.5,5.0,5.5,6.0]
    if logg not in allowed_logg:
        err_msg = """
        !!! ERROR !!!'
        'Logg value is {0}, but it must on of the following values :'
        {1}
        """.format(logg, ', '.join(np.array(allowed_logg,dtype='<U99')))
        raise ValueError(err_msg)

    # read s1d or template as a table
    tbl = Table.read(file)
    hdr = fits.getheader(file,ext=1)

    if smart:
        teff = int(np.round(hdr['OBJTEMP'], -2))

    # get the mask corresponding to the temperature of your object
    mask = Table.read(model_path+
                      '/lte0{0}-{1:.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.csv'.format(teff,logg))

    if 'Template' not in file:
        # read header to get BERV
        hdr = fits.getheader(file,ext=1)
        berv = hdr['BERV']*1000
    else:
        print('This is a template, no berv adjustment')
        berv = 0

    # wavelength and flux from table. Converting to numpy array to speed-up
    wave = np.array(tbl['wavelength'])
    sp = np.array((tbl['flux']))


    # save for the mask
    wave_mask = np.array(mask['wavelength'])
    weight_mask = np.array(mask['weight'])

    # speedier than keep all the mask lines
    good_mask = (wave_mask>np.min(wave))*(wave_mask<np.max(wave))
    weight_mask,wave_mask = weight_mask[good_mask], wave_mask[good_mask]

    # valid domain is where we have finite spectrum
    valid_sp = np.isfinite(sp)

    spl_valid = ius(wave,valid_sp,k=1,ext=1)

    good_line = np.ones_like(wave_mask,dtype = bool)

    # coarse scan to reject lines that are affected at some point by a gap in the spectrum
    dv = np.arange(-200,200,3,dtype = float)*1000

    for i in range(len(dv)):
        good_line*=(spl_valid(et.doppler(wave_mask,-dv[i]+berv))>0.99)

    # rejecting lines bad at least one
    wave_mask = wave_mask[good_line]
    weight_mask = weight_mask[good_line]

    spl = ius(wave[valid_sp],sp[valid_sp],k=1,ext=1)

    # finer scan
    dv = np.arange(-200,200,.5,dtype = float)*1000
    ccf = np.zeros_like(dv)

    for i in range(len(dv)):
        ccf[i] = np.sum(weight_mask*spl(et.doppler(wave_mask,-dv[i]+berv)))

    fit,_ = et.robust_polyfit(dv, ccf, 2, 3)

    ccf/=np.polyval(fit,dv)

    if doplot:
        plt.plot(dv/1000,ccf)
        plt.xlabel( 'velocity [km/s]')
        plt.ylabel('contrast')
        # plt.plot(dv/1000,et.super_gauss(dv,*fit))
        plt.title(hdr['OBJECT'])
        plt.show()
    #cen, ew, amp,expo, zp, slope
    p0 =dv[np.argmin(ccf)],5000,np.min(ccf)-np.nanmedian(ccf),1.4,np.nanmedian(ccf),0
    fit = et.fit_super_gauss(dv,ccf,p0)

    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('Object : {}'.format(hdr['OBJECT']))
    if 'OBJTEMP' in hdr:
        print('Header effective temperature {} K'.format(hdr['OBJTEMP']))
        print('Input effective temperature {} K'.format(teff))
    print()
    print('~~~~ Super-gaussian fit properties ~~~~')
    print(' velocity : {:.3f} km/s'.format(fit[0]/1000))

    fwhm_fit = 2*fit[1]*(-2*np.log(0.5))**(1/fit[3])
    print(' FWHM : {:.3f} km/s'.format(fwhm_fit/1000))

    # get bisector properties
    depth, bisector_position, width_ccf = bisector(dv/1000, ccf,doplot=doplot,low_high_cut=0.2 )

    v0 = ius(depth[np.argsort(depth)],bisector_position[np.argsort(depth)])(0.5)
    fwhm_bis = ius(depth[np.argsort(depth)],width_ccf[np.argsort(depth)])(0.5)

    print('~~~~ bisector properties ~~~~')
    print(' Radial velocity : {0:.3f} km/s'.format(v0))
    print(' FWHM : {0:.3f} km/s'.format(fwhm_bis))

    out_dict = dict()
    out_dict['FWHM_FIT'] = fwhm_fit
    out_dict['VELOCITY_FIT'] = fit[0]
    out_dict['FWHM_BISECTOR'] = fwhm_bis
    out_dict['VELOCITY_BISECTOR'] = v0
    out_dict['OBJECT'] = hdr['OBJECT']
    out_dict['LOGG'] = logg
    out_dict['TEFF'] = teff


    return out_dict

def mk_model_mask(model_path = '/Volumes/courlan/HiResFITS/'):
    # some parameters, don't worry
    dv = 0.00  # km/s -- width of the CCF box
    c = constants.c/1000.0 # speed of light

    # create directory if needed
    if not os.path.isdir(model_path):
        os.system('mkdir {0}'.format(model_path))

    # read wavelength and flux. The wavelength is expressed in Ang, we convert to µm
    ftp_link = 'ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/'
    wave_file = 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
    if not os.path.isfile(model_path+'/'+wave_file):
        os.system('wget {0}{1}'.format(ftp_link,wave_file) )
        os.system('mv {0} {1}'.format(wave_file,model_path))
    w = fits.getdata(model_path+'/'+wave_file) / 10

    good_domain = (w>350)*(w<2500)
    w = w[good_domain]

    for logg in ([4.0,5.0,5.5,6.0]):
        # get goettigen models if you don't have them.
        for temperature in np.arange(3500, 6100, 100):
            outname = '{0}/lte0{1}-{2:.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(model_path,temperature,logg)

            print('We are downloading/checking presence of model file {0}'.format(outname))

            # we need to remove the '.' in the outname_mask as it messes the name of the ccf file
            outname_mask = outname.split('.fits')[0]+'.csv'
            if os.path.isfile(outname_mask):
                print('File {0} exists, we are happy!'.format(outname_mask))
                continue

            if not os.path.isfile(outname):
                os.system(
                    'wget {0}PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte0{1}-{2:.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(ftp_link,temperature,logg))

                os.system('mv lte0{1:.0f}-{2:.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits {0}'.format(model_path,temperature,logg))
            else:
                print('File {0} exists, we are happy!'.format(outname))

            print('\tReading data and finding derivatives')
            f = fits.getdata(outname)
            f = f[good_domain]
            f*=w # express in photons, not in energy as models provide things normally
            f/=np.nanmedian(f)
            # find the first and second derivative of the flux. Derivative as a function of log(lambda) is
            # derivative in the velocity space
            df = np.gradient(f)/np.gradient(np.log(w ))
            ddf = np.gradient(df)

            # lines are regions there is a sign change in the derivative of the flux
            # we also have some checks for NaNs
            line = np.where((np.sign(df[1:]) != np.sign(df[:-1])) &
                            np.isfinite(ddf[1:])
                            & np.isfinite(df[1:])
                            & np.isfinite(df[:-1]) & (ddf[1:]>0))[0]

            weight = ddf[line]

            for ite in range(9):
                print(len(line))
                keep = np.ones_like(weight,dtype = bool)

                for i in range(5,len(line)-6):
                    # line is less than 1/5th of best neighbour from -5 to +5 in the list
                    if weight[i]/np.max(weight[i-5:i+5]) < 0.1:
                        keep[i] = False
                line = line[keep]
                weight = weight[keep]

            print('\tConstructing table')
            # create the output table
            tbl = Table()
            tbl['wavelength'] = np.zeros_like(line, dtype=float)
            tbl['weight'] = weight
            tbl['weight'] /=np.nanmedian(np.abs(tbl['weight']))
            tbl['weight'] = np.round(tbl['weight'],4)

            print('\tcomputing line positions')
            wave_cen = np.zeros_like(tbl,dtype = float)
            for i in tqdm(range(len(line))):
                # we perform a linear interpolation to find the exact wavelength
                # where the derivatives goes to zero
                wave_cen[i] = (np.polyfit(df[line[i]:line[i] + 2], w[line[i]:line[i] + 2], 1))[1]
            tbl['wavelength'] = np.round(wave_cen,6)

            tbl = tbl[np.isfinite(tbl['wavelength'])]

            print('\twriting {0}'.format(outname_mask))
            tbl.write(outname_mask, overwrite=True)

def epsilon(wave):
    # for a wavelength in nm, return the Epsilon value
    # of limb darkening. Here we have a place-holder law
    # It is 1.0 for 900 nm and 0.1 for 2500 nm
    #
    # it has to be between 0 and 1 for all wavelength

    eps = np.polyval(np.polyfit([900, 2500], [1.0, 0.1], 1), wave)

    return eps

def get_broad_lsf(wave0, epsilon, vsini):
    # create a rotation profile
    wave_tmp = (1 + np.arange(-np.int(vsini) - 1, np.int(vsini) + 1) / (constants.c / 1000)) * wave0 * 10
    f_tmp = np.zeros_like(wave_tmp)
    f_tmp[np.int(vsini)] = 1

    lsf = rotBroad(wave_tmp, f_tmp, epsilon, vsini)
    lsf = lsf[lsf != 0]
    # print(lsf)
    # print(wave0)
    return lsf

def get_lsf(beta=2.24, ew=2.20):
    """
    wave and flux of a spectrum

    we assume that the input spectrum is a SPIRou s1d_v
    the

    LSF shape : np.exp( -0.5 * np.abs( ( wave/wave0 - 1)*c/ew  )**beta )

    """

    dvs = np.arange(-np.ceil(ew * 5), np.ceil(ew * 5) + 1)

    lsf = np.exp(-0.5 * np.abs(dvs / ew) ** beta)

    lsf /= np.sum(lsf)

    return dvs, lsf

def model2spirougrid(spirou_file,teff = 3500, logg = 4.5, doplot = True,model_path = '/Volumes/courlan/HiResFITS/'):
    """
     please download a model and the corresponding wavelength grid from here:

    http://phoenix.astro.physik.uni-goettingen.de/?page_id=15

        these are expressed in energy while SPIRou spectra are expressed in
        photons, do not forget to multiply the model by lambda to get a spectrum
        proportional to photons

    """

    model_file = model_path+'/lte0{0}-{1:.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(
        teff,logg)

    outname = model_file.split('.fits')[0] + '_fake_sp.fits'
    if os.path.isfile(outname):
        print('{} exists'.format(outname))
        return outname

    # read a spirou template to get the wavelength grid. This has to be a _v_ file
    wave_spirou = np.array(Table.read(spirou_file)['wavelength'])
    tbl0 = Table.read(spirou_file)
    hdr0 = fits.getheader(spirou_file, ext=0)
    hdr1 = fits.getheader(spirou_file, ext=1)

    hdr1['OBJET'] = 'MODEL'
    hdr1['OBJTEMP'] = teff

    # we read the model wavelength grid
    # Careful with the various definition of wavelength, SPIRou uses
    # nm and many models use Angs. Here we added a /10 factor to convert to
    # nm.
    wave_model = fits.getdata(model_path+'/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits') / 10.0

    flux_model = fits.getdata(model_file)

    # velocity parameters
    dv0 = 0.0  # velocity difference in km/s
    vsini = 0.01  # in km/s

    # LSF parameters
    beta = 2.24  # shape factor exponent
    ew = 2.20  # e width factor

    # number of wavelength samplings you want for your limb darkening law
    nlimb = 10  # 10 seems like a fair number

    s1d = interp1d(wave_model, flux_model)

    flux_interpol = np.zeros_like(wave_spirou)

    dvs, lsf = get_lsf(beta=beta, ew=ew)

    wave_limb = np.arange(np.min(wave_spirou), np.max(wave_spirou), np.max(wave_spirou) / nlimb)
    wave_limb = np.append(wave_limb, np.max(wave_spirou))

    lsfs = []
    for i in range(len(wave_limb)):
        # get the stellar broadening LSF
        lsf2 = get_broad_lsf(wave_limb[i], epsilon(wave_limb[i]), vsini)

        # get the combined instrumental and broadening LSF
        lsf_full = np.convolve(lsf, lsf2, mode='full')
        lsf_full /= np.nansum(lsf_full)

        lsfs.append(lsf_full)

    lsfs = np.array(lsfs)

    # get the corresponding velocity grid
    dvs = np.arange(len(lsf_full), dtype=float)
    dvs -= np.mean(dvs)

    # perform the convolution
    for i in range(len(dvs)):
        print('Dv bin {0} in {1} of lsf, normalisation = {2}'.format(i + 1, len(dvs), lsf_full[i]))

        lsf_int = interp1d(wave_limb, lsfs[:, i])

        flux_interpol += s1d(wave_spirou * (1 + (dvs[i] - dv0) / (constants.c / 1000))) * lsf_int(wave_spirou)

    tbl0['flux'] = flux_interpol

    hdu0 = fits.PrimaryHDU()
    hdu1 = fits.BinTableHDU(tbl0)
    hdu0.header = hdr0
    hdu1.header = hdr1

    # convert back from dictionnary to table and save
    new_hdul = fits.HDUList([hdu0, hdu1])
    new_hdul.writeto(outname, overwrite=True)

    if doplot:
        xrange = [1800, 1805]
        # plot the resulting convolved spectrum
        plt.plot(wave_model, flux_model, 'k-', label='input spectrum')
        plt.xlim(xrange)
        g = (wave_model > xrange[0]) * (wave_model < xrange[1])
        plt.ylim([np.min(flux_model[g]), np.max(flux_model[g])])

        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Flux (arbitrary)')
        plt.plot(wave_spirou, flux_interpol, 'r-', label='convolved spectrum')
        plt.title(
            'star : vsini={0}km/s, epsilon={1}, dv={2}km/s\nLSF: beta={3}, ew={4}'.format(vsini, epsilon, dv0, beta, ew))
        plt.legend()
        plt.savefig('comp_convolve_sp.pdf')
        plt.show()

    return outname

def get_autocorrelation(file,doplot = True):

    ac_name = file.split('.fits')[0]+'_autocorrelation.csv'

    if not os.path.isfile(ac_name):
        hdr = fits.getheader(file,ext=1)
        tbl = Table.read(file)

        wave = np.array(tbl['wavelength'])
        flux = np.array(tbl['flux'])
        flux -= et.lowpassfilter(flux,1000)

        g = (wave>2250)*(wave<2350)
        wave = wave[g]
        flux = flux[g]


        ac_velo = np.arange(-200,201,1.0,dtype = float)
        ac = np.zeros_like(ac_velo)

        valid = np.isfinite(flux)
        spl_mask = ius(wave,np.array(valid,dtype = float),k=1,ext=1)
        spl = ius(wave[valid],flux[valid],k=2,ext=1)

        norm = np.mean(flux[valid]**2)

        print('Computing auto-correlation function')
        for i in tqdm(range(len(ac_velo))):
            if ac[i] ==0:
                sp1 = spl(et.doppler(wave,1000*ac_velo[i]))
                mask = spl_mask(et.doppler(wave,1000*ac_velo[i]))

                keep = (mask>0.99)*valid

                tmp = np.mean(sp1[keep]*flux[keep])/norm
                ac[np.abs(ac_velo) == np.abs(ac_velo[i])] = tmp

        tbl = Table()
        tbl['VELOCITY'] = ac_velo
        tbl['AUTOCORRELATION'] = ac

        tbl.write(ac_name)

    print('Reading {}'.format(ac_name))

    tbl = Table.read(ac_name)
    ac_velo = np.array(tbl['VELOCITY'])
    ac =np.array( tbl['AUTOCORRELATION'])


    #lsf = np.array(ac)

    #for ite in range(19):
    #    lsf = lsf/np.sum(lsf**2)
    #    lsf2 = np.convolve(lsf,lsf,mode = 'same')
    #
    #    corr = ac/np.sum(ac**2) - lsf2/np.sum(lsf2**2)
    #    lsf+=corr



    #if doplot:
    """
    fig, ax = plt.subplots(nrows = 1, ncols = 1,sharex = True)
    ax[0].plot(ac_velo,lsf2/np.sum(lsf2),'g.',alpha =0.5,label = 'reconstructed AC')
    ax[0].plot(ac_velo,ac/np.sum(ac),color = 'red',alpha=0.5,label = 'observed AC')
    ax[0].plot(ac_velo,ac/np.sum(ac) - lsf2/np.sum(lsf2),'.',alpha=0.5,label = 'residual AC')
    ax[0].set(title=file+'\nAC',xlabel = 'velocity [km/s]',ylabel = 'ac power')
    ax[0].legend()


    ax[1].plot(ac_velo,lsf,alpha = 0.5,color = 'orange', label = 'effective line profile')
    ax[1].set(title=file+'\neffective profile',xlabel = 'velocity [km/s]',ylabel = 'intensity')
    ax[1].legend()

    plt.tight_layout()
    plt.savefig(file.split('.fits')[0]+'_LINE.pdf')
    plt.show()
    """

    plt.clf()
    plt.plot(ac_velo/np.sqrt(2),ac)
    plt.xlabel('velocity [km/s]')
    plt.ylabel('autocorrelation')
    plt.title(file)
    plt.savefig(file.split('.fits')[0]+'_LINE.pdf')
    if doplot:
        plt.show()
    plt.clf()


    try:
        depth, bisector_position, width_ccf = bisector(ac_velo/np.sqrt(2),1-ac)
        fwhm_fit = ius(depth[np.argsort(depth)],width_ccf[np.argsort(depth)],)(0.5)
        print('FWHM = {0:.4f} km/s'.format(fwhm_fit))
    except:
        return -1

    return fwhm_fit




if True:
    files = glob.glob('templates/Template*s1d*fits')

    for file in files:
        try:
            tbl = Table(fits.getdata(file,ext=2))
            if len(tbl)<10:
                cmd = 'rm '+file
                print(cmd)
                os.system(cmd)
        except:
            print('Err')
    files1 = np.array(glob.glob('templates/Template*s1d*fits'))
    files2 = np.array(glob.glob('/Volumes/courlan/HiResFITS/*fake_sp.fits'))
    files = np.append(files1,files2)

    tbl = Table()
    tbl['OBJECT'] = np.zeros_like(files)
    tbl['file'] = files
    tbl['teff'] = np.zeros_like(files,dtype = float)
    tbl['FWHM'] = np.zeros_like(files,dtype = float)
    for i in range(len(files)):
        print(i,len(files))

        if 'FP' in files[i]:
            continue

        fwhm = get_autocorrelation(files[i],doplot = False)
        hdr = fits.getheader(files[i],ext=1)
        tbl['OBJECT'][i] = hdr['OBJECT']
        tbl['teff'][i] = hdr['OBJTEMP']
        tbl['FWHM'][i] = fwhm


    tbl['FWHM'][tbl['FWHM']<1] = np.nan
    tbl['teff'][tbl['teff']>4000] = np.nan

    tbl['type'] = 'observation'

    for i in range(len(tbl)):
        if '-4.50-' in tbl['file'][i]:
            tbl['type'][i] = 'logg = 4.5'
        if '-5.00-' in tbl['file'][i]:
            tbl['type'][i] = 'logg = 5.0'
        if '-5.50-' in tbl['file'][i]:
            tbl['type'][i] = 'logg = 5.5'

        if '699' in tbl['file'][i]:
            tbl['type'][i] = 'GL699'
        if 'TRAP' in tbl['file'][i]:
            tbl['type'][i] = 'Trappist-1'
        if 'AUM' in tbl['file'][i]:
            tbl['type'][i] = 'AUMic'


    for tp in  np.unique(tbl['type']):
        g = tbl['type'] == tp

        plt.plot(tbl['teff'][g],tbl['FWHM'][g],'o',label = tp)

    plt.ylabel('AC fwhm [km/s]')
    plt.xlabel('teff [K]')
    plt.legend()
    plt.savefig('teff_vsini.pdf')
    plt.show()