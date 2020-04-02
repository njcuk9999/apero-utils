from astropy.io import fits
import numpy as np
from photutils import DAOStarFinder, morphology,data_properties
import matplotlib.pyplot as plt
import os
from astropy.table import Table
from tqdm import tqdm
from nirps_pp import rot8

def bulk_orientation(file,force = False,nsigma_cut = 25, rotation = 0, parano = False):
    # pass a DARK_FP or FP_DARK file and determine if the
    # output file has the nominal SPIRou and HARPS-like orientation
    # with redder orders at the bottom (i.e., small y) and longer
    # wavelengths along ordres to the left (i.e. small x) of the
    # array.
    #
    # usage : bulk_orientation("sample_DARK_FP.fits")
    # returns : True/False that the orientation is nominal
    #
    # You can test a rot8 parameter for that file to see
    # which rotation should be used to match the nominal orientation
    # This "rotation" value is an integer running from 0 to 7 that defines
    # the rotation and/or symmetry to be applied to the image
    #
    # You may change the nsigma_cut = 25 default value, but with
    # a proper DARK_FP, this should not be necessary.
    #
    # Set to parano = True if you want to see an additional plot to double-check
    # that the fit of distance vs x and y is indeed improving the scatter.
    #

    tbl_name = file.split('.fits')[0]+'_rot'+str(rotation)+'.csv'
    outname_png =  file.split('.fits')[0]+'_rot'+str(rotation)+'.png'

    im = fits.getdata(file)
    im = rot8(im,rotation)

    # rotate the image according to the rot8 parameter
    print()
    print('Table name : ',tbl_name)
    print('File exists : ',os.path.isfile(tbl_name))
    print('Force : ',force)
    print()

    # we check if file exists or force = True
    if (os.path.isfile(tbl_name) == False) or (force == True):
        nsigma_cut = 25

        print('finding 1-sigma std')
        std = np.nanpercentile(np.abs(im), 68)

        print('finding PSFs in the image, {0} std from noise'.format(nsigma_cut))

        daofind = DAOStarFinder(fwhm=3.0, threshold=nsigma_cut * std)

        sources = Table(daofind(im))
        sources.write(tbl_name, overwrite = True)
    else:
        print('File {0} already exists, we read it'.format(tbl_name))

    # load source table
    sources = Table.read(tbl_name)

    # get numpy arrays instead of tables for x and y
    x = np.array(sources['xcentroid'])
    y = np.array(sources['ycentroid'])

    # find distance to nearest neighbours in x and y
    dx = np.zeros_like(x)
    dy = np.zeros_like(x)
    for i in tqdm(range(len(x))):
        dd = np.sqrt( (x[i]-x)**2+(y-y[i])**2 ) + 1e9*(x[i] == x)
        id_near = np.argmin(dd)
        dx[i] = x[i] - x[id_near]
        dy[i] = y[i] - y[id_near]

    # find abs distance to nearest neighbour
    dnear = np.sqrt( dx**2 + dy**2)

    # only keep points that are not isolate
    good = dnear < 20 # FP peaks must be within <20 pixels of each other
    x,y, dx, dy, dnear = x[good], y[good], dx[good], dy[good], dnear[good]

    # median median abs distance to neighbours and decide
    # if we have an x or y cross-dispersion
    med_dx = np.nanmedian(np.abs(dx))
    med_dy = np.nanmedian(np.abs(dy))

    nominal = True

    # we expected the vector between adjacent FP peaks to be along the 'fast' dispersion
    crossdispersed = ''
    if ((med_dx)>8) and (med_dy<3):
        crossdispersed = 'x'

    if ((med_dx)<3) and (med_dy>8):
        crossdispersed = 'y'
        nominal = False

    # we have a problem if the above conditions have not been satisfied. We
    # plot the x/y vectors for guidance
    if crossdispersed =='':
        plt.plot(x,y,'g.', xlabel = 'x pix', ylabel = 'ypix')
        plt.show(block = True)
        print('Houston, we have a problem!')
        return []

    # find the change in FP to FP peak distance along the X axis for 1/2 of the FOV
    slope_x = np.sum((dnear-np.mean(dnear))*(x-np.mean(x)))/np.sum((x-np.mean(x))**2)*2048
    # find the change in FP to FP peak distance along the Y axis for 1/2 of the FOV
    slope_y = np.sum((dnear-np.mean(dnear))*(y-np.mean(y)))/np.sum((y-np.mean(y))**2)*2048

    if np.abs(slope_x) < 1:
        print('We expect a slope of FP gap of >1 pixel over 1/2 of the FOV in X')
        print('we have a slope of {0} over the FOV and this is too low to conclude.'.format(slope_x))
        #return #[]

    if np.abs(slope_y) < 1:
        print('We expect a slope of FP gap of >1 pixel over 1/2 of the FOV in Y')
        print('we have a slope of {0} over the FOV and this is too low to conclude.'.format(slope_y))
        #return #[]

    # we find residuals in dnear after applying slopes
    dnear_corr = dnear-(x-np.mean(x))*slope_x/2048-(y-np.mean(y))*slope_y/2048

    good = dnear_corr < 20 # FP peaks must be within <20 pixels of each other
    x,y, dx, dy, dnear, dnear_corr = x[good], y[good], dx[good], dy[good], dnear[good], dnear_corr[good]

    if parano:
        print('Check that the dispersions is greatly reduced between the before/after')
        print(' graphs')
        fig, ax = plt.subplots(nrows = 2, ncols =2)
        ax[0,0].hist2d(x,dnear,bins=20)#,'g.')
        ax[0,0].set(xlabel = 'x pixel', ylabel = 'd near (pix)', title = 'Before model',ylim = [0,20])
        ax[0,1].hist2d(y,dnear,bins=20)#,'g.')
        ax[0,1].set(xlabel = 'y pixel', ylabel = 'd near (pix)', title = 'Before model',ylim = [0,20])

        ax[1,0].hist2d(np.append(x,0),np.append(dnear_corr,19.99),bins=20)#,'g.')
        ax[1,0].set(xlabel = 'x pixel', ylabel = 'd near (pix)', title = 'After model',ylim = [0,20])
        ax[1,1].hist2d(np.append(y,0),np.append(dnear_corr,19.99),bins=20)#,'g.')
        ax[1,1].set(xlabel = 'y pixel', ylabel = 'd near (pix)', title = 'After model',ylim = [0,20])
        plt.tight_layout()
        plt.show(block = True)

    plt.imshow(im, vmin = np.nanpercentile(im,5), vmax = np.nanpercentile(im,99),origin = 'lower',
               cmap = 'inferno')
    plt.xlabel('x pix')
    plt.ylabel('y pix')
    plt.title(file+'\nrotation = '+str(rotation))

    print()
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print()

    if (crossdispersed == 'x'):
        print('cross dispersion is along the X axis')
        print('fast dispersion is along the Y axis')

        if (slope_x > 0):
            print('Blue side of orders are at small X positions')
            plt.text(1024,2048,'Blue wavelengths\nalong orders', color = 'cyan', horizontalalignment='center')

            print('Red side of orders are at large X positions')
            plt.text(3072,2048,'Red wavelengths\nalong orders', color = 'orange', horizontalalignment='center')
            nominal = False
        else:
            print('Red side of orders are at small X positions')
            plt.text(1024,2048,'Red wavelengths\nalong orders', color = 'orange', horizontalalignment='center')

            print('Blue side of orders are at large X positions')
            plt.text(3072,2048,'Blue wavelengths\nalong orders', color = 'cyan', horizontalalignment='center')



        if  (slope_y < 0):
            print('Blue orders are at large Y positions')
            plt.text(2048,3072,'Blue\norders', color = 'cyan', horizontalalignment='center')

            print('Red orders are at small Y positions')
            plt.text(2048,1024,'Red\norders', color = 'orange', horizontalalignment='center')
        else:
            print('Red orders are at large Y positions')
            plt.text(2048,3072,'Red\norders', color = 'orange', horizontalalignment='center')

            print('Red orders are at small Y positions')
            plt.text(2048,1024,'Blue\norders', color = 'cyan', horizontalalignment='center')

            nominal = False


    if (crossdispersed == 'y'):
        print('cross dispersion is along the Y axis')
        print('fast dispersion is along the X axis')

        nominal = False

        if  (slope_x < 0):
            print('Blue orders are at large X positions')
            plt.text(3072,2048,'Blue\norders', color = 'cyan', horizontalalignment='center')

            print('Red orders are at small X positions')
            plt.text(1024,2048,'Red\norders', color = 'orange', horizontalalignment='center')
        else:
            print('Red orders are at large X positions')
            plt.text(3072,2048,'Red\norders', color = 'orange', horizontalalignment='center')

            print('Blue orders are at small X positions')
            plt.text(1024,2048,'Blue\norders', color = 'cyan', horizontalalignment='center')

        if (slope_y < 0):
            print('Blue side of orders are at large Y positions')
            plt.text(2048,3072,'Blue wavelengths\nalong orders', color = 'cyan', horizontalalignment='center')

            print('Red side of orders are at small Y positions')
            plt.text(2048,1024,'Red wavelengths\nalong orders', color = 'orange', horizontalalignment='center')
        else:
            print('Red side of orders are at large Y positions')
            plt.text(2048,3072,'Red wavelengths\nalong orders', color = 'orange', horizontalalignment='center')

            print('Blue side of orders are at small Y positions')
            plt.text(2048,1024,'Blue wavelengths\nalong orders', color = 'cyan', horizontalalignment='center')


    plt.text(2048,3072+512,'Nominal\nblue orders', color = 'yellow', horizontalalignment='center')
    plt.text(2048,1024-512,'Nominal\nred orders', color = 'yellow', horizontalalignment='center')

    plt.text(3072,2048-512,'Nominal\nBlue wavelengths\nalong orders', color = 'yellow', horizontalalignment='center')
    plt.text(1024,2048-512,'Nominal\nRed wavelengths\nalong orders', color = 'yellow', horizontalalignment='center')
    print()
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    if nominal == False:
        print('     current rotation passed to rot8 is {0}'.format(rotation))
        print('     orientation is *not* nominal')
    else:
        print('     current rotation passed to rot8 is {0}'.format(rotation))
        print('     orientation is nominal')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

    print('crossdispersed = ',crossdispersed)
    print('slope_x = ',slope_x)
    print('slope_y = ',slope_y)

    plt.savefig(outname_png)
    plt.show(block = True)


    return nominal