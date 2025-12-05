from astropy.io import fits
import glob
import numpy as np
import os
import sys
from scipy.signal import convolve2d
from astropy.table import Table

def mk_isolated_nans(image):
    # input image is known to have its pixels in the right position
    # without shifts

    nans = np.isfinite(image) == 0
    ker = np.ones([3,3],dtype = float)
    neighbours = convolve2d(nans,ker,mode='same')

    isolated = nans*(neighbours==1)

    y,x = np.where(isolated)

    tbl = Table()
    tbl['x'] = x
    tbl['y'] = y
    tbl.write('isolated_nans.csv',overwrite = True)

def check_hdr(hdr):
    """

    :param hdr:

    Header from a NIRPS image to be cheked against a number of QCs

    We'll need to add more checks with time.

    :return:

    """

    keys_to_check = ['OBJECT']

    output = ''
    for key in keys_to_check:
        if key not in hdr:

            output+='missing:'+key+' '
            print('Key {0} is missing'.format(key))
        else:

            print('Key {0} is present'.format(key))
    # we consider the header all fine only if output empty string
    return output

def check_image(image):
    """

    :param hdr:

    image to be cheked against a number of QCs

    We'll need to add more checks with time.

    :return:


    """

    output = ''

    # if more than 1% NaNs, something is wrong
    if np.mean(np.isfinite(image))<0.99:
        print('Fraction of NaNs : {0}'.format(1-np.mean(np.isfinite(image))))
        output += 'too many NaNs'


    # we consider the header all fine only if output empty string
    return output

def fix_shift2(image,threshold_valid = 0.9):
    """
    Read an image and fix a row of NaNS if required

    Status:
     status = 0 -> all good, no shit
     status = 1 -> normal shift in expected direction when we have a
                   controler issue, we apply correction
     status = -1 -> opposite direction of known shift pattern, we
                    apply correction
     status = -2 -> unknown shift, we need to look at the image

    threshold_valid
        is the fraction of NaNs that are in the expected position
        should be close to but not equal to 1. 0.9 seems OK, but this may require a
        confirmation with a large number of files

    :param image:

    :return:
    """

    tbl = Table.read('isolated_nans.csv')
    ref_nans = tbl['y'],tbl['x']

    #ref_nans = np.array(fits.getdata('isolated_nans.fits'), dtype = bool)

    # getting image size
    number_amps = 32
    width_image = image.shape[0]
    # get width of amplifier ribbon
    width_amp = int(width_image / number_amps)
    # look for a discontinuity between consecutive columns, this traces
    # reference vs science pixels

    frac_nans = np.nanmean(~np.isfinite(image[ref_nans]))

    print('no shift, frac of NaNs in right place: ',frac_nans)
    if frac_nans>threshold_valid:
        return image,0


    # try +1 and -1 roll
    image2 = np.array(image)

    # the usual correction
    for i in range(number_amps):
        offset = 1 - (i % 2) * 2

        start = i * width_amp
        end = i * width_amp + width_amp

        image2[:, start:end] = np.roll(image2[:, start:end], offset, axis=1)

    frac_nans = np.nanmean(~np.isfinite(image2[ref_nans]))

    print('usual shift, frac of NaNs in right place: ',frac_nans)
    if frac_nans>threshold_valid:

        # we correct the NaN colums that appear at 256*N by interpolating
        # neighbouring columns
        for i in range(16):
            prev = image2[:,i*256-2]
            next = image2[:,i*256+1]

            image2[:, i * 256 - 1] = 2/3*prev+1/3*next
            image2[:, i * 256 ] = 1/3*prev+2/3*next

        return image2,1

    image2 = np.array(image)
    for i in range(number_amps):
        offset = (i % 2) * 2 - 1

        start = i * width_amp
        end = i * width_amp + width_amp

        image[:, start:end] = np.roll(image[:, start:end], offset, axis=1)

    frac_nans = np.nanmean(~np.isfinite(image2[ref_nans]))
    print('opposite of usual shift, frac of NaNs in right place: ',frac_nans)
    if frac_nans>threshold_valid:
        return image2,-1

    return image, -2

"""
Code to parse files--

A night directory
Search for files matching a given pattern
Construct an 'out name' and check if it exists, skip if it does
if it exists
Open the fits file

sanity of header :
    Check if OBJECT in header
    checks for other things that could go wrong
if wrong: 
   skip, don't write to new folder unless you can guess
   -> add to table of bad files and why

Check for fraction of saturated pixels and skip if f>0.05
   -> add to table of bad files and why


Check for shifted pix and correct

Save to outname
"""

if __name__ == '__main__':


    tbl = Table()

    # get nightname from arguments
    args = sys.argv


    # two arguments are given, the first one is the path to
    # science files and second one is destination folder for
    # corrected and valid frames.

    # check that outdir exists, if not, create it
    if os.path.isdir(args[2]) == False:
        print('\n\tWe create {0}\n'.format(args[2]))
        os.mkdir(args[2])

    if len(args) == 3:
        # get all files
        files = glob.glob(os.path.join('.', args[1], '*.fits'))

        tbl['FILE'] = files
        tbl['STATUS'] = np.zeros_like(files,dtype = 'U999')
        tbl['OUTNAME'] = np.zeros_like(files, dtype = 'U999')

        # loop around files and see if we need to fix them
        for it, filename in enumerate(files):

            # get the output name from the filename.
            outname = args[2]+'/'+filename.split('/')[-1]

            if os.path.isfile(outname):
                tbl['OUTNAME'][it] = outname # it exists, therefore we have an 'outname' value
                print('File {0} already exists, we skip'.format(outname))
                tbl['STATUS'][it]+= 'exists '
                continue

            # print progress
            pargs = [os.path.basename(filename), it + 1, len(files)]
            print('Processing file {0} ({1}/{2})'.format(*pargs))
            # load data HDU
            hdu = fits.open(filename)

            try: # we attempt to read, it may be corrupted.

                # load the data into image 1
                image1 = np.array(hdu[1].data)
                # get header
                hdr = hdu[0].header
            except:
                tbl['STATUS'][it] += 'corrupted file'
                continue

            tbl['STATUS'][it] += check_hdr(hdr)
            tbl['STATUS'][it] += check_image(image1)

            if tbl['STATUS'][it] != '':
                print('Status = {0}, we skip\n'.format(tbl['STATUS'][it]))
                # we need to skip, something wrong with HDR or image
                continue



            # see / fix data
            image2, status = fix_shift2(image1)

            if status == 0:
                tbl['STATUS'][it]+=', no shift'

            if status == (1):
                tbl['STATUS'][it]+=', known pix shift'

            if status == (-1):
                tbl['STATUS'][it]+=', opposite to known pix shift'

            if status == (-2):
                tbl['STATUS'][it]+=', bad pix shift'
                continue

            # push data back into hdu
            hdu[1].data = image2
            # write hdu to file

            tbl['OUTNAME'][it] = outname # we keep track of output name
            print('\n\twe write : {0}\n'.format(outname))
            hdu.writeto(outname, overwrite=True)
            # close hdu
            hdu.close()

        tbl.write(args[2]+'/status.csv', overwrite=True)
