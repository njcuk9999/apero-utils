import glob
from astropy.io import fits
import os

def clean_directory(obj, delete = False):
    # checks that a folder with a target name only has files
    # with that object as a name.
    #   If delete == False
    #       only lists files that have the wrong object
    #   If delete == True
    #       deletes files with the wrong object

    files = glob.glob(obj+'/*.fits')

    for file in files:
        hdr = fits.getheader(file)

        if len(hdr) < 20:
            hdr = fits.getheader(file, ext=1)

        if hdr['OBJECT'] != obj:
            if delete == False:
                print(file+'\t'+hdr['OBJECT']+' != '+obj)
            else:
                cmd = 'rm '+file
                print(cmd)
                os.system(cmd)