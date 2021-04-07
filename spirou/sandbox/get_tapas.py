import os
import glob
from astropy.table import Table
import matplotlib.pyplot as plt

# Inputs for the query, update only the request ID as sent in the TAPAS email
request = 7381 # request

tapas_path = '/Volumes/courlan/TAPAS/' # big drive with space to spare
user = 'artigau@astro.umontreal.ca' # user ID


# get curret directory
current_dir = os.getcwd()

# create folder for run
if not os.path.isdir(str(request)):
    # create folder for TAPAS run
    os.system('mkdir {0}{1}'.format(tapas_path,request))
    # change to tapas directory
    os.chdir('{0}{1}'.format(tapas_path,request))

    # loop through 12 files in case this is a double request
    for id in range(1,12):
        url = 'ftp://tapas:tapas@ftp.ipsl.fr/{0}/Ether_TAPAS_{1}/tapas_0000{2}.ipac.gz'.format(user,
                                                                                               request,
                                                                                               str(id).zfill(2))
        print('wget '+url)
        # get the URL
        os.system('wget {}'.format(url))
    # unzip the ipac files
    os.system('gunzip *.gz')
    # find the ipac files in the folder
    files = glob.glob('tapas_??????.ipac')

    # loop through files and plot the spectra
    for file in files:
        outname = file.split('.')[0]+'.fits'
        print(outname)
        if not os.path.isfile(outname):
            tbl = Table.read(file,format='ascii')
            tbl.write(outname)

            plt.plot(tbl['wavelength'],tbl['transmittance'],label = file)
    plt.xlabel('wavelength')
    plt.ylabel('transmittance')
    plt.legend()
    plt.savefig('{}_transmittance.pdf'.format(request))
    plt.show()

# back to initial directory
os.chdir(current_dir)