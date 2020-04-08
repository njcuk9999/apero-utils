import numpy as np
from astropy.io import fits
import os
import sys
import matplotlib.pyplot as plt

info = """
Simple program to create a correlated double sampling image.
One passes two images as argument and a difference is made.
The two frames should be readouts from an infrared array. The
readouts do not need to be consecutive and can be any two readouts
within a sequence. 

Syntax:

python cds.py readout1.fits readout2.fits outname.fits -doplot -verbose

\tThe first two arguments (1st and 2nd readout) are
\tcompulsary. The 3rd argument (outname) is not, and if
\tit is not provided, then we will create a name for the
\toutput file from the names of the two input files.

\t-> The output name, if no value is provided, is :
\t\toutname = 'CDS_'+readout2+'_'+readout1+'.fits

\tOptional arguments: 

\t\t-doplot 
\t\tgenerates some plots on the statistical properties of the
\t\tCDS. Saves the plot file as outname.pdf

\t\t-verbose
\t\tsome stats on the CDS
"""


arg=np.asarray(sys.argv)

if len(arg) == 1:
    print(info)
    sys.exit()


arg=arg[1:] # first argument is simply the name of the program and needs to be removed

if len(arg) < 2:
    print('We need at least two arguments !')
    sys.exit()

arg = np.array(arg)

file1 = arg[0]
file2 = arg[1]

# if set, we will give more printouts
verbose = '-verbose' in arg
# if set, we will generate graphs on stats
doplot = '-doplot' in arg

if verbose:
    arg = arg[arg != '-verbose']
if doplot:
    arg = arg[arg != '-doplot']

if len(arg) == 2:
    outname = 'CDS_'+file2.split('.fits')[0]+'_'+file1.split('.fits')[0]+'.fits'
else:
    outname = arg[2]


if os.path.isfile(file1) == False:
    print('File {0} does not exist'.format(file1))
    sys.exit()

if os.path.isfile(file2) == False:
    print('File {0} does not exist'.format(file2))
    sys.exit()

im2 = np.array(fits.getdata(file2), dtype = float)
im1 = np.array(fits.getdata(file1), dtype = float)
diff = (im2 - im1)

if verbose:
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('Verbose = True')
    percentiles = [0.1,1,5,10,50,90,95,99,99.9]
    print('\tPercentiles\tCDS value')
    for percentile in percentiles:
        print('\t{0:3.1f}th\t\t{1:6.1f}'.format(percentile, np.nanpercentile(diff, percentile)))

if doplot:
    v1 = []
    v2 = []
    v3 = []

    fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = [12,6])
    percentiles = [0.1,1,2,3,5,10,20,30,50,70,80,90,95,97,98,99,99.9]
    for percentile in percentiles:
        v1.append(np.nanpercentile(im1, percentile))
        v2.append(np.nanpercentile(im2, percentile))
        v3.append(np.nanpercentile(diff, percentile))
    v1 = np.array(v1)
    v2 = np.array(v2)
    v3 = np.array(v3)

    ax[0].plot(v1, percentiles, 'g-', label ='Image 1')
    ax[0].plot(v2, percentiles, 'r-', label ='Image 2')
    ax[0].plot(v3, percentiles, 'b-', label ='CDS image')
    ax[0].legend(loc = 0)
    ax[0].set(xlabel = 'flux (ADU)', ylabel = 'Percentile', ylim = [0,100])
    ax[1].imshow(diff, vmin = np.min(v3), vmax = np.max(v3),origin = 'lower')
    ax[1].set(xlabel = 'x pix', ylabel = 'ypix', title = 'CDS image')

    outfig = outname.split('.fits')[0]+'.pdf'
    print('\t we write {0}'.format(outfig))
    plt.tight_layout()
    plt.savefig(outfig)

print('\t We write {0} '.format(outname))
fits.writeto(outname, diff, fits.getheader(file2), overwrite = True)
