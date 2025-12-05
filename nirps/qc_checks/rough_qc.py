from astropy.io import fits
import numpy as np
import glob

# to be changed to proper directory
files = np.array(glob.glob('/cosmos99/nirps/raw-data/nirps_he/2023-0*/*.fits'))
files = files[::-10]

################################################################################
# constraints
constraints = dict()

# minimum 99th, maximum 99th, max fraction of saturated pixels
constraints['WAVE,FP,FP'] = 10,10000,0.0005
constraints['OBJECT,SKY'] = 0,10000,0.0005
constraints['OBJECT,FP'] = 0,10000,0.0005
constraints['FLUX,STD,SKY'] = 0,10000,0.0005
constraints['ORDERDEF,DARK,LAMP'] = 0,10000,0.0005
constraints['ORDERDEF,LAMP,DARK'] = 0,10000,0.0005
constraints['TELLURIC,SKY'] = 0,10000,0.0005
constraints['DARK'] = 0,1.0,0.0005
constraints['LED,LAMP'] = 100,10000,0.0005
constraints['FLAT,LAMP,DARK'] = 100,10000,0.0005
constraints['FLAT,DARK,LAMP'] = 100,10000,0.0005
constraints['FLAT,LAMP,LAMP'] = 100,10000,0.0005
constraints['CONTAM,DARK,FP'] = 10,10000,0.0005
constraints['WAVE,UN1,UN1'] = 10,100,0.01
constraints['WAVE,FP,UN1'] = 10,100,0.01
constraints['WAVE,UN1,FP'] = 10,100,0.01

################################################################################

for i in range(1, len(files)):
    # fetch header
    h1 = fits.getheader(files[i])
    # get image
    im = fits.getdata(files[i],ext=1)
    # get the number of valid reads
    nread = fits.getdata(files[i],ext=3)

    # get the type of image
    DPRTYPE = h1['HIERARCH ESO DPR TYPE']
    # find the 99th percentile of image
    p99 = np.nanpercentile(im,99)
    # fraction of saturated pixels
    fsat = np.mean(nread != np.max(nread))

    # outputs
    out = '{}/{}\t{}\t{}\t{:5.1f}\t{:.4f}'
    keys = i,len(files), h1['DATE']   ,DPRTYPE,p99,fsat
    print(out.format(*keys))

    # retrieve from constraint dictionary the constraints for that file
    constraint_file = constraints[DPRTYPE]


    err_msg = '\t'+files[i] + ' has the following issues:\n'
    flag_qc = False
    # check against 3 constraints
    if p99 < constraint_file[0]:
        flag_qc = True
        err_msg+='\t\t p99 is too low, limit at {}, value at {:.2f}\n'.format(constraint_file[0], p99)
    if p99 > constraint_file[1]:
        flag_qc = True
        err_msg+='\t\t p99 is too high, limit at {}, value at {:.2f}\n'.format(constraint_file[1], p99)
    if fsat > constraint_file[2]:
        flag_qc = True
        err_msg+='\t\t fsat is too high, limit at {}, value at {:.2f}\n'.format(constraint_file[2], fsat)

    if flag_qc:
        # print error message, to be included in an email at some point
        print(err_msg)

