import numpy as np
from astropy.table import Table

def per_epoch_table(tbl_input,nMAD_cut = np.inf):
    #
    # Takes a table with many visits for each epoch
    # and computes stats to get a per-epoch mean velocity
    # One can provide a maximum median absolute deviation
    # to reject outliers. Reasonable values for the
    # nMAD_cut range from 3 (agressive clipping)
    # to 7 (pretty moderate clipping)
    #
    tbl = Table(tbl_input)

    udates = np.unique(tbl['DATE-OBS'])

    # Table with epoch-binned data
    tbl_bin = Table()

    tbl_bin['RV'] = np.zeros(len(udates), dtype = float)
    tbl_bin['RV_SIG'] = np.zeros(len(udates), dtype = float)
    tbl_bin['FORMAL_SIG'] = np.zeros(len(udates), dtype = float)
    tbl_bin['SIG_FROM_MAD'] = np.zeros(len(udates), dtype = float)

    tbl_bin['RV_N'] = np.zeros(len(udates), dtype = int)
    tbl_bin['MJDATE_MEAN'] = np.zeros(len(udates), dtype = float)
    tbl_bin['MJDATE_MIN'] = np.zeros(len(udates), dtype = float)
    tbl_bin['MJDATE_MAX'] = np.zeros(len(udates), dtype = float)

    tbl['KEEP'] = np.ones_like(tbl,dtype = bool)
    # Table per night
    tbl2 = Table()
    for i in range(len(udates)):
        g = (udates[i] == tbl['DATE-OBS'])

        nMAD = (tbl['RV'][g] - np.nanmedian(tbl['RV'][g]))/np.nanmedian(np.abs(tbl['RV'][g] - np.nanmedian(tbl['RV'][g])))
        nMAD = np.abs(nMAD)
        print(np.max(nMAD))

        tbl['KEEP'][g] = nMAD<nMAD_cut

        tbl_bin['RV'][i] = np.mean(tbl['RV'][g])
        tbl_bin['RV_SIG'][i] = np.std(tbl['RV'][g])
        tbl_bin['RV_N'][i] = np.sum(g)

        # normalized to the equivalent of 1 sigma
        tbl_bin['SIG_FROM_MAD'][i] = np.nanmedian(np.abs(tbl['RV'][g]-np.nanmedian(tbl['RV'][g])))/0.681

        tbl_bin['MJDATE_MEAN'][i] = np.mean(tbl['MJDATE'][g])
        tbl_bin['MJDATE_MIN'][i] = np.min(tbl['MJDATE'][g])
        tbl_bin['MJDATE_MAX'][i] = np.max(tbl['MJDATE'][g])

        tbl_bin['FORMAL_SIG'][i] =  tbl_bin['RV_SIG'][i]/np.sqrt(tbl_bin['RV_N'][i])

    print(len(tbl))
    tbl = tbl[tbl['KEEP']]

    return(tbl_bin)