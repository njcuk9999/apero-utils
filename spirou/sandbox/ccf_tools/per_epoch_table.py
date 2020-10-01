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
    tbl_bin['ERROR_RV'] = np.zeros(len(udates), dtype = float)
    tbl_bin['FORMAL_SIG'] = np.zeros(len(udates), dtype = float)
    tbl_bin['SIG_FROM_MAD'] = np.zeros(len(udates), dtype = float)

    tbl_bin['RV_N'] = np.zeros(len(udates), dtype = int)
    tbl_bin['MJDATE_MEAN'] = np.zeros(len(udates), dtype = float)
    tbl_bin['MJDATE_MIN'] = np.zeros(len(udates), dtype = float)
    tbl_bin['MJDATE_MAX'] = np.zeros(len(udates), dtype = float)

    tbl_bin['RV_MED'] = np.zeros(len(udates), dtype = float)
    tbl_bin['ERROR_RV_MED'] = np.zeros(len(udates), dtype = float)

    tbl['KEEP'] = np.ones_like(tbl,dtype = bool)
    # Table per night
    tbl2 = Table()
    print('\n')
    for i in range(len(udates)):
        g = (udates[i] == tbl['DATE-OBS'])
        weights = 1 / tbl[g]['ERROR_RV'] ** 2
        weights /= np.nansum(weights)

        nMAD = (tbl['RV'][g] - np.nanmedian(tbl['RV'][g]))/np.nanmedian(np.abs(tbl['RV'][g] - np.nanmedian(tbl['RV'][g])))
        nMAD = np.abs(nMAD)


        tbl['KEEP'][g] = nMAD<nMAD_cut

        tbl_bin['RV'][i] = np.sum(weights*tbl['RV'][g])

        # formal sigma for a weighted mean
        tbl_bin['ERROR_RV'][i] = 1/np.sqrt(np.nansum(1/tbl[g]['ERROR_RV']**2))

        # median RV for epoch
        tbl_bin['RV_MED'][i] = np.nanmedian(tbl['RV'][g])

        # median error for points at that epoch
        tbl_bin['ERROR_RV_MED'][i] = np.nanmedian(tbl['ERROR_RV'][g])

        tbl_bin['RV_N'][i] = np.sum(g)

        # normalized to the equivalent of 1 sigma
        tbl_bin['SIG_FROM_MAD'][i] = np.nanmedian(np.abs(tbl['RV'][g]-np.nanmedian(tbl['RV'][g])))/0.681

        tbl_bin['MJDATE_MEAN'][i] = np.mean(tbl['MJDATE'][g])
        tbl_bin['MJDATE_MIN'][i] = np.min(tbl['MJDATE'][g])
        tbl_bin['MJDATE_MAX'][i] = np.max(tbl['MJDATE'][g])

        tbl_bin['FORMAL_SIG'][i] =  tbl_bin['ERROR_RV'][i]/np.sqrt(tbl_bin['RV_N'][i])

        print('For epoch MJD = {0:5.2f}, the biggest outlier is at {1:.3f} median abs deviation'.format(tbl_bin['MJDATE_MEAN'][i],
                                                                                                  np.max(nMAD)))


    print('\nThe binned table has {0} distinct epochs\n'.format(len(tbl_bin)))
    return(tbl_bin)