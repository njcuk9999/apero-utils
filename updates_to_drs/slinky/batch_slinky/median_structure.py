from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# weighted median
def weighted_median(data0, weights0, axis = None):
    """
    :param data: 1D numpy array
    :param weights: 1D numpy array
    :return: float
    """
    # Sort the data

    data = np.array(data0)
    weights =  np.array(weights0)
    bad = ~np.isfinite(data + weights)

    data[bad] = 0
    weights[bad] = 0

    if axis == None:
        data = data.flatten()
        weights = weights.flatten()

        ind_sorted = np.argsort(data)
        data_sorted = data[ind_sorted]
        weights_sorted = weights[ind_sorted]
        # Calculate the cumulative sum of weights
        cumsum = np.cumsum(weights_sorted)
        # Find the median
        median = np.interp(0.5, cumsum / cumsum[-1], data_sorted)

        return median

    if axis == 0:

        median = np.zeros(data.shape[0])
        for ribbon in range(data.shape[0]):
            tmp_data = np.array(data[ribbon])
            tmp_weights = np.array(weights[ribbon])

            ind_sorted = np.argsort(tmp_data)
            data_sorted = tmp_data[ind_sorted]
            weights_sorted = tmp_weights[ind_sorted]
            # Calculate the cumulative sum of weights
            cumsum = np.cumsum(weights_sorted)
            # Find the median

            median[ribbon] = np.interp(0.5, cumsum / cumsum[-1], data_sorted)

        return median

    if axis == 1:

        median = np.zeros(data.shape[1])
        for ribbon in range(data.shape[1]):
            tmp_data = np.array(data[:,ribbon])
            tmp_weights = np.array(weights[:,ribbon])

            ind_sorted = np.argsort(tmp_data)
            data_sorted = tmp_data[ind_sorted]
            weights_sorted = tmp_weights[ind_sorted]
            # Calculate the cumulative sum of weights
            cumsum = np.cumsum(weights_sorted)
            # Find the median

            median[ribbon] = np.interp(0.5, cumsum / cumsum[-1], data_sorted)

        return median



sp = fits.getdata('faint_residual_NIRPS_0_17.fits')

med1 = np.median(sp, axis = 0)

weights = np.ones_like(sp)
for isp in range(sp.shape[0]):
    weights[isp] = np.abs(sp[isp] - med1)

for ite in range(4):
    med2 = np.nanmedian(sp,axis=0)#weighted_median(sp, weights, axis = 1)
    #med2 = weighted_median(sp, weights, axis = 1)
    med2/=np.nansum(med2**2)

    amps = np.zeros(sp.shape[0])
    for isp in range(sp.shape[0]):
        amps[isp] = np.nansum(sp[isp]*med2)
        sp[isp] = sp[isp]/amps[isp]
        weights[isp] = weights[isp]/amps[isp]

    plt.plot(med2)
plt.show()

sp = fits.getdata('faint_residual_NIRPS_0_17.fits')

sp3 = np.array(sp)
for isp in range(sp.shape[0]):
    amp = np.nansum(sp[isp]*med2)
    print(amp)
    sp3[isp] = sp[isp] - med2*amp

fits.writeto('faint_residual_NIRPS_0_17_median.fits',sp3,overwrite=True)