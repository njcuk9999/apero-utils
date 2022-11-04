#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
APERO fast math functions

Usually replacing a defined numpy function

Created on 2019-09-18 at 10:53

@author: cook
"""
import time
from typing import Tuple, Union

from astropy.io import fits
import bottleneck as bn
import numpy as np
from scipy import signal

# whether we are using bottleneck (should be true always here)
#   we compare to numpy separately
HAS_BOTTLENECK = True
# Whether to recreate vector
MAKE_VECTOR = False

# =============================================================================
# Define variables
# =============================================================================
# alias to non changed nan functions
nanpercentile = np.nanpercentile


# =============================================================================
# Define functions
# =============================================================================
def nanargmax(a: Union[list, np.ndarray],
              axis: Union[None, int, Tuple[int]] = None
              ) -> Union[int, float, np.ndarray]:
    """
    Bottleneck or numpy implementation of nanargmax depending on imports

    :param a: numpy array, Input array. If `a` is not an array, a conversion
              is attempted.
    :param axis: {int, None}, optional, Axis along which the median is computed.
                 The default (axis=None) is to compute the median of the
                 flattened array.

    :type a: np.ndarray
    :type axis: int

    :return: the argument maximum of array `a` (int, float or np.ndarray)
    """
    # set function name
    # _ = display_func('nanargmax', __NAME__)
    # check bottleneck functionality
    if HAS_BOTTLENECK:
        # return bottleneck function
        return bn.nanargmax(a, axis=axis)
    else:
        # return numpy function
        return np.nanargmax(a, axis=axis)


def nanargmin(a: Union[list, np.ndarray],
              axis: Union[None, int, Tuple[int]] = None
              ) -> Union[int, float, np.ndarray]:
    """
    Bottleneck or numpy implementation of nanargmin depending on imports

    :param a: numpy array, Input array. If `a` is not an array, a conversion
              is attempted.
    :param axis: {int, None}, optional, Axis along which the median is computed.
                 The default (axis=None) is to compute the median of the
                 flattened array.

    :type a: np.ndarray
    :type axis: int

    :return: the argument minimum of array `a` (int, float or np.ndarray)
    """
    # set function name
    # _ = display_func('nanargmin', __NAME__)
    # check bottleneck functionality
    if HAS_BOTTLENECK:
        # return bottleneck function
        return bn.nanargmin(a, axis=axis)
    else:
        # return numpy function
        return np.nanargmin(a, axis=axis)


def nanmax(a: Union[list, np.ndarray],
           axis: Union[None, int, Tuple[int]] = None,
           **kwargs) -> Union[int, float, np.ndarray]:
    """
    Bottleneck or numpy implementation of nanmax depending on imports

    :param a: numpy array, Input array. If `a` is not an array, a conversion
              is attempted.
    :param axis: {int, None}, optional, Axis along which the median is computed.
                 The default (axis=None) is to compute the median of the
                 flattened array.
    :param kwargs: keyword arguments passed to numpy function only

    :type a: np.ndarray
    :type axis: int

    :return: the maximum of array `a` (int, float or np.ndarray)
    """
    # set function name
    # _ = display_func('nanmax', __NAME__)
    # check bottleneck functionality
    if HAS_BOTTLENECK and len(kwargs) == 0:
        # return bottleneck function
        return bn.nanmax(a, axis=axis)
    else:
        # return numpy function
        return np.nanmax(a, axis=axis, **kwargs)


def nanmin(a: Union[list, np.ndarray],
           axis: Union[None, int, Tuple[int]] = None,
           **kwargs) -> Union[int, float, np.ndarray]:
    """
    Bottleneck or numpy implementation of nanmin depending on imports

    :param a: numpy array, Input array. If `a` is not an array, a conversion
              is attempted.
    :param axis: {int, None}, optional, Axis along which the median is computed.
                 The default (axis=None) is to compute the median of the
                 flattened array.
    :param kwargs: keyword arguments passed to numpy function only

    :type a: np.ndarray
    :type axis: int

    :return: the minimum of array `a` (int, float or np.ndarray)
    """
    # set function name
    # _ = display_func('nanmin', __NAME__)
    # check bottleneck functionality
    if HAS_BOTTLENECK and len(kwargs) == 0:
        # return bottleneck function
        return bn.nanmin(a, axis=axis)
    else:
        # return numpy function
        return np.nanmin(a, axis=axis, **kwargs)


def nanmean(a: Union[list, np.ndarray],
            axis: Union[None, int, Tuple[int]] = None,
            **kwargs) -> Union[int, float, np.ndarray]:
    """
    Bottleneck or numpy implementation of nanmean depending on imports

    :param a: numpy array, Input array. If `a` is not an array, a conversion
              is attempted.
    :param axis: {int, None}, optional, Axis along which the median is computed.
                 The default (axis=None) is to compute the median of the
                 flattened array.
    :param kwargs: keyword arguments passed to numpy function only

    :type a: np.ndarray
    :type axis: int

    :return: the mean of array `a` (int, float or np.ndarray)
    """
    # set function name
    # _ = display_func('nanmin', __NAME__)
    # check bottleneck functionality
    if HAS_BOTTLENECK and len(kwargs) == 0:
        # return bottleneck function
        return bn.nanmean(a, axis=axis)
    else:
        # return numpy function
        return np.nanmean(a, axis=axis, **kwargs)


def nanmedian(a: Union[list, np.ndarray],
              axis: Union[None, int, Tuple[int]] = None,
              **kwargs) -> Union[int, float, np.ndarray]:
    """
    Bottleneck or numpy implementation of nanmedian depending on imports

    :param a: numpy array, Input array. If `a` is not an array, a conversion
              is attempted.
    :param axis: {int, None}, optional, Axis along which the median is computed.
                 The default (axis=None) is to compute the median of the
                 flattened array.
    :param kwargs: keyword arguments passed to numpy function only

    :type a: np.ndarray
    :type axis: int

    :return: the median of array `a` (int, float or np.ndarray)
    """
    # set function name
    # _ = display_func('nanmedian', __NAME__)
    # check bottleneck functionality
    if HAS_BOTTLENECK and len(kwargs) == 0:
        # return bottleneck function
        return bn.nanmedian(a, axis=axis)
    else:
        # return numpy function
        return np.nanmedian(a, axis=axis, **kwargs)


def nanstd(a: Union[list, np.ndarray],
           axis: Union[None, int, Tuple[int]] = None, ddof: int = 0,
           **kwargs) -> Union[int, float, np.ndarray]:
    """
    Bottleneck or numpy implementation of nanstd depending on imports

    :param a: numpy array, Input array. If `a` is not an array, a conversion
              is attempted.
    :param axis: {int, None}, optional, Axis along which the median is computed.
                 The default (axis=None) is to compute the median of the
                 flattened array.
    :param ddof: int, optional. Means Delta Degrees of Freedom. The divisor
                 used in calculations is ``N - ddof``, where ``N`` represents
                 the number of non-NaN elements. By default `ddof` is zero.
    :param kwargs: keyword arguments passed to numpy function only

    :type a: np.ndarray
    :type axis: int
    :type ddof: int

    :return: the standard deviation of array `a` (int, float or np.ndarray)
    """
    # set function name
    # _ = display_func('nanstd', __NAME__)
    # check bottleneck functionality
    if HAS_BOTTLENECK and len(kwargs) == 0:
        # return bottleneck function
        return bn.nanstd(a, axis=axis, ddof=ddof)
    else:
        # return numpy function
        return np.nanstd(a, axis=axis, ddof=ddof, **kwargs)


def nansum(a: Union[list, np.ndarray],
           axis: Union[None, int, Tuple[int]] = None,
           **kwargs) -> Union[int, float, np.ndarray]:
    """
    Bottleneck or numpy implementation of nansum depending on imports

    :param a: numpy array, Input array. If `a` is not an array, a conversion
              is attempted.
    :param axis: {int, None}, optional, Axis along which the median is computed.
                 The default (axis=None) is to compute the median of the
                 flattened array.
    :param kwargs: keyword arguments passed to numpy function only

    :type a: np.ndarray
    :type axis: int

    :return: the sum of array `a` (int, float or np.ndarray)
    """
    # set function name
    # _ = display_func('nansum', __NAME__)
    # check bottleneck functionality
    if HAS_BOTTLENECK and len(kwargs) == 0:
        # make sure vector a is an array
        if not isinstance(a, np.ndarray):
            a1 = np.array(a)
        else:
            a1 = a
        # bottle neck return in type given for bool array this is not
        #  what we want
        if a1.dtype == bool:
            a1 = a1.astype(int)
        # return bottleneck function
        return bn.nansum(a1, axis=axis)
    else:
        # return numpy function
        return np.nansum(a, axis=axis, **kwargs)


def median(a: Union[list, np.ndarray],
           axis: Union[None, int, Tuple[int]] = None,
           **kwargs) -> Union[int, float, np.ndarray]:
    """
    Bottleneck or numpy implementation of median depending on imports

    :param a: numpy array, Input array. If `a` is not an array, a conversion
              is attempted.
    :param axis: {int, None}, optional, Axis along which the median is computed.
                 The default (axis=None) is to compute the median of the
                 flattened array.
    :param kwargs: keyword arguments passed to numpy function only

    :type a: np.ndarray
    :type axis: int

    :return: the median of array `a` (int, float or np.ndarray)
    """
    # set function name
    # _ = display_func('nansum', __NAME__)
    # check bottleneck functionality
    if HAS_BOTTLENECK and len(kwargs) == 0:
        # return bottleneck function
        return bn.median(a, axis=axis)
    else:
        # return numpy function
        return np.median(a, axis=axis, **kwargs)


def medfilt_1d(a: Union[list, np.ndarray],
               window: Union[None, int] = None) -> np.ndarray:
    """
    Bottleneck or scipy.signal implementation of medfilt depending on imports

    :param a: numpy array, Input array. If `a` is not an array, a conversion
              is attempted.
    :param window: int, The number of elements in the moving window.

    :type a: np.ndarray
    :type window: int

    :return: the 1D median filtered array of `a` (int, float or np.ndarray)
    """
    # set function name
    # _ = display_func('medfilt_1d', __NAME__)
    # check bottleneck functionality
    if HAS_BOTTLENECK:
        # get half window size
        half_window = window // 2
        # need to shift
        a1 = np.append(a, [np.nan] * half_window)
        # median filter (via bottleneck function)
        y = bn.move_median(a1, window=window, min_count=half_window)
        # return shifted bottleneck function
        return y[half_window:]
    else:
        # return scipy function
        return signal.medfilt(a, kernel_size=window)


def test(hdr=None):
    """
    This is where the testing is done

    :param hdr: either a fits Header or None
    :return:
    """
    func_names = ['nanargmax', 'nanargmin', 'nanmax', 'nanmin', 'nanmean',
                  'nanmedian', 'nanstd', 'nansum']

    bottleneck_funcs = [nanargmax, nanargmin, nanmax, nanmin, nanmean,
                        nanmedian, nanstd, nansum]

    numpy_funcs = [np.nanargmax, np.nanargmin, np.nanmax, np.nanmin, np.nanmean,
                   np.nanmedian, np.nanstd, np.nansum]

    storage = dict()

    for it in range(len(bottleneck_funcs)):

        name = func_names[it]
        bfunc = bottleneck_funcs[it]
        nfunc = numpy_funcs[it]

        start1 = time.time()
        bvalue = bfunc(vector)
        end1 = time.time()

        start2 = time.time()
        nvalue = nfunc(vector)
        end2 = time.time()

        time1 = end1 - start1
        time2 = end2 - start2
        speed = time2 / time1
        diff = bvalue - nvalue
        rdiff = diff / bvalue

        print(func_names[it])
        key = f'F{it:02d}'

        if hdr is not None:
            bvaluef = f'\t(comp:{hdr[key + "_BVAL"]:.5e})'
            nvaluef = f'\t(comp:{hdr[key + "_NVAL"]:.5e})'
            difff = f'\t(comp:{hdr[key + "_DIFF"]:.5e})'
            rdifff = f'\t(comp:{hdr[key + "_RDFF"]:.5e})'
        else:
            bvaluef = ''
            nvaluef = ''
            difff = ''
            rdifff = ''

        print(f'\tbottleneck time: {time1:.5e}')
        print(f'\tnumpy time: {time2:.5e}')
        print(f'\tSpeed up: {speed:.5e}')

        print(f'\tBottleneck value: {bvalue:.5e} {bvaluef}')
        print(f'\tNumpy value: {nvalue:.5e} {nvaluef}')
        print(f'\tdiff: {diff:.5e} {difff}')
        print(f'\tdiff/value: {rdiff:.5e} {rdifff}')

        storage[key + '_NAME'] = (name, '{key} = {name}')
        storage[key + '_BVAL'] = (bvalue, f'bottleneck value for {name}')
        storage[key + '_NVAL'] = (nvalue, f'numpy value for {name}')
        storage[key + '_DIFF'] = (diff, f'diff value for {name}')
        storage[key + '_RDFF'] = (rdiff, f'diff/value for {name}')

    return storage


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # define vector
    np.random.seed(1)

    if MAKE_VECTOR:
        size = 4096
        vector = np.random.random(size=size ** 2).reshape((size, size))
        # add nan values
        nanvalues = np.random.randint(size, size=size * 2).reshape(size, 2)
        vector[nanvalues[:, 0], nanvalues[:, 1]] = np.nan
        # add small values
        smallpos = np.random.randint(size, size=size * 2).reshape(size, 2)
        smallvalue = np.random.random(size=size) / 1000
        vector[smallpos] = smallvalue
        # add large values
        largepos = np.random.randint(size, size=size * 2).reshape(size, 2)
        largevalue = np.random.random(size=size) * 1000
        vector[largepos] = largevalue
        hdr = None
    else:
        vector = fits.getdata('fake_vector.fits')
        hdr = fits.getheader('fake_vector.fits')
    # ----------------------------------------------------------------------
    # tests

    _storage = test(hdr)

    if MAKE_VECTOR:

        header = fits.Header()
        for key in _storage:
            header[key] = _storage[key]
        fits.writeto('fake_vector.fits', data=vector, header=header,
                     overwrite=True)

# =============================================================================
# End of code
# =============================================================================
