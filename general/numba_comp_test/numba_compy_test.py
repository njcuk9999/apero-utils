#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare jit vs jit[fastmath=True] vs jit[fastmath=False]

for lin_mini and odd_ratio_mean

Created on 2022-11-01 at 13:32

@author: cook
"""
import time
from typing import Tuple

import numpy as np

# try to import bottleneck module
# noinspection PyBroadException
try:
    import bottleneck as bn

    HAS_BOTTLENECK = True
except Exception as e:
    HAS_BOTTLENECK = False
# try to import numba module
# noinspection PyBroadException
try:
    from numba import jit

    HAS_NUMBA = True
except Exception as _:
    jit = None
    HAS_NUMBA = False


# =============================================================================
# Define functions
# =============================================================================
# Set "nopython" mode for best performance, equivalent to @nji
@jit(nopython=True)
def lin_mini1(vector: np.ndarray, sample: np.ndarray, mm: np.ndarray,
              v: np.ndarray, sz_sample: Tuple[int], case: int,
              recon: np.ndarray, amps: np.ndarray,
              no_recon: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """
    Linear minimization of sample with vector

    Used internally in math.genearl.linear_minimization - you probably should
    use the linear_minization function instead of this directly

    :param vector: vector of N elements
    :param sample: sample: matrix N * M each M column is adjusted in
                   amplitude to minimize the chi2 according to the input vector
    :param mm: zero filled vector for filling size = M
    :param v: zero filled vector for filling size = M
    :param sz_sample: tuple the shape of the sample (N, M)
    :param case: int, if case = 1 then vector.shape[0] = sample.shape[1],
                 if case = 2 vector.shape[0] = sample.shape[0]
    :param recon: zero filled vector size = N, recon output
    :param amps: zero filled vector size = M, amplitudes output
    :param no_recon: boolean if True does not calculate recon
                     (output = input for recon)

    :returns: amps, recon
    """
    # do not set function name here -- cannot use functions here
    # case 1
    if case == 1:
        # fill-in the co-variance matrix
        for i in range(sz_sample[0]):
            for j in range(i, sz_sample[0]):
                mm[i, j] = np.sum(sample[i, :] * sample[j, :])
                # we know the matrix is symetric, we fill the other half
                # of the diagonal directly
                mm[j, i] = mm[i, j]
            # dot-product of vector with sample columns
            v[i] = np.sum(vector * sample[i, :])
        # if the matrix cannot we inverted because the determinant is zero,
        # then we return a NaN for all outputs
        if np.linalg.det(mm) == 0:
            amps = np.zeros(sz_sample[0]) + np.nan
            recon = np.zeros_like(v)
            return amps, recon
        # invert coveriance matrix
        inv = np.linalg.inv(mm)
        # retrieve amplitudes
        for i in range(len(v)):
            for j in range(len(v)):
                amps[i] += inv[i, j] * v[j]
        # reconstruction of the best-fit from the input sample and derived
        # amplitudes
        if not no_recon:
            for i in range(sz_sample[0]):
                recon += amps[i] * sample[i, :]
        return amps, recon
    # same as for case 1 but with axis flipped
    if case == 2:
        # fill-in the co-variance matrix
        for i in range(sz_sample[1]):
            for j in range(i, sz_sample[1]):
                mm[i, j] = np.sum(sample[:, i] * sample[:, j])
                # we know the matrix is symetric, we fill the other half
                # of the diagonal directly
                mm[j, i] = mm[i, j]
            # dot-product of vector with sample columns
            v[i] = np.sum(vector * sample[:, i])
        # if the matrix cannot we inverted because the determinant is zero,
        # then we return a NaN for all outputs
        if np.linalg.det(mm) == 0:
            return amps, recon
        # invert coveriance matrix
        inv = np.linalg.inv(mm)
        # retrieve amplitudes
        for i in range(len(v)):
            for j in range(len(v)):
                amps[i] += inv[i, j] * v[j]
        # reconstruction of the best-fit from the input sample and derived
        # amplitudes
        if not no_recon:
            for i in range(sz_sample[1]):
                recon += amps[i] * sample[:, i]
        return amps, recon


# Set "nopython" mode for best performance, equivalent to @nji
@jit(nopython=True)
def odd_ratio_mean1(value: np.ndarray, error: np.ndarray,
                    odd_ratio: float = 2e-4, nmax: int = 10,
                    conv_cut=1e-2) -> Tuple[float, float]:
    """
    Provide values and corresponding errors and compute a weighted mean

    :param value: np.array (1D), value array
    :param error: np.array (1D), uncertainties for value array
    :param odd_ratio: float, the probability that the point is bad
                    Recommended value in Artigau et al. 2021 : f0 = 0.002
    :param nmax: int, maximum number of iterations to pass through
    :param conv_cut: float, the convergence cut criteria - how precise we have
                     to get

    :return: tuple, 1. the weighted mean, 2. the error on weighted mean
    """
    # deal with NaNs in value or error
    keep = np.isfinite(value) & np.isfinite(error)
    # deal with no finite values
    if np.sum(keep) == 0:
        return np.nan, np.nan
    # remove NaNs from arrays
    value, error = value[keep], error[keep]
    # work out some values to speed up loop
    error2 = error ** 2
    # placeholders for the "while" below
    guess_prev = np.inf
    # the 'guess' must be started as close as we possibly can to the actual
    # value. Starting beyond ~3 sigma (or whatever the odd_ratio implies)
    # would lead to the rejection of pretty much all points and would
    # completely mess the convergence of the loop
    guess = np.nanmedian(value)
    bulk_error = 1.0
    ite = 0
    # loop around until we do all required iterations
    while (np.abs(guess - guess_prev) / bulk_error > conv_cut) and (ite < nmax):
        # store the previous guess
        guess_prev = float(guess)
        # model points as gaussian weighted by likelihood of being a valid point
        # nearly but not exactly one for low-sigma values
        gfit = (1 - odd_ratio) * np.exp(-0.5 * ((value - guess) ** 2 / error2))
        # find the probability that a point is bad
        odd_bad = odd_ratio / (gfit + odd_ratio)
        # find the probability that a point is good
        odd_good = 1 - odd_bad
        # calculate the weights based on the probability of being good
        weights = odd_good / error2
        # update the guess based on the weights
        if np.sum(np.isfinite(weights)) == 0:
            guess = np.nan
        else:
            guess = np.nansum(value * weights) / np.nansum(weights)
            # work out the bulk error
            bulk_error = np.sqrt(1.0 / np.nansum(odd_good / error2))
        # keep track of the number of iterations
        ite += 1

    # return the guess and bulk error
    return guess, bulk_error


# Set "nopython" mode for best performance, equivalent to @nji
@jit(nopython=True, fastmath=True)
def lin_mini2(vector: np.ndarray, sample: np.ndarray, mm: np.ndarray,
              v: np.ndarray, sz_sample: Tuple[int], case: int,
              recon: np.ndarray, amps: np.ndarray,
              no_recon: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """
    Linear minimization of sample with vector

    Used internally in math.genearl.linear_minimization - you probably should
    use the linear_minization function instead of this directly

    :param vector: vector of N elements
    :param sample: sample: matrix N * M each M column is adjusted in
                   amplitude to minimize the chi2 according to the input vector
    :param mm: zero filled vector for filling size = M
    :param v: zero filled vector for filling size = M
    :param sz_sample: tuple the shape of the sample (N, M)
    :param case: int, if case = 1 then vector.shape[0] = sample.shape[1],
                 if case = 2 vector.shape[0] = sample.shape[0]
    :param recon: zero filled vector size = N, recon output
    :param amps: zero filled vector size = M, amplitudes output
    :param no_recon: boolean if True does not calculate recon
                     (output = input for recon)

    :returns: amps, recon
    """
    # do not set function name here -- cannot use functions here
    # case 1
    if case == 1:
        # fill-in the co-variance matrix
        for i in range(sz_sample[0]):
            for j in range(i, sz_sample[0]):
                mm[i, j] = np.sum(sample[i, :] * sample[j, :])
                # we know the matrix is symetric, we fill the other half
                # of the diagonal directly
                mm[j, i] = mm[i, j]
            # dot-product of vector with sample columns
            v[i] = np.sum(vector * sample[i, :])
        # if the matrix cannot we inverted because the determinant is zero,
        # then we return a NaN for all outputs
        if np.linalg.det(mm) == 0:
            amps = np.zeros(sz_sample[0]) + np.nan
            recon = np.zeros_like(v)
            return amps, recon
        # invert coveriance matrix
        inv = np.linalg.inv(mm)
        # retrieve amplitudes
        for i in range(len(v)):
            for j in range(len(v)):
                amps[i] += inv[i, j] * v[j]
        # reconstruction of the best-fit from the input sample and derived
        # amplitudes
        if not no_recon:
            for i in range(sz_sample[0]):
                recon += amps[i] * sample[i, :]
        return amps, recon
    # same as for case 1 but with axis flipped
    if case == 2:
        # fill-in the co-variance matrix
        for i in range(sz_sample[1]):
            for j in range(i, sz_sample[1]):
                mm[i, j] = np.sum(sample[:, i] * sample[:, j])
                # we know the matrix is symetric, we fill the other half
                # of the diagonal directly
                mm[j, i] = mm[i, j]
            # dot-product of vector with sample columns
            v[i] = np.sum(vector * sample[:, i])
        # if the matrix cannot we inverted because the determinant is zero,
        # then we return a NaN for all outputs
        if np.linalg.det(mm) == 0:
            return amps, recon
        # invert coveriance matrix
        inv = np.linalg.inv(mm)
        # retrieve amplitudes
        for i in range(len(v)):
            for j in range(len(v)):
                amps[i] += inv[i, j] * v[j]
        # reconstruction of the best-fit from the input sample and derived
        # amplitudes
        if not no_recon:
            for i in range(sz_sample[1]):
                recon += amps[i] * sample[:, i]
        return amps, recon


# Set "nopython" mode for best performance, equivalent to @nji
@jit(nopython=True, fastmath=True)
def odd_ratio_mean2(value: np.ndarray, error: np.ndarray,
                    odd_ratio: float = 2e-4, nmax: int = 10,
                    conv_cut=1e-2) -> Tuple[float, float]:
    """
    Provide values and corresponding errors and compute a weighted mean

    :param value: np.array (1D), value array
    :param error: np.array (1D), uncertainties for value array
    :param odd_ratio: float, the probability that the point is bad
                    Recommended value in Artigau et al. 2021 : f0 = 0.002
    :param nmax: int, maximum number of iterations to pass through
    :param conv_cut: float, the convergence cut criteria - how precise we have
                     to get

    :return: tuple, 1. the weighted mean, 2. the error on weighted mean
    """
    # deal with NaNs in value or error
    keep = np.isfinite(value) & np.isfinite(error)
    # deal with no finite values
    if np.sum(keep) == 0:
        return np.nan, np.nan
    # remove NaNs from arrays
    value, error = value[keep], error[keep]
    # work out some values to speed up loop
    error2 = error ** 2
    # placeholders for the "while" below
    guess_prev = np.inf
    # the 'guess' must be started as close as we possibly can to the actual
    # value. Starting beyond ~3 sigma (or whatever the odd_ratio implies)
    # would lead to the rejection of pretty much all points and would
    # completely mess the convergence of the loop
    guess = np.nanmedian(value)
    bulk_error = 1.0
    ite = 0
    # loop around until we do all required iterations
    while (np.abs(guess - guess_prev) / bulk_error > conv_cut) and (ite < nmax):
        # store the previous guess
        guess_prev = float(guess)
        # model points as gaussian weighted by likelihood of being a valid point
        # nearly but not exactly one for low-sigma values
        gfit = (1 - odd_ratio) * np.exp(-0.5 * ((value - guess) ** 2 / error2))
        # find the probability that a point is bad
        odd_bad = odd_ratio / (gfit + odd_ratio)
        # find the probability that a point is good
        odd_good = 1 - odd_bad
        # calculate the weights based on the probability of being good
        weights = odd_good / error2
        # update the guess based on the weights
        if np.sum(np.isfinite(weights)) == 0:
            guess = np.nan
        else:
            guess = np.nansum(value * weights) / np.nansum(weights)
            # work out the bulk error
            bulk_error = np.sqrt(1.0 / np.nansum(odd_good / error2))
        # keep track of the number of iterations
        ite += 1

    # return the guess and bulk error
    return guess, bulk_error


# Set "nopython" mode for best performance, equivalent to @nji
@jit(nopython=True, fastmath=False)
def lin_mini3(vector: np.ndarray, sample: np.ndarray, mm: np.ndarray,
              v: np.ndarray, sz_sample: Tuple[int], case: int,
              recon: np.ndarray, amps: np.ndarray,
              no_recon: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """
    Linear minimization of sample with vector

    Used internally in math.genearl.linear_minimization - you probably should
    use the linear_minization function instead of this directly

    :param vector: vector of N elements
    :param sample: sample: matrix N * M each M column is adjusted in
                   amplitude to minimize the chi2 according to the input vector
    :param mm: zero filled vector for filling size = M
    :param v: zero filled vector for filling size = M
    :param sz_sample: tuple the shape of the sample (N, M)
    :param case: int, if case = 1 then vector.shape[0] = sample.shape[1],
                 if case = 2 vector.shape[0] = sample.shape[0]
    :param recon: zero filled vector size = N, recon output
    :param amps: zero filled vector size = M, amplitudes output
    :param no_recon: boolean if True does not calculate recon
                     (output = input for recon)

    :returns: amps, recon
    """
    # do not set function name here -- cannot use functions here
    # case 1
    if case == 1:
        # fill-in the co-variance matrix
        for i in range(sz_sample[0]):
            for j in range(i, sz_sample[0]):
                mm[i, j] = np.sum(sample[i, :] * sample[j, :])
                # we know the matrix is symetric, we fill the other half
                # of the diagonal directly
                mm[j, i] = mm[i, j]
            # dot-product of vector with sample columns
            v[i] = np.sum(vector * sample[i, :])
        # if the matrix cannot we inverted because the determinant is zero,
        # then we return a NaN for all outputs
        if np.linalg.det(mm) == 0:
            amps = np.zeros(sz_sample[0]) + np.nan
            recon = np.zeros_like(v)
            return amps, recon
        # invert coveriance matrix
        inv = np.linalg.inv(mm)
        # retrieve amplitudes
        for i in range(len(v)):
            for j in range(len(v)):
                amps[i] += inv[i, j] * v[j]
        # reconstruction of the best-fit from the input sample and derived
        # amplitudes
        if not no_recon:
            for i in range(sz_sample[0]):
                recon += amps[i] * sample[i, :]
        return amps, recon
    # same as for case 1 but with axis flipped
    if case == 2:
        # fill-in the co-variance matrix
        for i in range(sz_sample[1]):
            for j in range(i, sz_sample[1]):
                mm[i, j] = np.sum(sample[:, i] * sample[:, j])
                # we know the matrix is symetric, we fill the other half
                # of the diagonal directly
                mm[j, i] = mm[i, j]
            # dot-product of vector with sample columns
            v[i] = np.sum(vector * sample[:, i])
        # if the matrix cannot we inverted because the determinant is zero,
        # then we return a NaN for all outputs
        if np.linalg.det(mm) == 0:
            return amps, recon
        # invert coveriance matrix
        inv = np.linalg.inv(mm)
        # retrieve amplitudes
        for i in range(len(v)):
            for j in range(len(v)):
                amps[i] += inv[i, j] * v[j]
        # reconstruction of the best-fit from the input sample and derived
        # amplitudes
        if not no_recon:
            for i in range(sz_sample[1]):
                recon += amps[i] * sample[:, i]
        return amps, recon


# Set "nopython" mode for best performance, equivalent to @nji
@jit(nopython=True, fastmath=False)
def odd_ratio_mean3(value: np.ndarray, error: np.ndarray,
                    odd_ratio: float = 2e-4, nmax: int = 10,
                    conv_cut=1e-2) -> Tuple[float, float]:
    """
    Provide values and corresponding errors and compute a weighted mean

    :param value: np.array (1D), value array
    :param error: np.array (1D), uncertainties for value array
    :param odd_ratio: float, the probability that the point is bad
                    Recommended value in Artigau et al. 2021 : f0 = 0.002
    :param nmax: int, maximum number of iterations to pass through
    :param conv_cut: float, the convergence cut criteria - how precise we have
                     to get

    :return: tuple, 1. the weighted mean, 2. the error on weighted mean
    """
    # deal with NaNs in value or error
    keep = np.isfinite(value) & np.isfinite(error)
    # deal with no finite values
    if np.sum(keep) == 0:
        return np.nan, np.nan
    # remove NaNs from arrays
    value, error = value[keep], error[keep]
    # work out some values to speed up loop
    error2 = error ** 2
    # placeholders for the "while" below
    guess_prev = np.inf
    # the 'guess' must be started as close as we possibly can to the actual
    # value. Starting beyond ~3 sigma (or whatever the odd_ratio implies)
    # would lead to the rejection of pretty much all points and would
    # completely mess the convergence of the loop
    guess = np.nanmedian(value)
    bulk_error = 1.0
    ite = 0
    # loop around until we do all required iterations
    while (np.abs(guess - guess_prev) / bulk_error > conv_cut) and (ite < nmax):
        # store the previous guess
        guess_prev = float(guess)
        # model points as gaussian weighted by likelihood of being a valid point
        # nearly but not exactly one for low-sigma values
        gfit = (1 - odd_ratio) * np.exp(-0.5 * ((value - guess) ** 2 / error2))
        # find the probability that a point is bad
        odd_bad = odd_ratio / (gfit + odd_ratio)
        # find the probability that a point is good
        odd_good = 1 - odd_bad
        # calculate the weights based on the probability of being good
        weights = odd_good / error2
        # update the guess based on the weights
        if np.sum(np.isfinite(weights)) == 0:
            guess = np.nan
        else:
            guess = np.nansum(value * weights) / np.nansum(weights)
            # work out the bulk error
            bulk_error = np.sqrt(1.0 / np.nansum(odd_good / error2))
        # keep track of the number of iterations
        ite += 1

    # return the guess and bulk error
    return guess, bulk_error


def linear_minimization(vector: np.ndarray, sample: np.ndarray,
                        no_recon: bool = False, linfunc=None
                        ) -> Tuple[np.ndarray, np.ndarray]:
    """
    wrapper function that sets everything for the @jit later
    In particular, we avoid the np.zeros that are not handled
    by numba, size of input vectors and sample to be adjusted

    :param vector: 2d matrix that is (N x M) or (M x N)
    :param sample: 1d vector of length N
    :param no_recon: bool, if True does not calculate recon
    :param linfunc: special one time thing where we pass a lin_mini code
                    to here
    :return:
    """
    # set function name
    func_name = 'linear_minimization()'
    # get sample and vector shapes
    sz_sample = sample.shape  # 1d vector of length N
    sz_vector = vector.shape  # 2d matrix that is N x M or M x N
    # define which way the sample is flipped relative to the input vector
    if sz_vector[0] == sz_sample[0]:
        case = 2
    elif sz_vector[0] == sz_sample[1]:
        case = 1
    else:
        emsg = ('Neither vector[0]==sample[0] nor vector[0]==sample[1] '
                '(function = {0})')
        print(emsg)
        raise ValueError(emsg.format(func_name))
    # ----------------------------------------------------------------------
    # Part A) we deal with NaNs
    # ----------------------------------------------------------------------
    # set up keep vector
    keep = None
    # we check if there are NaNs in the vector or the sample
    # if there are NaNs, we'll fit the rest of the domain
    isnan = (np.sum(np.isnan(vector)) != 0) or (np.sum(np.isnan(sample)) != 0)
    # ----------------------------------------------------------------------
    # case 1: sample is not flipped relative to the input vector
    if case == 1:
        if isnan:
            # we create a mask of non-NaN
            keep = np.isfinite(vector) * np.isfinite(np.sum(sample, axis=0))
            # redefine the input vector to avoid NaNs
            vector = vector[keep]
            sample = sample[:, keep]
            # re-find shapes
            sz_sample = sample.shape  # 1d vector of length N
        # matrix of covariances
        mm = np.zeros([sz_sample[0], sz_sample[0]])
        # cross-terms of vector and columns of sample
        vec = np.zeros(sz_sample[0])
        # reconstructed amplitudes
        amps = np.zeros(sz_sample[0])
        # reconstruted fit
        recon = np.zeros(sz_sample[1])
    # ----------------------------------------------------------------------
    # case 2: sample is flipped relative to the input vector
    elif case == 2:
        # same as for case 1, but with axis flipped
        if isnan:
            # we create a mask of non-NaN
            keep = np.isfinite(vector) * np.isfinite(np.sum(sample, axis=1))
            vector = vector[keep]
            sample = sample[keep, :]
            # re-find shapes
            sz_sample = sample.shape  # 1d vector of length N
        mm = np.zeros([sz_sample[1], sz_sample[1]])
        vec = np.zeros(sz_sample[1])
        amps = np.zeros(sz_sample[1])
        recon = np.zeros(sz_sample[0])
    # ----------------------------------------------------------------------
    # should not get here (so just repeat the raise from earlier)
    else:
        emsg = ('Neither vector[0]==sample[0] nor vector[0]==sample[1] '
                '(function = {0})')
        raise ValueError(emsg.format(func_name))

    # ----------------------------------------------------------------------
    # Part B) pass to optimized linear minimization
    # ----------------------------------------------------------------------
    # pass all variables and pre-formatted vectors to the @jit part of the code
    amp_out, recon_out = linfunc(vector, sample, mm, vec, sz_sample,
                                 case, recon, amps, no_recon=no_recon)
    # ----------------------------------------------------------------------
    # if we had NaNs in the first place, we create a reconstructed vector
    # that has the same size as the input vector, but pad with NaNs values
    # for which we cannot derive a value
    if isnan:
        recon_out2 = np.zeros_like(keep) + np.nan
        recon_out2[keep] = recon_out
        recon_out = recon_out2

    return amp_out, recon_out


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # TODO: @etienne fill these in with real values
    npix = 1000
    nvec = 7
    # length of input vector
    vector = np.arange(npix)
    # sample of dummy vectors for linear minimization
    sample = np.random.normal(size=[nvec, npix]) + 1.0
    # dummy bad
    bad = np.random.normal(size=[nvec, npix]) > 3  # more than 3-sigma -> NaN
    sample[bad] = np.nan

    nvalue = 1e5
    value = np.random.normal(size=int(nvalue))
    # errors are between 0.5 and 1.5
    err = 1 + np.random.random(size=int(nvalue)) * 0.5
    bad = (np.random.normal(size=int(nvalue)) > 3)
    value[bad] = np.nan
    err[bad] = np.nan
    odd_ratio = 1e-5


    lin_mini_kwargs = dict(vector=vector, sample=sample)
    odd_ratio_mean_args = dict(value=value, error=err, odd_ratio=odd_ratio)

    # ----------------------------------------------------------------------
    lfunc = [lin_mini1, lin_mini2, lin_mini3]
    lname = ['lin_mini1', 'lin_mini2', 'lin_mini3']
    lout = dict()
    # test these
    for name, func in zip(lname, lfunc):
        # run once to build
        _ = linear_minimization(linfunc=func, **lin_mini_kwargs)
        # run timed with results
        start = time.time()
        lout[name] = linear_minimization(linfunc=func, **lin_mini_kwargs)
        end = time.time()
        print(f'{name} took: {end - start} s')

    ofunc = [odd_ratio_mean1, odd_ratio_mean2, odd_ratio_mean3]
    oname = ['odd_ratio_mean1', 'odd_ratio_mean2', 'odd_ratio_mean3']
    oout = dict()
    for name, func in zip(oname, ofunc):
        # run once to build
        try:
            _ = func(**odd_ratio_mean_args)
        except:
            pass
        # run timed with results
        start = time.time()
        try:
            oout[name] = func(**odd_ratio_mean_args)
        except Exception as e:
            oout[name] = np.nan, np.nan
            print(f'ERROR {type(e)}: {str(e)}')
        end = time.time()
        print(f'{name} took: {end - start} s')
    # ----------------------------------------------------------------------
    # TODO: @etienne results for lin mini
    # ----------------------------------------------------------------------
    name_ref = 'lin_mini1'
    for name in lout:
        if name == name_ref:
            continue
        # get difference and ratio compared to "name_ref"
        diff_amp = lout[name][0] - lout[name_ref][0]
        diff_recon = lout[name][1] - lout[name_ref][1]
        # ratio_amp = lout[name][0] / lout[name_ref][0]
        # ratio_recon = lout[name][1] / lout[name_ref][1]
        #
        print(f'{name_ref}-{name} Amp std: {np.nanstd(diff_amp)}   Recon std: {np.nanstd(diff_recon)}')


    # ----------------------------------------------------------------------
    # TODO: @etienne results for lin mini
    # ----------------------------------------------------------------------
    name_ref = 'odd_ratio_mean1'
    for name in oout:
        # get difference and ratio compared to "name_ref"
        guess = oout[name][0]
        bulk_error = oout[name][1]
        # ratio_guess = oout[name][0] / oout[name_ref][0]
        # ratio_bulk_error = oout[name][1] / oout[name_ref][1]
        print(f'{name_ref}-{name}:   {guess}+/-{bulk_error}')


# =============================================================================
# End of code
# =============================================================================
