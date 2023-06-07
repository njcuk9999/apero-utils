#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_apero_eso.py

Compare APERO and ESO radial velocities using LBL or CCF.

@author: lim
"""
# Imports
import os.path
import matplotlib.pyplot as plt
import glob
from astropy.io import ascii, fits
import h5py
from astropy.time import Time
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from dace_query.spectroscopy import Spectroscopy

# For uniform fonts in plots
fontsize = 20
markersize = 12
capsize = 2
plt.rcParams["font.size"] = fontsize
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.family"] = "STIXGeneral"
plt.rcParams["axes.linewidth"] = 1
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.major.width"] = 1
plt.rcParams["ytick.major.width"] = 1
plt.rcParams["xtick.major.size"] = 4
plt.rcParams["ytick.major.size"] = 4
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True
# plt.style.use('dark_background')


# =====================================================================================================================
# Define variables
# =====================================================================================================================
# APERO reduction info
mode = "he"  # "ha": high accuracy | "he": high efficiency
version = "07276"  # version number I assume?
onoff = "online"  # "online" | "offline"

# Paths to LBL rdb or CCF files
# path_1 = "/cosmos99/nirps/lbl-data/nirps_{}_{}_{}/lblrdb/".format(mode, version, onoff)  # APERO LBL
path_1 = "/cosmos99/nirps/apero-data/nirps_{}_{}_{}/red/".format(mode, version, onoff)  # APERO CCF
path_2 = ''  # ESO (this path is not used; the code will query DACE for CCF, and crash for LBL)
# path_2 = "/cosmos99/nirps/apero-data/nirps_{}_{}_{}/red/".format(mode, version, onoff)  # APERO CCF

# Target & template
target = "Proxima"  # name of target
template = "Proxima"  # name of template; only used if on of the two dataset is from LBL

# Pipelines to compare (APERO/ESO)
pipeline_1 = "APERO"
pipeline_2 = "ESO"

# RV methods to compare (LBL/CCF)
rv_method_1 = "CCF"
rv_method_2 = "CCF"

# Path to save figures (and ccf dictionary)
path_savefig = ""  # TODO Where?
# apero_ccf_dict_fname = "apero_ccf_dict.h5"

# Number of sigmas to reject outliers
nsig = 10


# =====================================================================================================================
# Define functions
# =====================================================================================================================
def loadhdf5(filename):
    """
    Load dictionary or object from hdf5 file.

    :param filename:    (str) Name of the hdf5 file.
    :return:            obj
    """
    obj = dict()
    hf = h5py.File(filename, 'r')
    for key in hf.keys():
        obj[key] = hf[key][()]
    hf.close()
    return obj


def savehdf5(obj, filename, verbose=False):
    """
    Save dictionary or object to machine-readable hdf5 file.

    :param obj:         (dict or object) Dictionary or object to save.
    :param filename:    (str) Name of the hdf5 file.
    :param verbose:     (bool, optional) Set to True to print status as the dictionary/object is being saved.
                        Default: False.
    :return:
    """
    hf = h5py.File(filename,'w')
    for key in obj.keys():
        try:
            if type(obj[key]) is dict:
                grp = hf.create_group(key)
                if verbose:
                    print('creating group: '+str(key))
                for keyLevel2 in obj[key].keys():
                    try:
                        if type(obj[key][keyLevel2]) is dict:
                            print('key = '+key+'/'+keyLevel2 +'is a dict and will not be saved.')
                        else:
                            grp.create_dataset("keyLevel2",data=obj[key][keyLevel2])
                        if verbose:
                            print(key+'/'+keyLevel2)
                    except:
                        print('Element could not be save: key = '+key+'/'+keyLevel2)
            else:
                hf.create_dataset(key,data=obj[key])
                if verbose:
                    print(key)
        except:
            print('Element could not be save: key = '+str(key))
    hf.close()
    return hf


def sigma(im):
    """"""  # TODO
    return (np.nanpercentile(im, 85) - np.nanpercentile(im, 15)) / 2


def lowpassfilter(input_vect, width=101):
    """"""  # TODO
    # Computes a low-pass filter of an input vector. This is done while properly handling
    # NaN values, but at the same time being reasonably fast.
    # Algorithm:
    #
    # provide an input vector of an arbtrary length and compute a running NaN median over a
    # box of a given length (width value). The running median is NOT computed at every pixel
    # but at steps of 1/4th of the width value. This provides a vector of points where
    # the nan-median has been computed (ymed) and mean position along the input vector (xmed)
    # of valid (non-NaN) pixels. This xmed/ymed combination is then used in a spline to
    # recover a vector for all pixel positions within the input vector.
    #
    # When there are no valid pixel in a 'width' domain, the value is skipped in the creation
    # of xmed and ymed, and the domain is splined over.

    # indices along input vector
    index = np.arange(len(input_vect))

    # placeholders for x and y position along vector
    xmed = []
    ymed = []

    # loop through the lenght of the input vector
    for i in np.arange(-width // 2, len(input_vect) + width // 2, width // 4):

        # if we are at the start or end of vector, we go 'off the edge' and
        # define a box that goes beyond it. It will lead to an effectively
        # smaller 'width' value, but will provide a consistent result at edges.
        low_bound = i
        high_bound = i + int(width)

        if low_bound < 0:
            low_bound = 0
        if high_bound > (len(input_vect) - 1):
            high_bound = (len(input_vect) - 1)

        pixval = index[low_bound:high_bound]

        if len(pixval) < 3:
            continue

        # if no finite value, skip
        if np.max(np.isfinite(input_vect[pixval])) == 0:
            continue

        # mean position along vector and NaN median value of
        # points at those positions
        xmed.append(np.nanmean(pixval))
        ymed.append(np.nanmedian(input_vect[pixval]))

    xmed = np.array(xmed, dtype=float)
    ymed = np.array(ymed, dtype=float)

    # we need at least 3 valid points to return a
    # low-passed vector.
    if len(xmed) < 3:
        return np.zeros_like(input_vect) + np.nan

    if len(xmed) != len(np.unique(xmed)):
        xmed2 = np.unique(xmed)
        ymed2 = np.zeros_like(xmed2)
        for i in range(len(xmed2)):
            ymed2[i] = np.mean(ymed[xmed == xmed2[i]])
        xmed = xmed2
        ymed = ymed2

    # splining the vector
    spline = ius(xmed, ymed, k=1, ext=3)
    lowpass = spline(np.arange(len(input_vect)))

    return lowpass


def compute_ptpscatter(v=None, method='linear_sigma'):
    """"""  # TODO
    if v is None:
        return ['naive_sigma', 'linear_sigma', 'lowpass_sigma', 'quadratic_sigma']

    if method == 'naive_sigma':
        # we don't know that there's a transit
        return (np.nanpercentile(v, 84) - np.nanpercentile(v, 16)) / 2.0

    if method == 'linear_sigma':
        # difference to immediate neighbours
        return sigma(v[1:-1] - (v[2:] + v[0:-2]) / 2) / np.sqrt(1 + 0.5)
    if method == 'lowpass_sigma':
        return sigma(v - lowpassfilter(v, width=15))

    if method == 'quadratic_sigma':
        v2 = -np.roll(v, -2) / 3 + np.roll(v, -1) + np.roll(v, 1) / 3
        return sigma(v - v2) / np.sqrt(20 / 9.0)


def get_target_template(fname):
    """
    Get target name and template name from file name.
    :param fname:   (str) File name.
    :return:        target (str), template (str)
    """
    target_template = fname.split('/')[-1][:-len(".rdb")]  # should look something like lbl_TARGETNAME_TEMPLATENAME
    target_template = target_template.split('_')[1:]  # should look something like [TARGETNAME, TEMPLATENAME]
    target = target_template[0]
    template = target_template[1]
    return target, template


def open_rdb(fname):
    """
    Load LBL rdb file to get JD, velocity, and velocity error bar.

    :param fname:   (str) Absolute path to the rdb file.
    :return:        rjd (1d array), rv (1d array), rv_err (1d array)
    """
    # Open file
    tbl_i = ascii.read(fname)

    # Get date, vrad, svrad
    rjd, rv, rv_err = tbl_i["rjd"], tbl_i["vrad"], tbl_i["svrad"]

    return rjd, rv, rv_err


def query_dace(target, mode="HE"):
    """
    Query DACE to get RV time series
    :param target:  (str) Target name
    :param mode:    (str, optional) "HE" for "high efficiency", "HA" for "high accuracy". Default: "HE".
    :return:        rdj (1d array), rv (1d array), rv_err (1d array)
    """
    # Retrieve radial velocity timeseries from DACE
    dict_dace = Spectroscopy.get_timeseries(target=target, sorted_by_instrument=True, output_format='numpy')
    dict_dace = dict_dace["NIRPS"]["3.0.0"][mode.upper()]

    # Get date, rv, rv_err
    rjd, rv, rv_err = dict_dace["rjd"], dict_dace["rv"], dict_dace["rv_err"]

    return rjd, rv, rv_err


def process_rvs(rjd, rv, rv_err, pipeline, rv_method, nsig=10):
    """
    Compute median, standard deviation, and median error bar, and remove outliers.

    :param rjd:         (1d array) Reduced Julian day array.
    :param rv:          (1d array) Radial velocity array.
    :param rv_err:      (1d array) Radial velocity uncertainty array.
    :param pipeline:    (str) "APERO" or "ESO".
    :param rv_method    (str) "LBL" or "CCF".
    :param nsig:        (float, optional) Number of sigmas away from median velocity to reject outliers. Default: 10.
    :return:            rv_dict (dict)
    """
    # Convert date from RJD to JD
    jd = rjd + 2.4e6  # [RJD -> JD]
    # date_apero = Time(jdmid, format="jd").to_datetime()

    # Compute std, median, and median error bar before outlier rejection
    std_full = np.nanpercentile(rv, q=[16, 84])
    std_full = (std_full[-1] - std_full[0]) / 2
    if std_full == 0.:
        assert rv_err.size == 1
        std_full = rv_err.copy()
    median_full = np.nanmedian(rv)
    mederr_full = np.nanmedian(rv_err)
    ptpscat_full = compute_ptpscatter(rv, method="linear_sigma")

    # Make a copy of the full time series
    jd_full, rv_full, rv_err_full = jd.copy(), rv.copy(), rv_err.copy()

    # Remove nsig * sigma outliers
    # TODO: Median-filter before finding outliers?
    # TODO: Compare with how outliers are removed in ARI plots.
    # inonout = np.where(np.abs(rv - median_full) / rv_err < nsig)  # non-outliers indices
    inonout = np.where((rv_err < 50 * mederr_full) * (np.abs(rv - median_full) / std_full < nsig))
    # iout_up = np.where((rv - median_full) / rv_err > nsig)  # outliers indices (> nsig)
    iout_up = np.where((rv_err > 50 * mederr_full) * (rv > median_full) + ((rv - median_full) / std_full > nsig))
    # iout_lo = np.where((rv - median_full) / rv_err < -nsig)  # outliers indices (< -nsig)
    iout_lo = np.where((rv_err > 50 * mederr_full) * (rv < median_full) + ((rv - median_full) / std_full < -nsig))
    jd_out_up = jd[iout_up]
    jd_out_lo = jd[iout_lo]
    jd, rv, rv_err = jd[inonout], rv[inonout], rv_err[inonout]

    # Compute std, median and median error bar, after outlier rejection
    std = np.nanstd(rv)
    median = np.nanmedian(rv)
    mederr = np.nanmedian(rv_err)
    ptpscat = compute_ptpscatter(rv, method="linear_sigma")

    # Make dictionary to store all arrays
    rv_dict = {"full": {"jd": jd_full,
                        "rv": rv_full,
                        "rv_err": rv_err_full,
                        "std": std_full,
                        "median": median_full,
                        "mederr": mederr_full,
                        "ptpscat": ptpscat_full
                        },
               "no_outliers": {"jd": jd,
                               "rv": rv,
                               "rv_err": rv_err,
                               "std": std,
                               "median": median,
                               "mederr": mederr,
                               "indices": inonout,
                               "ptpscat": ptpscat
                               },
               "outliers_up": {"jd": jd_out_up,
                               # "rv": rv_out_up,  # obsolete
                               "indices": iout_up},
               "outliers_lo": {"jd": jd_out_lo,
                               # "rv": rv_out_lo,  # obsolete
                               "indices": iout_lo},
               "pipeline": pipeline,
               "rv_method": rv_method
               }

    return rv_dict


def compare_lbl(path_apero, nsig=10, path_savefig=''):
    """
    OBSOLETE
    Compare APERO LBL radial velocities to ESO CCF radial velocities.

    :param path_apero:      (str) Absolute path to APERO LBL rdb folder.
    :param nsig:            (float, optional) Number of sigmas away from median velocity to reject outliers. Default: 10.
    :param path_savefig:    (str, optional) Path where to save figures. Default: ''.
    :return: 0
    """
    if path_apero[-1] != '/':
        path_apero = path_apero + '/'

    # List of APERO rdb files
    rdb_list = glob.glob(path_apero + "lbl_*.rdb")

    # Sort the list of APERO rdb file
    rdb_list.sort()

    # Iterate over files
    for i_rbd, rdb_i in enumerate(rdb_list):
        # Get target name and template name
        target, template = get_target_template(rdb_i)

        # Put APERO data into a dictionary ----------------------------------------------------------------------------
        # Get date, rv, rv_err from rdb file
        rjd_apero, rv_apero, rv_err_apero = open_rdb(rdb_i)

        # Process rvs (remove outliers, get JD, compute median and std, etc.)
        dict_apero = process_rvs(rjd_apero, rv_apero, rv_err_apero, pipeline="APERO", rv_method="LBL", nsig=nsig)

        # Put ESO data into a dictionary ------------------------------------------------------------------------------
        # Get date, rv, rv_err from DACE
        rjd_eso, rv_eso, rv_err_eso = query_dace(target, mode=mode)

        # Process rvs (remove outliers, get JD, compute median and std, etc.)
        dict_eso = process_rvs(rjd_eso, rv_eso, rv_err_eso, pipeline="ESO", rv_method="LBL", nsig=nsig)

        # Compare RVs -------------------------------------------------------------------------------------------------
        compare_rvs(dict_1=dict_apero, dict_2=dict_eso, target=target, template=template, nsig=nsig,
                    path_savefig=path_savefig)

        exit()  # TODO remove this to loop through all rdb files

    return 0


def compare_ccf(path_apero, compare_target=None, apero_ccf_dict_fname="apero_ccf_dict.h5", nsig=10, path_savefig=''):
    """
    OBSOLETE
    Compare APERO to ESO CCF radial velocities.

    :param path_apero:              (str) Absolute path to APERO reduced data folder.
    :param compare_target:          (str or list of str, optional) Name of a/multiple specific target(s) to compare.
                                    Default: None (iterate through all targets).
    :param apero_ccf_dict_fname:    (str, optional) Name of the dictionary to save or to read targets, nights, and CCF
                                    file names.
                                    Default: "apero_ccf_dict.h5".
    :param nsig:                    (float, optional) Number of sigmas away from median velocity to reject outliers.
                                    Default: 10.
    :param path_savefig:            (str, optional) Path where to save figures. Default: ''.
    :return: 0
    """

    # APERO -----------------------------------------------------------------------------------------------------------
    if path_apero[-1] != '/':
        path_apero = path_apero + '/'

    # List of APERO nights
    night_list = glob.glob(path_apero + "2*")

    # Sort the list of APERO nights
    night_list.sort()

    if os.path.isfile(path_savefig + apero_ccf_dict_fname):
        targets = loadhdf5(path_savefig + apero_ccf_dict_fname)
    else:
        # Make a dictionary of targets, their respective list of nights, and their respective CCF file names
        targets = {}
        # Iterate over nights
        for i_night, night_i in enumerate(night_list):
            night_short = night_i.split('/')[-1]  # YYYY-MM-DD

            # Get all CCF files in this night
            ccf_files = glob.glob(night_i + "/NIRPS*_pp_e2dsff_tcorr_A_ccf_*_neg.fits_A.fits")

            # Iterate over CCF files to get target
            for ccf_file in ccf_files:
                target = fits.getheader(ccf_file)["DRSOBJN"]
                if target not in targets.keys():
                    targets[target] = {night_short: [ccf_file]}
                else:
                    if night_short not in targets[target].keys():
                        targets[target][night_short] = [ccf_file]
                    else:
                        targets[target][night_short].append(ccf_file)
        # Save dictionary to avoid re-building at every run
        # _ = savehdf5(targets, path_savefig + apero_ccf_dict_fname)

    # Read CCF RVs
    rjd_apero, rv_apero, rv_err_apero = {}, {}, {}
    # Iterate over targets
    if compare_target is None:
        compare_target = list(targets.keys())
    else:
        if type(compare_target) == str:
            compare_target = [compare_target]
    for i_target, target_i in enumerate(compare_target):
        if target_i not in targets.keys():
            print("{} not found. We skip it.".format(target_i))
            continue
        # Iterate over nights
        for j_night, night_j in enumerate(targets[target_i].keys()):
            # Iterate over CCF files
            for k_ccffile, ccffile_k in enumerate(targets[target_i][night_j]):
                if j_night == 0:
                    rjd_apero[target_i] = []
                    rv_apero[target_i] = []
                    rv_err_apero[target_i] = []
                # Get rjd, rv, rv_err
                hdr = fits.getheader(ccffile_k)
                rjd_apero[target_i].append(hdr["MJDMID"] + .5)  # [MJD -> RJD]
                rv_apero[target_i].append(hdr["RV_CORR"] * 1e3)  # [km/s -> m/s]
                rv_err_apero[target_i].append(hdr["DVRMS_CC"])  # [m/s]
        rjd_apero[target_i] = np.array(rjd_apero[target_i])
        rv_apero[target_i] = np.array(rv_apero[target_i])
        rv_err_apero[target_i] = np.array(rv_err_apero[target_i])

    # Put ESO ---------------------------------------------------------------------------------------------------------
    # Iterate over targets
    for i_target, target_i in enumerate(rjd_apero.keys()):
        # APERO
        dict_apero_i = process_rvs(rjd_apero[target_i], rv_apero[target_i], rv_err_apero[target_i], pipeline="APERO",
                                   rv_method="CCF", nsig=nsig)

        # ESO
        # Get date, rv, rv_err from DACE
        rjd_eso, rv_eso, rv_err_eso = query_dace(target_i, mode=mode)
        # Process rvs (remove outliers, get JD, compute median and std, etc.)
        dict_eso_i = process_rvs(rjd_eso, rv_eso, rv_err_eso, pipeline="ESO", rv_method="CCF", nsig=nsig)

        compare_rvs(dict_1=dict_apero_i, dict_2=dict_eso_i, target=target_i, template=None, nsig=10,
                    path_savefig=path_savefig)
    return 0


def compare_rvs(dict_1, dict_2, target, template=None, nsig=10, path_savefig=''):
    """
    Compare radial velocities from two reductions.

    :param dict_1:          (dict) Dictionary of APERO RVs.
    :param dict_2:          (dict) Dictionary of ESO RVs.
    :param target:          (str) Target name.
    :param template:        (str, optional) Template name. Default: None.
    :param nsig:            (float, optional) Number of sigmas away from median velocity to reject outliers.
                            Default: 10.
    :param path_savefig:    (str, optional) Path where to save figures. Default: ''.
    :return: 0
    """
    # Convert target and template names to upper case
    target = target.upper()
    if template is not None:
        template = template.upper()
    print("TARGET: {}".format(target))

    # Compare ---------------------------------------------------------------------------------------------------------
    # Search for missing data in both pipelines - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Find all observation dates (JDs rounded to integers)
    lasilla_midnight = .5 + 4 / 24  # midnight in La Silla (roughly, depends on hour change), in JD
    nights_1 = np.unique(np.round(dict_1["full"]["jd"] - lasilla_midnight))
    nights_2 = np.unique(np.round(dict_2["full"]["jd"] - lasilla_midnight))
    # Convert rounded JDs to JDs at midnight in La Silla
    nights_1 = nights_1 + lasilla_midnight  # midnight in La Silla (roughly)
    nights_2 = nights_2 + lasilla_midnight  # midnight in La Silla (roughly)
    # Convert array of JDs to nightly bins
    nights_1_bins = [[n_i - .5, n_i + .5] for n_i in nights_1]
    nights_2_bins = [[n_i - .5, n_i + .5] for n_i in nights_2]

    # Find missing data points
    discard_in_1 = np.array([], dtype=int)
    discard_in_2 = np.array([], dtype=int)
    # Iterate over APERO nights
    for i_night_bin, night_bin_i in enumerate(nights_1_bins):

        # Get night start and end times
        night_start, night_end = night_bin_i[0], night_bin_i[1]

        # Get indices of points in this night
        i_pts_1 = (night_start < dict_1["full"]["jd"]) * (dict_1["full"]["jd"] <= night_end)
        i_pts_2 = (night_start < dict_2["full"]["jd"]) * (dict_2["full"]["jd"] <= night_end)

        # Get number of points in this night
        npts_1 = np.sum(i_pts_1)
        npts_2 = np.sum(i_pts_2)

        if npts_1 != npts_2:
            night_start_ymd = Time(night_start, format="jd").to_datetime()
            night_end_ymd = Time(night_end, format="jd").to_datetime()
            print("Night {} - {}: {} pts in {} {}, {} pts in {} {}".format(night_start_ymd, night_end_ymd,
                                                                           npts_1, dict_1["pipeline"],
                                                                           dict_1["rv_method"],
                                                                           npts_2, dict_2["pipeline"],
                                                                           dict_2["rv_method"]))
            # Discard the entire night for both reductions (cannot compare)
            if npts_1 > 0:
                discard_in_1 = np.concatenate((discard_in_1, np.where(i_pts_1)[0]))
            if npts_2 > 0:
                discard_in_2 = np.concatenate((discard_in_2, np.where(i_pts_2)[0]))
    # Iterate over ESO nights
    for i_night_bin, night_bin_i in enumerate(nights_2_bins):

        # Get night start and end times
        night_start, night_end = night_bin_i[0], night_bin_i[1]

        # Get indices of points in this night
        i_pts_1 = (night_start < dict_1["full"]["jd"]) * (dict_1["full"]["jd"] <= night_end)
        i_pts_2 = (night_start < dict_2["full"]["jd"]) * (dict_2["full"]["jd"] <= night_end)

        # Get number of points in this night
        npts_1 = np.sum(i_pts_1)
        npts_2 = np.sum(i_pts_2)

        if npts_2 != npts_1:
            night_start_ymd = Time(night_start, format="jd").to_datetime()
            night_end_ymd = Time(night_end, format="jd").to_datetime()
            print("Night {} - {}: {} pts in {} {}, {} pts in {} {}".format(night_start_ymd, night_end_ymd,
                                                                           npts_2, dict_2["pipeline"],
                                                                           dict_2["rv_method"],
                                                                           npts_1, dict_1["pipeline"],
                                                                           dict_1["rv_method"]))
            # Discard the entire night for both reductions (cannot compare)
            if npts_2 > 0:
                discard_in_2 = np.concatenate((discard_in_2, np.where(i_pts_2)[0]))
            if npts_1 > 0:
                discard_in_1 = np.concatenate((discard_in_1, np.where(i_pts_1)[0]))
    # Keep only unique indices
    discard_in_1 = np.unique(discard_in_1)
    discard_in_2 = np.unique(discard_in_2)

    # For comparisons, keep only dates where we have both APERO and ESO data
    igood_1 = np.arange(dict_1["full"]["jd"].size)
    igood_1 = np.delete(igood_1, discard_in_1)
    igood_2 = np.arange(dict_2["full"]["jd"].size)
    igood_2 = np.delete(igood_2, discard_in_2)

    diff, jd_diff, std_diff, med_diff, i_diff_outliers_up, diff_outliers_up, jd_diff_outliers_up, i_diff_outliers_lo, \
        diff_outliers_lo, jd_diff_outliers_lo, diff_no_outliers, jd_diff_no_outliers = None, None, None, None, None, \
        None, None, None, None, None, None, None
    if igood_1.size != 0 and igood_2.size != 0:
        # Compute difference
        diff = dict_1["full"]["rv"][igood_1] - dict_2["full"]["rv"][igood_2]

        # Make a JD array for the differences
        jd_diff = dict_1["full"]["jd"][igood_1]

        # Compute difference standard deviation in differences
        std_diff = np.nanpercentile(diff, q=[16, 84])
        std_diff = (std_diff[-1] - std_diff[0]) / 2

        # Compute median diff
        med_diff = np.nanmedian(diff)

        # Remove outliers (but store them to plot them)
        i_diff_outliers = np.where(np.abs(diff - med_diff) / std_diff > nsig)
        i_diff_outliers_up = (diff - med_diff) / std_diff > nsig
        i_diff_outliers_lo = (diff - med_diff) / std_diff < -nsig
        diff_no_outliers = np.delete(diff, i_diff_outliers)
        jd_diff_no_outliers = np.delete(jd_diff, i_diff_outliers)
        diff_outliers_up, diff_outliers_lo, jd_diff_outliers_up, jd_diff_outliers_lo = None, None, None, None
        if np.sum(i_diff_outliers_up) > 0:
            diff_outliers_up = np.ones(np.sum(i_diff_outliers_up)) * np.nanmax(diff_no_outliers)
            jd_diff_outliers_up = jd_diff[i_diff_outliers_up]
        if np.sum(i_diff_outliers_lo) > 0:
            diff_outliers_lo = np.ones(np.sum(i_diff_outliers_lo)) * np.nanmin(diff_no_outliers)
            jd_diff_outliers_lo = jd_diff[i_diff_outliers_lo]

    # Plot --------------------------------------------------------------------------------------------------------
    # Set up figure/axes
    fig, axs = plt.subplots(nrows=2, ncols=2, sharex="col", sharey="row", figsize=(20, 12),
                            gridspec_kw={"width_ratios": [2, 1]})

    axs[0, 0].set(ylabel="Relative velocity [m/s]")
    axs[1, 0].set(ylabel="Velocity difference [m/s]")
    axs[-1, 0].set(xlabel="Date")

    # Plot APERO data
    # Plot data, median-subtracted
    axs[0, 0].errorbar(dict_1["no_outliers"]["jd"], dict_1["no_outliers"]["rv"]
                       - dict_1["no_outliers"]["median"], dict_1["no_outliers"]["rv_err"],
                       ls='', marker='.', alpha=.5, color="C0", capsize=capsize,
                       label=("{} {} ({} pts)\nmed. = {:.2f} m/s\nmed.err. = {:.2f} m/s\nstd = {:.2f} m/s\nptp scat. = {:.2f} m/s"
                              .format(dict_1["pipeline"], dict_1["rv_method"], dict_1["full"]["jd"].size,
                                      dict_1["no_outliers"]["median"],
                                      dict_1["no_outliers"]["mederr"],
                                      dict_1["no_outliers"]["std"],
                                      dict_1["no_outliers"]["ptpscat"]
                                      )
                              )
                       )
    axs[0, 0].axhline(0, color="k", ls=':')

    # Plot ESO data
    # Plot data, median-subtracted
    axs[0, 0].errorbar(dict_2["no_outliers"]["jd"], dict_2["no_outliers"]["rv"]
                       - dict_2["no_outliers"]["median"], dict_2["no_outliers"]["rv_err"],
                       ls='', marker='.', alpha=.5, color="C1", capsize=capsize,
                       label=("\n{} {} ({} pts)\nmed. = {:.2f} m/s\nmed.err. = {:.2f} m/s\nstd = {:.2f} m/s\nptp scat. = {:.2f} m/s"
                              .format(dict_2["pipeline"], dict_2["rv_method"], dict_2["full"]["jd"].size,
                                      dict_2["no_outliers"]["median"],
                                      dict_2["no_outliers"]["mederr"],
                                      dict_2["no_outliers"]["std"],
                                      dict_2["no_outliers"]["ptpscat"]
                                      )
                              )
                       )

    # Plot outliers
    ax0_ymin, ax0_ymax = axs[0, 0].get_ylim()
    # In APERO
    if dict_1["outliers_up"]["jd"].size > 0:
        axs[0, 0].plot(dict_1["outliers_up"]["jd"], ax0_ymax * np.ones_like(dict_1["outliers_up"]["jd"]),
                       ls='', marker=r'$\uparrow$', alpha=.2, color="C0", markersize=markersize)
    if dict_1["outliers_lo"]["jd"].size > 0:
        axs[0, 0].plot(dict_1["outliers_lo"]["jd"], ax0_ymin * np.ones_like(dict_1["outliers_lo"]["jd"]),
                       ls='', marker=r'$\downarrow$', alpha=.2, color="C0", markersize=markersize)
    # In ESO
    if dict_2["outliers_up"]["jd"].size > 0:
        axs[0, 0].plot(dict_2["outliers_up"]["jd"], ax0_ymax * np.ones_like(dict_2["outliers_up"]["jd"]),
                       ls='', marker=r'$\uparrow$', alpha=.2, color="C1", markersize=markersize)
    if dict_2["outliers_lo"]["jd"].size > 0:
        axs[0, 0].plot(dict_2["outliers_lo"]["jd"], ax0_ymin * np.ones_like(dict_2["outliers_lo"]["jd"]),
                       ls='', marker=r'$\downarrow$', alpha=.2, color="C1", markersize=markersize)

    if igood_1.size != 0 and igood_2.size != 0:
        # Plot difference between APERO and ESO
        axs[1, 0].plot(jd_diff_no_outliers, diff_no_outliers,
                       ls='', marker='.', color='k',
                       label="{} {} - {} {}, std = {:.2f} m/s".format(dict_1["pipeline"], dict_1["rv_method"],
                                                                      dict_2["pipeline"], dict_2["rv_method"],
                                                                      std_diff))
        axs[1, 0].axhline(med_diff, color='k', ls=':')
        # Plot outliers
        if np.sum(i_diff_outliers_up) > 0:
            axs[1, 0].plot(jd_diff_outliers_up, diff_outliers_up, ls='', marker=r"$\uparrow$", alpha=.2, color='k',
                           markersize=markersize, label="Outliers (up/down)")
        if np.sum(i_diff_outliers_lo) > 0:
            axs[1, 0].plot(jd_diff_outliers_lo, diff_outliers_lo, ls='', marker=r"$\downarrow$", alpha=.2, color='k',
                           markersize=markersize)
        # Plot a line at missing or discarded data
        for discard in discard_in_1:
            axs[1, 0].axvline(dict_1["full"]["jd"][discard], color="C0", lw=.5,
                              label="{} {} pt".format(dict_1["pipeline"], dict_1["rv_method"])
                              if discard == discard_in_1[0] else '')
        for discard in discard_in_2:
            axs[1, 0].axvline(dict_2["full"]["jd"][discard], color="C1", lw=.5,
                              label="{} {} pt".format(dict_2["pipeline"], dict_2["rv_method"])
                              if discard == discard_in_2[0] else '')

        # Plot histogram of differences between APERO and ESO
        axs[1, 1].set(xlabel="Count")
        axs[1, 1].hist(diff_no_outliers, bins=9, histtype="step", color='k', orientation="horizontal")
        # Plot single point to show error bar size
        axs[1, 1].errorbar(1, med_diff, dict_1["no_outliers"]["mederr"],
                           ls='', marker='', color="C0", capsize=capsize,
                           label="APERO {} med. err. bar".format(dict_1["rv_method"]))
        axs[1, 1].errorbar(2, med_diff, dict_2["no_outliers"]["mederr"],
                           ls='', marker='', color="C1", capsize=capsize,
                           label="ESO {} med. err. bar".format(dict_2["rv_method"]))

    # Final edits to the plot -------------------------------------------------------------------------------------
    # Put legend outside subplots for row 1
    h, l = axs[0, 0].get_legend_handles_labels()
    axs[0, 1].legend(h, l, borderaxespad=0, fontsize=fontsize - 4, loc=2)
    axs[0, 1].axis("off")
    axs[1, 0].legend(frameon=False, fontsize=fontsize - 4, ncols=3)
    axs[1, 1].legend(frameon=False, fontsize=fontsize - 4)
    axs[0, 0].set(title=("Target: {} | {} {}{} | {} {}{}"
                         .format(target,
                                 dict_1["pipeline"],
                                 dict_1["rv_method"],
                                 " (template: {})".format(template) if dict_1["rv_method"] == "LBL" else '',
                                 dict_2["pipeline"],
                                 dict_2["rv_method"],
                                 " (template: {})".format(template) if dict_2["rv_method"] == "LBL" else '')
                         )
                  )
    fig.tight_layout()

    # Make a YYYY-MM-DD x axis
    xlims = axs[-1, 0].get_xlim()
    xlims = (int(xlims[0]), int(xlims[1]))
    deltax = np.max((int((xlims[1] - xlims[0]) // 4), 1))
    if deltax > 1:
        xticks = np.arange(xlims[0], xlims[1] + deltax / 2, deltax)
        xlabels = Time(xticks, format="jd").to_datetime()
        xlabels = ["{:4d}-{:02d}-{:02d}".format(xlab_i.year, xlab_i.month, xlab_i.day) for xlab_i in xlabels]
    else:
        xticks = axs[-1, 0].get_xticks()
        xticks = xticks[::xticks.size // 4]
        xlabels = Time(xticks, format="jd").to_datetime()
        xlabels = ["{:4d}-{:02d}-{:02d} {:02d}:{:02d}"
                   .format(xlab_i.year, xlab_i.month, xlab_i.day, xlab_i.hour, xlab_i.minute) for xlab_i in xlabels]
    axs[-1, 0].set_xticks(ticks=xticks, labels=xlabels)

    # Save plot ---------------------------------------------------------------------------------------------------
    # Make a file name
    if path_savefig != '':
        if path_savefig[-1] != '/':
            path_savefig = path_savefig + '/'
    savefname = path_savefig + "{}_{}-{}{}_{}-{}{}".format(target,
                                                           dict_1["pipeline"], dict_1["rv_method"],
                                                           "-{}".format(template) if dict_1["rv_method"] == "LBL"
                                                           else '',
                                                           dict_2["pipeline"], dict_2["rv_method"],
                                                           "-{}".format(template) if dict_2["rv_method"] == "LBL"
                                                           else ''
                                                           )
    # Save
    # fig.savefig(path_savefig + savefname + ".png")
    # fig.savefig(path_savefig + savefname + ".pdf")
    # print("Saved in {}.png/pdf".format(savefname))

    plt.show()

    return 0


def run_comparison(target, template=None,
                   path_1='', path_2='', pipeline_1="APERO", pipeline_2="ESO", rv_method_1="LBL", rv_method_2="CCF",
                   mode="he", version="07276", onoff="online", path_savefig='', nsig=10):
    """
    Compare any set of 2 RV time series.

    :param target:          (str) Target name.
    :param template:        (str, optional) Template name if one of the two or both RV methods is/are LBL.
                            Default: None.
    :param path_1:          (str, optional) Absolute path to the data for time series 1. Default: ''.
    :param path_2:          (str, optional) Absolute path to the data for time series 2. Default: ''.
    :param pipeline_1:      (str, optional) Pipeline name. Can be either "APERO" or "ESO". Default: "APERO".
    :param pipeline_2:      (str, optional) Pipeline name. Can be either "APERO" or "ESO". Default: "ESO".
    :param rv_method_1:     (str, optional) Name of RV method for time series 1. Can be either "LBL" or "CCF".
                            Default: "LBL".
    :param rv_method_2:     (str, optional) Name of RV method for time series 2. Can be either "LBL" or "CCF".
                            Default: "CCF".
    :param mode:            (str, optional) Instrument mode. Can be either "he" (high efficiency) or "ha" (high
                            accuracy). Default: "he".
    :param version:         (str, optional) Reduction version. Default: "07276".
    :param onoff:           (str, optional) Reduction mode. Can be either "online" or "offline". Default: "online".
    :param path_savefig:    (str, optional) Absolute path where figures will be saved. Default: ''.
    :param nsig:            (float, optional) Number of sigmas to reject outliers. Default: 10.
    :return:                0
    """
    # Sanity checks
    pipeline_1, pipeline_2 = pipeline_1.upper(), pipeline_2.upper()
    assert mode in ["he", "ha"], "ERROR! 'mode' must be either 'he' (high efficiency) or 'ha' (high accuracy)."
    assert version in ["07276"], "ERROR! 'version' must be '07276'."  # add options when more than 1 are availabe
    assert onoff in ["online", "offline"], "ERROR! 'onoff' must be either 'online' or 'offline'."
    assert pipeline_1 in ["APERO", "ESO"], "ERROR! 'pipeline_1' must be either 'APERO' or 'ESO'."
    assert pipeline_2 in ["APERO", "ESO"], "ERROR! 'pipeline_2' must be either 'APERO' or 'ESO'."
    assert rv_method_1 in ["CCF", "LBL"], "ERROR! 'rv_method_1' must be either 'CCF' or 'LBL'."
    assert rv_method_2 in ["CCF", "LBL"], "ERROR! 'rv_method_2' must be either 'CCF' or 'LBL'."
    if rv_method_1 == "LBL" or rv_method_2 == "LBL":
        assert template is not None, "ERROR! 'template' must be provided to load LBL data."

    # Put data into a dictionaries ------------------------------------------------------------------------------------
    paths = [path_1, path_2]
    pipelines = [pipeline_1, pipeline_2]
    rv_methods = [rv_method_1, rv_method_2]
    dicts = []
    for i, (path, pipeline, rv_method) in enumerate(zip(paths, pipelines, rv_methods)):
        if rv_method == "LBL":
            if pipeline == "APERO":
                # Get date, rv, rv_err from rdb file
                rjd, rv, rv_err = open_rdb(path + "lbl_{}_{}.rdb".format(target.upper(), template.upper()))
            else:  # ESO: No LBL at the moment
                print("ERROR! LBL RVs have not been computed with ESO reduction yet.")
                exit()
        else:  # CCF
            if pipeline == "APERO":
                # Lists to store times, RVs, RV_ERRs
                rjd, rv, rv_err = [], [], []

                # List of APERO nights
                night_list = glob.glob(path + "2*")

                # Sort the list of APERO nights
                night_list.sort()

                # Iterate over nights
                for i_night, night_i in enumerate(night_list):
                    # Get all CCF files in this night
                    ccf_files = glob.glob(night_i + "/NIRPS*_pp_e2dsff_tcorr_A_ccf_*_neg.fits_A.fits")

                    # Iterate over CCF files to get target
                    for ccf_file in ccf_files:
                        if fits.getheader(ccf_file)["DRSOBJN"] != target.upper():
                            continue
                        # Read RJD, RV, RV_ERR
                        hdr = fits.getheader(ccf_file)
                        rjd.append(hdr["MJDMID"] + .5)  # [MJD -> RJD]
                        rv.append(hdr["RV_CORR"] * 1e3)  # [km/s -> m/s]
                        rv_err.append(hdr["DVRMS_CC"])  # [m/s]

                # Convert lists to arrays
                rjd, rv, rv_err = np.array(rjd), np.array(rv), np.array(rv_err)

            else:  # ESO
                # Get date, rv, rv_err from DACE
                rjd, rv, rv_err = query_dace(target, mode=mode)

        # Process rvs (remove outliers, get JD, compute median and std, etc.)
        dicts.append(process_rvs(rjd=rjd, rv=rv, rv_err=rv_err, pipeline=pipeline, rv_method=rv_method, nsig=nsig))

    # Compare RVs -----------------------------------------------------------------------------------------------------
    compare_rvs(dict_1=dicts[0], dict_2=dicts[1], target=target, template=template, nsig=nsig,
                path_savefig=path_savefig)

    return 0


# =====================================================================================================================
# Start of code
# =====================================================================================================================
if __name__ == '__main__':
    run_comparison(target=target, template=template,
                   path_1=path_1, path_2=path_2, pipeline_1=pipeline_1, pipeline_2=pipeline_2,
                   rv_method_1=rv_method_1, rv_method_2=rv_method_2,
                   mode=mode, version=version, onoff=onoff, path_savefig=path_savefig, nsig=nsig)


# =====================================================================================================================
# End of code
# =====================================================================================================================
