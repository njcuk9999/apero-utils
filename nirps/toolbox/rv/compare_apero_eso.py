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
from dace_query.spectroscopy import Spectroscopy

# For uniform fonts in plots
fontsize = 16
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
mode = "he"  # "ha": high accuracy | "he": high efficiency
version = "07276"  # version number I assume?
onoff = "online"  # "online" | "offline"

# Path to lbl rdb files
path_apero_lbl = "/cosmos99/nirps/lbl-data/nirps_{}_{}_{}/lblrdb/".format(mode, version, onoff)  # APERO

# Path to ccf files
path_apero_ccf = "/cosmos99/nirps/apero-data/nirps_{}_{}_{}/red/".format(mode, version, onoff)  # APERO

# Target
target = "Proxima"
template = "Proxima"

# RV method
rv_method_apero = "LBL"
rv_method_eso = "CCF"

# Path to save figures (and ccf dictionary)
path_savefig = ""  # TODO Where?
# apero_ccf_dict_fname = "apero_ccf_dict.h5"

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


def query_dace(target):
    """
    Query DACE to get RV time series
    :param target:  (str) Target name
    :return: rdj    (1d array), rv (1d array), rv_err (1d array)
    """
    # Retrieve radial velocity timeseries from DACE
    dict_eso_dace = Spectroscopy.get_timeseries(target=target, sorted_by_instrument=True, output_format='numpy')
    dict_eso_dace = dict_eso_dace["NIRPS"]["3.0.0"][mode.upper()]

    # Get date, rv, rv_err
    rjd, rv, rv_err = dict_eso_dace["rjd"], dict_eso_dace["rv"], dict_eso_dace["rv_err"]

    return rjd, rv, rv_err


def process_rvs(rjd, rv, rv_err, rv_method, nsig=10):
    """
    Compute median, standard deviation, and median error bar, and remove outliers.

    :param rjd:         (1d array) Reduced Julian day array.
    :param rv:          (1d array) Radial velocity array.
    :param rv_err:      (1d array) Radial velocity uncertainty array.
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

    # Make dictionary to store all arrays
    rv_dict = {"full": {"jd": jd_full,
                        "rv": rv_full,
                        "rv_err": rv_err_full,
                        "std": std_full,
                        "median": median_full,
                        "mederr": mederr_full},
               "no_outliers": {"jd": jd,
                               "rv": rv,
                               "rv_err": rv_err,
                               "std": std,
                               "median": median,
                               "mederr": mederr,
                               "indices": inonout},
               "outliers_up": {"jd": jd_out_up,
                               # "rv": rv_out_up,  # obsolete
                               "indices": iout_up},
               "outliers_lo": {"jd": jd_out_lo,
                               # "rv": rv_out_lo,  # obsolete
                               "indices": iout_lo},
               "rv_method": rv_method
               }

    return rv_dict


def compare_lbl(path_apero, nsig=10):
    """
    Compare APERO to ESO LBL radial velocities.

    :param path_apero:  (str) Absolute path to APERO LBL rdb folder.
    :param nsig:        (float, optional) Number of sigmas away from median velocity to reject outliers. Default: 10.
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
        dict_apero = process_rvs(rjd_apero, rv_apero, rv_err_apero, "LBL", nsig=nsig)

        # Put ESO data into a dictionary ------------------------------------------------------------------------------
        # Get date, rv, rv_err from DACE
        rjd_eso, rv_eso, rv_err_eso = query_dace(target)

        # Process rvs (remove outliers, get JD, compute median and std, etc.)
        dict_eso = process_rvs(rjd_eso, rv_eso, rv_err_eso, "LBL", nsig=nsig)

        # Compare RVs -------------------------------------------------------------------------------------------------
        compare_rvs(dict_apero=dict_apero, dict_eso=dict_eso, target=target, template=template, nsig=nsig)

        exit()  # TODO remove this to loop through all rdb files

    return 0


def compare_ccf(path_apero, compare_target=None, apero_ccf_dict_fname="apero_ccf_dict.h5", nsig=10):
    """
    Compare APERO to ESO CCF radial velocities.

    :param path_apero:              (str) Absolute path to APERO reduced data folder.
    :param compare_target:          (str or list of str, optional) Name of a/multiple specific target(s) to compare.
                                    Default: None (iterate through all targets).
    :param apero_ccf_dict_fname:    (str, optional) Name of the dictionary to save or to read targets, nights, and CCF
                                    file names.
                                    Default: "apero_ccf_dict.h5".
    :param nsig:                    (float, optional) Number of sigmas away from median velocity to reject outliers.
                                    Default: 10.
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
                target = fits.getheader(ccf_file)["OBJECT"]
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
        dict_apero_i = process_rvs(rjd_apero[target_i], rv_apero[target_i], rv_err_apero[target_i], rv_method="CCF",
                                   nsig=nsig)

        # ESO
        # Get date, rv, rv_err from DACE
        rjd_eso, rv_eso, rv_err_eso = query_dace(target_i)
        # Process rvs (remove outliers, get JD, compute median and std, etc.)
        dict_eso_i = process_rvs(rjd_eso, rv_eso, rv_err_eso, rv_method="CCF", nsig=nsig)

        compare_rvs(dict_apero=dict_apero_i, dict_eso=dict_eso_i, target=target_i, template=None,
                    nsig=10)
    return 0


def compare_rvs(dict_apero, dict_eso, target, template=None, nsig=10, path_savefig=''):
    """
    Compare radial velocities from two reductions.

    :param dict_apero:      (dict) Dictionary of APERO RVs.
    :param dict_eso:        (dict) Dictionary of ESO RVs.
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

    # Compare APERO and ESO ---------------------------------------------------------------------------------------
    # Search for missing data in both pipelines - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Find all observation dates (JDs rounded to integers)
    lasilla_midnight = .5 + 4 / 24  # midnight in La Silla (roughly, depends on hour change), in JD
    nights_apero = np.unique(np.round(dict_apero["full"]["jd"] - lasilla_midnight))
    nights_eso = np.unique(np.round(dict_eso["full"]["jd"] - lasilla_midnight))
    # Convert rounded JDs to JDs at midnight in La Silla
    nights_apero = nights_apero + lasilla_midnight  # midnight in La Silla (roughly)
    nights_eso = nights_eso + lasilla_midnight  # midnight in La Silla (roughly)
    # Convert array of JDs to nightly bins
    nights_apero_bins = [[n_i - .5, n_i + .5] for n_i in nights_apero]
    nights_eso_bins = [[n_i - .5, n_i + .5] for n_i in nights_eso]

    # Find missing data points
    discard_in_apero = np.array([], dtype=int)
    discard_in_eso = np.array([], dtype=int)
    # Iterate over APERO nights
    for i_night_bin, night_bin_i in enumerate(nights_apero_bins):

        # Get night start and end times
        night_start, night_end = night_bin_i[0], night_bin_i[1]

        # Get indices of points in this night
        i_pts_apero = (night_start < dict_apero["full"]["jd"]) * (dict_apero["full"]["jd"] <= night_end)
        i_pts_eso = (night_start < dict_eso["full"]["jd"]) * (dict_eso["full"]["jd"] <= night_end)

        # Get number of points in this night
        npts_apero = np.sum(i_pts_apero)
        npts_eso = np.sum(i_pts_eso)

        if npts_apero != npts_eso:
            night_start_ymd = Time(night_start, format="jd").to_datetime()
            night_end_ymd = Time(night_end, format="jd").to_datetime()
            print("Night {} - {}: {} pts in APERO, {} pts in ESO".format(night_start_ymd, night_end_ymd,
                                                                         npts_apero, npts_eso))
            # Discard the entire night for both reductions (cannot compare)
            if npts_apero > 0:
                discard_in_apero = np.concatenate((discard_in_apero, np.where(i_pts_apero)[0]))
            if npts_eso > 0:
                discard_in_eso = np.concatenate((discard_in_eso, np.where(i_pts_eso)[0]))
    # Iterate over ESO nights
    for i_night_bin, night_bin_i in enumerate(nights_eso_bins):

        # Get night start and end times
        night_start, night_end = night_bin_i[0], night_bin_i[1]

        # Get indices of points in this night
        i_pts_eso = (night_start < dict_eso["full"]["jd"]) * (dict_eso["full"]["jd"] <= night_end)
        i_pts_apero = (night_start < dict_apero["full"]["jd"]) * (dict_apero["full"]["jd"] <= night_end)

        # Get number of points in this night
        npts_eso = np.sum(i_pts_eso)
        npts_apero = np.sum(i_pts_apero)

        if npts_eso != npts_apero:
            night_start_ymd = Time(night_start, format="jd").to_datetime()
            night_end_ymd = Time(night_end, format="jd").to_datetime()
            print("Night {} - {}: {} pts in ESO, {} pts in APERO".format(night_start_ymd, night_end_ymd,
                                                                         npts_eso, npts_apero))
            # Discard the entire night for both reductions (cannot compare)
            if npts_eso > 0:
                discard_in_eso = np.concatenate((discard_in_eso, np.where(i_pts_eso)[0]))
            if npts_apero > 0:
                discard_in_apero = np.concatenate((discard_in_apero, np.where(i_pts_apero)[0]))
    # Keep only unique indices
    discard_in_apero = np.unique(discard_in_apero)
    discard_in_eso = np.unique(discard_in_eso)

    # For comparisons, keep only dates where we have both APERO and ESO data
    igood_apero = np.arange(dict_apero["full"]["jd"].size)
    igood_apero = np.delete(igood_apero, discard_in_apero)
    igood_eso = np.arange(dict_eso["full"]["jd"].size)
    igood_eso = np.delete(igood_eso, discard_in_eso)

    diff, jd_diff, std_diff, med_diff, i_diff_outliers_up, diff_outliers_up, jd_diff_outliers_up, i_diff_outliers_lo, \
        diff_outliers_lo, jd_diff_outliers_lo, diff_no_outliers, jd_diff_no_outliers = None, None, None, None, None, \
        None, None, None, None, None, None, None
    if igood_apero.size != 0 and igood_eso.size != 0:
        # Compute difference
        diff = dict_apero["full"]["rv"][igood_apero] - dict_eso["full"]["rv"][igood_eso]

        # Make a JD array for the differences
        jd_diff = dict_apero["full"]["jd"][igood_apero]

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
    fig, axs = plt.subplots(nrows=3, ncols=2, sharex="col", sharey="row", figsize=(20, 12),
                            gridspec_kw={"width_ratios": [2, 1], "height_ratios":[1, 2, 2]})
    # axs[0, 1].remove()
    # axs[1, 1].remove()
    axs[0, 0].set(ylabel="Velocity [m/s]")
    axs[1, 0].set(ylabel="Relative velocity [m/s]")
    axs[2, 0].set(ylabel="Velocity difference [m/s]")
    axs[-1, 0].set(xlabel="Date")

    # Plot APERO data
    # Plot data
    axs[0, 0].errorbar(dict_apero["no_outliers"]["jd"], dict_apero["no_outliers"]["rv"],
                       dict_apero["no_outliers"]["rv_err"],
                       ls='', marker='.', alpha=.5, color="C0", capsize=capsize,
                       label="APERO ({})".format(dict_apero["rv_method"]))
    # Plot median
    axs[0, 0].axhline(dict_apero["no_outliers"]["median"], color="C0", ls=':')
    # Plot data, median-subtracted
    axs[1, 0].errorbar(dict_apero["no_outliers"]["jd"], dict_apero["no_outliers"]["rv"]
                       - dict_apero["no_outliers"]["median"], dict_apero["no_outliers"]["rv_err"],
                       ls='', marker='.', alpha=.5, color="C0", capsize=capsize,
                       label=("APERO ({}, {} pts)\nmed.err. = {:.2f} m/s\nstd = {:.2f} m/s"
                              .format(dict_apero["rv_method"], dict_apero["full"]["jd"].size,
                                      dict_apero["no_outliers"]["mederr"], dict_apero["no_outliers"]["std"])))
    axs[1, 0].axhline(0, color="k", ls=':')

    # Plot ESO data
    # Plot data
    axs[0, 0].errorbar(dict_eso["no_outliers"]["jd"], dict_eso["no_outliers"]["rv"],
                       dict_eso["no_outliers"]["rv_err"],
                       ls='', marker='.', alpha=.5, color="C1", capsize=capsize,
                       label="ESO ({})".format(dict_eso["rv_method"]))
    # Plot median
    axs[0, 0].axhline(dict_eso["no_outliers"]["median"], color="C1", ls=':')
    # Plot data, median-subtracted
    axs[1, 0].errorbar(dict_eso["no_outliers"]["jd"], dict_eso["no_outliers"]["rv"]
                       - dict_eso["no_outliers"]["median"], dict_eso["no_outliers"]["rv_err"],
                       ls='', marker='.', alpha=.5, color="C1", capsize=capsize,
                       label=("\nESO ({}, {} pts)\nmed.err. = {:.2f} m/s\nstd = {:.2f} m/s"
                              .format(dict_eso["rv_method"], dict_eso["full"]["jd"].size,
                                      dict_eso["no_outliers"]["mederr"], dict_eso["no_outliers"]["std"])))

    # Plot outliers
    ax0_ymin, ax0_ymax = axs[0, 0].get_ylim()
    ax1_ymin, ax1_ymax = axs[1, 0].get_ylim()
    # In APERO
    if dict_apero["outliers_up"]["jd"].size > 0:
        axs[0, 0].plot(dict_apero["outliers_up"]["jd"], ax0_ymax * np.ones_like(dict_apero["outliers_up"]["jd"]),
                       ls='', marker=r'$\uparrow$', alpha=.2, color="C0", markersize=markersize,
                       label="Outliers (any color, up/down)")
        axs[1, 0].plot(dict_apero["outliers_up"]["jd"], ax1_ymax * np.ones_like(dict_apero["outliers_up"]["jd"]),
                       ls='', marker=r'$\uparrow$', alpha=.2, color="C0", markersize=markersize)
    if dict_apero["outliers_lo"]["jd"].size > 0:
        axs[0, 0].plot(dict_apero["outliers_lo"]["jd"], ax0_ymin * np.ones_like(dict_apero["outliers_lo"]["jd"]),
                       ls='', marker=r'$\downarrow$', alpha=.2, color="C0", markersize=markersize)
        axs[1, 0].plot(dict_apero["outliers_lo"]["jd"], ax1_ymin * np.ones_like(dict_apero["outliers_lo"]["jd"]),
                       ls='', marker=r'$\downarrow$', alpha=.2, color="C0", markersize=markersize)
    # In ESO
    if dict_eso["outliers_up"]["jd"].size > 0:
        axs[0, 0].plot(dict_eso["outliers_up"]["jd"], ax0_ymax * np.ones_like(dict_eso["outliers_up"]["jd"]),
                       ls='', marker=r'$\uparrow$', alpha=.2, color="C1", markersize=markersize)
        axs[1, 0].plot(dict_eso["outliers_up"]["jd"], ax1_ymax * np.ones_like(dict_eso["outliers_up"]["jd"]),
                       ls='', marker=r'$\uparrow$', alpha=.2, color="C1", markersize=markersize)
    if dict_eso["outliers_lo"]["jd"].size > 0:
        axs[0, 0].plot(dict_eso["outliers_lo"]["jd"], ax0_ymin * np.ones_like(dict_eso["outliers_lo"]["jd"]),
                       ls='', marker=r'$\downarrow$', alpha=.2, color="C1", markersize=markersize)
        axs[1, 0].plot(dict_eso["outliers_lo"]["jd"], ax1_ymin * np.ones_like(dict_eso["outliers_lo"]["jd"]),
                       ls='', marker=r'$\downarrow$', alpha=.2, color="C1", markersize=markersize)

    if igood_apero.size != 0 and igood_eso.size != 0:
        # Plot difference between APERO and ESO
        axs[2, 0].plot(jd_diff_no_outliers, diff_no_outliers,
                       ls='', marker='.', color='k',
                       label="APERO ({}) - ESO ({}), std = {:.2f} m/s".format(dict_apero["rv_method"],
                                                                              dict_eso["rv_method"],
                                                                              std_diff))
        axs[2, 0].axhline(med_diff, color='k', ls=':')
        # Plot outliers
        if np.sum(i_diff_outliers_up) > 0:
            axs[2, 0].plot(jd_diff_outliers_up, diff_outliers_up, ls='', marker=r"$\uparrow$", alpha=.2, color='k',
                           markersize=markersize, label="Outliers (up/down)")
        if np.sum(i_diff_outliers_lo) > 0:
            axs[2, 0].plot(jd_diff_outliers_lo, diff_outliers_lo, ls='', marker=r"$\downarrow$", alpha=.2, color='k',
                           markersize=markersize)
        # Plot a line at missing or discarded data
        for discard in discard_in_apero:
            axs[2, 0].axvline(dict_apero["full"]["jd"][discard], color="C0", lw=.5,
                              label="APERO pt" if discard == discard_in_apero[0] else '')
        for discard in discard_in_eso:
            axs[2, 0].axvline(dict_eso["full"]["jd"][discard], color="C1", lw=.5,
                              label="ESO pt" if discard == discard_in_eso[0] else '')

        # Plot histogram of differences between APERO and ESO
        axs[2, 1].set(xlabel="Count")
        axs[2, 1].hist(diff_no_outliers, bins=9, histtype="step", color='k', orientation="horizontal")
        # Plot single point to show error bar size
        axs[2, 1].errorbar(1, med_diff, dict_apero["no_outliers"]["mederr"],
                           ls='', marker='', color="C0", capsize=capsize,
                           label="APERO ({}) med. err. bar".format(dict_apero["rv_method"]))
        axs[2, 1].errorbar(2, med_diff, dict_eso["no_outliers"]["mederr"],
                           ls='', marker='', color="C1", capsize=capsize,
                           label="ESO ({}) med. err. bar".format(dict_eso["rv_method"]))

    # Final edits to the plot -------------------------------------------------------------------------------------
    # Put legend outside subplots for rows 1 & 2
    h, l = axs[0, 0].get_legend_handles_labels()
    axs[0, 1].legend(h, l, borderaxespad=0, fontsize=fontsize - 4, loc=2)
    axs[0, 1].axis("off")
    h, l = axs[1, 0].get_legend_handles_labels()
    axs[1, 1].legend(h, l, borderaxespad=0, fontsize=fontsize - 4, loc=2)
    axs[1, 1].axis("off")
    # axs[0, 0].legend(frameon=False, fontsize=fontsize - 4)
    # axs[1, 0].legend(frameon=False, fontsize=fontsize - 4)
    axs[2, 0].legend(frameon=False, fontsize=fontsize - 4, ncols=3)
    axs[2, 1].legend(frameon=False, fontsize=fontsize - 4)
    axs[0, 0].set(title=("Target: {} | APERO: {}{} | ESO: {}{}"
                         .format(target,
                                 dict_apero["rv_method"],
                                 " (template: {})".format(template) if dict_apero["rv_method"] == "LBL" else '',
                                 dict_eso["rv_method"],
                                 " (template: {})".format(template) if dict_eso["rv_method"] == "LBL" else '')
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
    savefname = path_savefig + "{}_APERO-{}{}_ESO-{}{}".format(target,
                                                               dict_apero["rv_method"],
                                                               "-{}".format(template) if dict_apero["rv_method"] == "LBL"
                                                               else '',
                                                               dict_eso["rv_method"],
                                                               "-{}".format(template) if dict_eso["rv_method"] == "LBL"
                                                               else ''
                                                               )
    # Save
    # fig.savefig(path_savefig + savefname + ".png")
    # fig.savefig(path_savefig + savefname + ".pdf")
    # print("Saved in {}.png/pdf".format(savefname))

    plt.show()

    return 0


def run_comparison(path_apero, target, template=None, rv_method_apero="LBL", rv_method_eso="CCF", mode="he",
                   version="07276", onoff="online", path_savefig='', nsig=10):
    """"""  # TODO add docstring
    # Sanity checks
    assert mode in ["he", "ha"], "ERROR! 'mode' must be either 'he' (high efficiency) or 'ha' (high accuracy)."
    assert version in ["07276"], "ERROR! 'version' must be '07276'."  # add options when more than 1 are availabe
    assert onoff in ["online", "offline"], "ERROR! 'onoff' must be either 'online' or 'offline'."
    assert rv_method_apero in ["CCF", "LBL"], "ERROR! 'rv_method_apero' must be either 'CCF' or 'LBL'."
    assert rv_method_eso in ["CCF", "LBL"], "ERROR! 'rv_method_eso' must be either 'CCF' or 'LBL'."
    if rv_method_apero == "LBL" or rv_method_eso == "LBL":
        assert template is not None, "ERROR! 'template' must be provided to load LBL data."

    # Put APERO data into a dictionary --------------------------------------------------------------------------------
    if rv_method_apero == "LBL":
        # Get date, rv, rv_err from rdb file
        rjd_apero, rv_apero, rv_err_apero = open_rdb(path_apero + "lbl_{}_{}.rdb".format(target.upper(),
                                                                                         template.upper()))

    else:  # CCF
        # Lists to store times, RVs, RV_ERRs
        rjd_apero, rv_apero, rv_err_apero = [], [], []

        # List of APERO nights
        night_list = glob.glob(path_apero + "2*")

        # Sort the list of APERO nights
        night_list.sort()

        # Iterate over nights
        for i_night, night_i in enumerate(night_list):
            night_short = night_i.split('/')[-1]  # YYYY-MM-DD

            # Get all CCF files in this night
            ccf_files = glob.glob(night_i + "/NIRPS*_pp_e2dsff_tcorr_A_ccf_*_neg.fits_A.fits")

            # Iterate over CCF files to get target
            for ccf_file in ccf_files:
                if fits.getheader(ccf_file)["OBJECT"] != target:
                    continue
                # Read RJD, RV, RV_ERR
                hdr = fits.getheader(ccf_file)
                rjd_apero.append(hdr["MJDMID"] + .5)  # [MJD -> RJD]
                rv_apero.append(hdr["RV_CORR"] * 1e3)  # [km/s -> m/s]
                rv_err_apero.append(hdr["DVRMS_CC"])  # [m/s]

        # Convert lists to arrays
        rjd_apero, rv_apero, rv_err_apero = np.array(rjd_apero), np.array(rv_apero), np.array(rv_err_apero)

    # Process rvs (remove outliers, get JD, compute median and std, etc.)
    dict_apero = process_rvs(rjd_apero, rv_apero, rv_err_apero, rv_method_apero, nsig=nsig)

    # Put ESO data into a dictionary ----------------------------------------------------------------------------------
    if rv_method_eso == "LBL":
        '''
        # Get date, rv, rv_err from rdb file
        rjd_eso, rv_eso, rv_err_eso = open_rdb(path_eso_lbl + "lbl_{}_{}.rdb".format(target.upper(), template.upper()))

        # Process rvs (remove outliers, get JD, compute median and std, etc.)
        dict_eso = process_rvs(rjd_eso, rv_eso, rv_err_eso, "LBL", nsig=nsig)
        '''
        print("ERROR! LBL RVs have not been computed with ESO reduction yet.")
        exit()
    else:
        # Get date, rv, rv_err from DACE
        rjd_eso, rv_eso, rv_err_eso = query_dace(target)

    # Process rvs (remove outliers, get JD, compute median and std, etc.)
    dict_eso = process_rvs(rjd_eso, rv_eso, rv_err_eso, rv_method=rv_method_eso, nsig=nsig)

    # Compare RVs -----------------------------------------------------------------------------------------------------
    compare_rvs(dict_apero=dict_apero, dict_eso=dict_eso, target=target, template=template, nsig=nsig,
                path_savefig=path_savefig)


# =====================================================================================================================
# Start of code
# =====================================================================================================================
if __name__ == '__main__':
    run_comparison(path_apero=path_apero_lbl, target=target, template=template, rv_method_apero=rv_method_apero,
                   rv_method_eso=rv_method_eso, mode=mode, version=version, onoff=onoff, path_savefig=path_savefig,
                   nsig=nsig)


# =====================================================================================================================
# End of code
# =====================================================================================================================
