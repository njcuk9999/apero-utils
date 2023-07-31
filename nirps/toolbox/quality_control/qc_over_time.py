# Code for plotting various qc checks over time, for different types of apero reduced files

from astropy.io import fits
import numpy as np
from glob import glob
import matplotlib.pyplot as plt

# directory for the plots
fig_path = '/home/nirps-client/Frederic/qc_over_time/'

# dictionaries with useful information for relevant file types
# commented out for file types with nothing/no need to plot

# pp
# name and value are reversed for older files. Fixed in 0.7.284+
# 3 qc checks, all numerical
pp = {
    'file_extension': 'pp',
    'folder': ['tmp'],
    'fibers': False,
    'n_qcc': 3,
    'numerical': ['1','2','3'],
    'qcc_names': ['snr_hotpix', 'max(rms_list)', 'EXPTIME']
}
# e2ds
# only one qc check, not numerical, same as e2dsff
# e2ds = {
#     'file_extension': 'pp_e2ds_',
#     'folder': 'red',
#     'fibers': True
# }
# e2dsff
# only one qc check, not numerical
# e2dsff = {
#     'file_extension': 'pp_e2dsff_',
#     'folder': 'red',
#     'fibers': True
# }
# tcorr
# 8 qc checks, all numerical
# no files for fiber B
tcorr = {
    'file_extension': 'pp_e2dsff_tcorr_A',
    'folder': ['red'],
    'fibers': False,
    'n_qcc': 8,
    'numerical': ['1','2','3','4','5','6','7','8'],
    'qcc_names': ['MED[EXTSNR]', 'NUM_NAN_CCF', 'EXPO_OTHERS L', 'EXPO_OTHERS U', 'EXPO_WATER L', 'EXPO_WATER U',
                  'ITERATIONS', 'SNR[64]']
}
# recon
# same qc checks as tcorr, no need to plot twice
# no files for fiber B
# recon = {
#     'file_extension': 'pp_e2dsff_recon_A',
#     'folder': 'red',
#     'fibers': False
# }
# pclean
# no qc checks
# no files for fiber B
# pclean = {
#     'file_extension': 'pp_tellu_pclean_A',
#     'folder': 'red',
#     'fibers': False
# }
# wave
# 6 qc checks, all numerical
# put * to encompass all wave file types
# could implement different types separately, i.e. wave_night, wavesol, waveref, etc.
wave = {
    'file_extension': 'wave*',
    'folder': ['calib', 'red'],
    'fibers': True,
    'n_qcc': 6,
    'numerical': ['1','2','3','4','5','6'],
    'qcc_names': ['MIN YWAVE DIFF A', 'MIN XWAVE DIFF A', 'MIN YWAVE DIFF B', 'MIN XWAVE DIFF B', 'DV[A - B]',
                  'CCFRV[A] - CCFRV[B]']
}
# blaze
# 1 qc check, numerical
blaze = {
    'file_extension': 'pp_blaze_',
    'folder': ['calib','red'],
    'fibers': True,
    'n_qcc': 1,
    'numerical': ['1'],
    'qcc_names': ['max_rms']
}
# shapel
# 3 qcc, 2 numerical (2 and 3)
shapel = {
    'file_extension': 'pp_shapel',
    'folder': ['calib','red'],
    'fibers': False,
    'n_qcc': 3,
    'numerical': ['2','3'],
    'qcc_names': ['XRES', 'YRES']
}
# flat
# 1 qcc, numerical
flat = {
    'file_extension': 'pp_flat_',
    'folder': ['calib','red'],
    'fibers': True,
    'n_qcc': 1,
    'numerical': ['1'],
    'qcc_names': ['max_rms']
}
# loco
# 2 qcc, 1 numerical (1)
loco = {
    'file_extension': 'pp_loco_',
    'folder': ['calib','red'],
    'fibers': True,
    'n_qcc': 2,
    'numerical': ['1'],
    'qcc_names': ['rorder_num']
}
# other calib file types:
# bmap: no qcc
# badpixel: no qcc
# order_profile: 2 qcc, same as loco
# fwhm-order_A: 2 qcc, same as loco
# with-order_A: 2 qcc, same as loco
# dark: 1 file, no qcc

# other red file types?
# hclines, fplines included in wave
# shapel_orderps: 2 qcc, same values as loco

# files in the tellu folder?
# tellu_obj : 8 qcc, numerical or bool, but only values for 1, 2 and 3?
# not sure about these

# s1d files: all the same QCCs and values as e2ds counterparts

# waveref, waveref_res_e2ds, wavered_cav__
# all included in wave

# list of file types this code needs to be run for
# may be missing some for now
files_to_check = ['pp', 'tcorr', 'wave', 'blaze', 'shapel', 'flat', 'loco']


def pp_get_n_v(header):
    """
    Helper function to get the right qcc header keys for older pp files.
    Issue fixed in 0.7.284+

    :param header: the header of the fits file of interest
    :return: tuple (n, v), the proper letters N and V for the names and values of QCCs for the file.
    """
    # oldest nirps files have PVERSION = 0.7.282
    if header['PVERSION'] == '0.7.282' or header['PVERSION'] == '0.7.283':
        v = 'N'
        n = 'V'
    else:
        v = 'V'
        n = 'N'
    return n, v


def list_files(file_type: str, fiber='A'):
    """
    Function to browse the relevant directories and return all files of a given type. Browses calib, red and tmp in
    both HE and HA reduced data directories.

    :param file_type: str, the file type for which to list all files, as defined by the dictionaries and the list of
                      files_to_check
    :param fiber: 'A' or 'B'
    :return:
    """

    file_extension = eval(file_type)['file_extension']
    # check whether file_type has files for both A and B or not
    a_b = eval(file_type)['fibers']

    files = glob('/cosmos99/nirps/apero-data/nirps_h*online/calib/*' + file_extension + a_b * fiber + '.fits')
    files.extend(glob('/cosmos99/nirps/apero-data/nirps_h*online/tmp/*/*' + file_extension + a_b * fiber + '.fits'))
    files.extend(glob('/cosmos99/nirps/apero-data/nirps_h*online/red/*/*' + file_extension + a_b * fiber + '.fits'))

    return files, a_b


def find_limit(logic: str):
    """
    Helper function to find the QCC pass/fail limit value from the logic string (QCCXXXL header key).

    :param logic: str, the value for header key QCCXXXL
    :return: limit, float, the numerical lower or upper limit for the QCC
    """

    # usually the logic string ends with the number/limit, so split string by white space and then
    # try to float from last word, catch ValueError if not a number and proceed until the limit is found
    # tested, should work in general

    found_limit = False
    first_guess = -1
    while not found_limit:
        try:
            limit = float(logic.split()[first_guess])
            found_limit = True
        except ValueError:
            first_guess -= 1

    return limit


def plot_qcc(qcc_value, qcc_value_failed, limit, mjd, mjd_failed, file_type, qcc_name, a_b=False, fiber='A'):
    """
    Function to make and save the figure for an individual QCC, to be called at the end of the other functions.
    """

    qcc_value = np.array(qcc_value)
    qcc_value_failed = np.array(qcc_value_failed)
    mjd = np.array(mjd)
    mjd_failed = np.array(mjd_failed)

    # display the number of files that pass and fail QCC
    n_failed = len(qcc_value_failed)
    n_passed = len(qcc_value)
    results_str = 'Passed: ' + str(n_passed) + '  Failed: ' + str(n_failed)

    # plotting
    # plot passed and failed separately to show with different colors
    # currently green and red, subject to change
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.plot(mjd, qcc_value, 'go')
    ax.plot(mjd_failed, qcc_value_failed, 'ro')
    ax.plot(mjd, limit * np.ones_like(mjd))
    ax.set_xlabel('MJD')
    ax.set_ylabel(qcc_name)
    ax.set_title(file_type + a_b * fiber)
    ax.text(0.99, 0.99, results_str, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes,
            fontsize=10)
    plt.savefig(fig_path + file_type + a_b * fiber + '_' + qcc_name + '.pdf')
    plt.close()

    return


# function to plot a specific qc check for a specific type of file, for all files of that type on disk
def qcc_by_number(file_type: str, qcc_number: str, fiber='A', log=True, plot=True):
    """
    Function to plot a specific qc check for a specific type of file, for all files of that type on disk. Using
    qcc_by_name may be preferred.

    :param file_type: str, the file type for which to list all files, as defined by the dictionaries and the list of
                      files_to_check
    :param qcc_number: str of format 'XXX', e.g. '001'
    :param fiber: 'A' or 'B', if relevant
    :param log: bool; option, to print some progress information
    :param plot: bool; option, in case the function is run only to get the failed files
    :return: list of files that failed the QCC
    """
    # probably better to generate all plots first and then check individual QCCs by name directly for failed files
    # keeping this function for now just in case

    files, a_b = list_files(file_type, fiber=fiber)

    h1 = fits.getheader(files[0])
    # only for PVERSION = 0.7.283 or lower, fixed in 0.7.284
    # (PVERSION is a header key)
    if file_type == 'pp':
        n, v = pp_get_n_v(h1)
    else:
        v = 'V'
        n = 'N'

    try:
        qcc_name = h1['QCC'+qcc_number + n]
    except KeyError:
        print('Unable to find the specified QC check. Input QCC number as a string of three numbers, e.g. \'001\'')
        return

    value = h1['QCC' + qcc_number + v]
    if type(value) != float and type(value) != int:
        print('This QC check does not have a numerical value, therefore it cannot be plotted.')
        return

    logic = h1['QCC' + qcc_number + 'L']
    limit = find_limit(logic)

    qcc_value = []
    qcc_value_failed = []
    mjd = []
    mjd_failed = []
    n_files = len(files)
    failed_files = []
    count = 0
    for file in files:
        h1 = fits.getheader(file)
        try:
            if file_type == 'pp':
                n, v = pp_get_n_v(h1)
            value = float(h1['QCC' + qcc_number + v])
            passed = h1['QCC' + qcc_number + 'P']
            if passed:
                qcc_value.append(value)
                mjd.append(h1['MJD-OBS'])
            else:
                qcc_value_failed.append(value)
                mjd_failed.append(h1['MJD-OBS'])
                failed_files.append(file)
        except KeyError:
            # sometimes files don't have QCC for some reason, ignore those
            if log:
                print('An error occurred with the following file: ' + file)
                print('Skipping to the next file.')
        count += 1
        if log and count % 100 == 0:
            # sanity check for longer run times/more files
            print(str(count)+' out of '+str(n_files)+' files read.')

    if plot:
        plot_qcc(qcc_value, qcc_value_failed, limit, mjd, mjd_failed, file_type, qcc_name, a_b=a_b, fiber=fiber)

    return failed_files


# plotting is an option in case the function is run only to get the failed files
def qcc_by_name(file_type: str, qcc_name: str, fiber='A', log=False, plot=False):
    """
    Function to go back and plot an individual qcc by name, intended to be used after an initial run of plot_all_qcc
    (for one file type) or main (all file types in files_to_check).

    :param file_type: str, the file type for which to list all files, as defined by the dictionaries and the list of
                      files_to_check
    :param qcc_name: str, the exact name of the QCC, as seen on the y-axis of existing plots
    :param fiber: 'A' or 'B', if relevant
    :param log: bool; option, to print some progress information
    :param plot: bool; option, in case the function is run only to get the failed files
    :return: list of files that failed QCCs
    """

    files, a_b = list_files(file_type, fiber=fiber)

    # find the appropriate QCC code/number from its name only
    qcc_code = ''
    h1 = fits.getheader(files[0])
    h1_qcc = h1['QCC*']
    for item in h1_qcc.items():
        if item[1] == qcc_name:
            qcc_code = item[0][:-1]

    logic = h1[qcc_code + 'L']
    limit = find_limit(logic)

    qcc_value = []
    qcc_value_failed = []
    mjd = []
    mjd_failed = []
    n_files = len(files)
    failed_files = []
    count = 0

    for file in files:
        h1 = fits.getheader(file)
        try:
            if file_type == 'pp':
                n, v = pp_get_n_v(h1)
            else:
                v = 'V'
            value = float(h1[qcc_code + v])
            passed = h1[qcc_code + 'P']
            if passed:
                qcc_value.append(value)
                mjd.append(h1['MJD-OBS'])
            else:
                qcc_value_failed.append(value)
                mjd_failed.append(h1['MJD-OBS'])
                failed_files.append(file)
        except KeyError:
            if log:
                print('An error occurred with the following file: ' + file)
                print('Skipping to the next file.')
        count += 1
        if log and count % 100 == 0:
            print(str(count) + ' out of ' + str(n_files) + ' files read.')

    if plot:
        plot_qcc(qcc_value, qcc_value_failed, limit, mjd, mjd_failed, file_type, qcc_name, a_b=a_b, fiber=fiber)

    return failed_files


def plot_all_qcc(file_type: str, fiber='A', log=False):
    """
    Function to plot all the (numerical) QCCs for a specified file type. Always produces plots. The destination of
    plots can be changed at the beginning of this file (qc_over_time.py).

    :param file_type: str, the file type for which to list all files, as defined by the dictionaries and the list of
                      files_to_check
    :param fiber: 'A' or 'B', if relevant.
    :param log: bool; option, to print some progress information
    :return: list of all QCC names (qcc_list_names), lists of all failed files for every QCC (failed_files_global)
    """
    qcc_list_num = eval(file_type)['numerical']
    qcc_list_names = eval(file_type)['qcc_names']

    files, a_b = list_files(file_type, fiber=fiber)
    n_files = len(files)

    failed_files_global = []
    qcc_value_global = []
    qcc_value_failed_global = []
    mjd_global = []
    mjd_failed_global = []
    limit_global = []

    h1 = fits.getheader(files[0])
    for num in qcc_list_num:
        qcc_value_global.append([])
        qcc_value_failed_global.append([])
        mjd_global.append([])
        mjd_failed_global.append([])
        failed_files_global.append([])

        logic = h1['QCC00' + num + 'L']
        limit = find_limit(logic)
        limit_global.append(limit)

    count = 0
    for file in files:
        h1 = fits.getheader(file)
        if file_type == 'pp':
            n, v = pp_get_n_v(h1)
        else:
            n = 'N'
            v = 'V'

        qcc_count = 0
        for num in qcc_list_num:
            qcc_number = '00' + num
            try:
                value = float(h1['QCC' + qcc_number + v])
                passed = h1['QCC' + qcc_number + 'P']
                if passed:
                    qcc_value_global[qcc_count].append(value)
                    mjd_global[qcc_count].append(h1['MJD-OBS'])
                else:
                    qcc_value_failed_global[qcc_count].append(value)
                    mjd_failed_global[qcc_count].append(h1['MJD-OBS'])
                    failed_files_global[qcc_count].append(file)
            except KeyError:
                if log:
                    print('An error occurred with the following file: ' + file)
                    print('Skipping to the next file.')
            qcc_count += 1

        count += 1
        if log and count % 100 == 0:
            print(str(count) + ' out of ' + str(n_files) + ' files read.')

    i = 0
    for qcc_name in qcc_list_names:
        plot_qcc(qcc_value_global[i], qcc_value_failed_global[i], limit_global[i], mjd_global[i], mjd_failed_global[i],
                 file_type, qcc_name, a_b=a_b, fiber=fiber)
        i += 1

    return qcc_list_names, failed_files_global


if __name__ == "__main__":
    # plot all QCCs for all files with log on with one short line of code
    # for fiber A

    print('Plotting all QCC for relevant file types')
    for file_type in files_to_check:
        print('File type: {}'.format(file_type))
        plot_all_qcc(file_type, fiber='A', log=True)
