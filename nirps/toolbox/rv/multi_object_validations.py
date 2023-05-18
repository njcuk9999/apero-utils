"""
Plot RV rms vs median uncertainties and BERV coverage for NIRPS HA and NIRPS HE lbl data products
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

ha_files = glob.glob('/cosmos99/nirps/lbl-data/nirps_ha_07276_online/lblrdb/lbl_*.rdb')
he_files = glob.glob('/cosmos99/nirps/lbl-data/nirps_he_07276_online/lblrdb/lbl_*.rdb')


def multi_object_validation_plot(files, fiber):

	# Empty arrays
	rv_rms = np.zeros(len(files))
	e_rv_median = np.zeros(len(files))
	berv_cov = np.zeros(len(files))
	color = np.zeros(len(files)).astype(bool)
	name = np.zeros(len(files)).astype(str)

	i = 0
	for file in files:
		tbl = Table.read(file)
		# Relevant time series
		rv = tbl['vrad']
		rv -= np.median(rv)
		e_rv = tbl['svrad']
		berv = tbl['BERV']

		# Store str for target_template to label the points in the figure
		identification = file.split('lbl_')[-1]
		identification = identification.split('.rdb')[0]
		name[i] = identification

		# Verify that there is no inf of NaN values in the timeseries. If yes the marker will be red
		valid = np.isfinite(rv)
		if np.sum(np.invert(valid)) > 0: 
			color[i] = True
		else: 
			color[i] = False

		# Calculate the quantities of interest
		rv_rms[i] = np.sqrt(np.mean(rv[valid]**2))    # RV RMS
		e_rv_median[i] = np.nanmedian(e_rv)           # Median RV uncertainty
		berv_cov[i] = np.max(berv) - np.min(berv)     # BERV coverage
		i += 1

	# Plot figure
	fig, ax = plt.subplots(nrows=2, ncols=2, sharey='row', sharex='col', figsize=(16, 13))

	# RV rms vs median uncertainties (ALL DATA)
	ax[0, 0].errorbar(e_rv_median, rv_rms, ls='None', marker='o', color='k', capsize=3)
	ax[0, 0].errorbar(e_rv_median[color], rv_rms[color], ls='None', marker='o', color='r', capsize=3)  # time series with non-finite values
	for k in range(len(name)):
		ax[0, 0].text(e_rv_median[k], rv_rms[k], name[k], color='grey', alpha=.5, zorder=0)
	ax[0, 0].plot(e_rv_median, e_rv_median, ls='-', marker='None', color='C1')
	ax[0, 0].set_ylabel('RV RMS (m/s)', fontsize=17)

	# RV rms vs median uncertainties (ZOOM-IN / WITHOUT OUTLIERS)
	median_rv_rms = np.median(rv_rms)
	std_rv_rms = np.std(rv_rms, ddof=1)
	mask_outliers = np.abs(rv_rms-median_rv_rms) < std_rv_rms

	ax[1, 0].errorbar(e_rv_median[mask_outliers], rv_rms[mask_outliers], ls='None', marker='o', color='k', capsize=3)
	ax[1, 0].errorbar(e_rv_median[mask_outliers*color], rv_rms[mask_outliers*color], ls='None', marker='o', color='r', capsize=3)
	for k in range(len(name[mask_outliers])):
		ax[1, 0].text(e_rv_median[mask_outliers][k], rv_rms[mask_outliers][k], name[mask_outliers][k], color='grey', alpha=.5, zorder=0)
	ax[1, 0].plot(e_rv_median[mask_outliers], e_rv_median[mask_outliers], ls='-', marker='None', color='C1')
	ax[1, 0].set_ylabel('RV RMS (m/s)', fontsize=17)
	ax[1, 0].set_xlabel('Median RV uncertainty (m/s)', fontsize=17)

	# RV rms vs BERV coverage (ALL DATA)
	ax[0, 1].errorbar(berv_cov, rv_rms, ls='None', marker='o', color='k', capsize=3)
	ax[0, 1].errorbar(berv_cov[color], rv_rms[color], ls='None', marker='o', color='r', capsize=3)  # time series with non-finite values
	for k in range(len(name)):
		ax[0, 1].text(berv_cov[k], rv_rms[k], name[k], color='grey', alpha=.5, zorder=0)

	# RV rms vs BERV coverage (ZOOM-IN / WITHOUT OUTLIERS)
	ax[1, 1].errorbar(berv_cov[mask_outliers], rv_rms[mask_outliers], ls='None', marker='o', color='k', capsize=3)    # , xerr=e_berv_median, yerr=e_rv_rms
	ax[1, 1].errorbar(berv_cov[mask_outliers*color], rv_rms[mask_outliers*color], ls='None', marker='o', color='r', capsize=3)  # time series with non-finite values  , xerr=e_berv_median[color], yerr=e_rv_rms[color]
	for k in range(len(name[mask_outliers])):
		ax[1, 1].text(berv_cov[mask_outliers][k], rv_rms[mask_outliers][k], name[mask_outliers][k], color='grey', alpha=.5, zorder=0)
	ax[1, 1].set_xlabel(r'$\Delta$ BERV (km/s)', fontsize=17)

	plt.tick_params(axis='both', direction='in', bottom=True, top=True, left=True, right=True)
	plt.tight_layout()
	fig.suptitle(fiber, x=0.97, fontsize=15)



# HA
multi_object_validation_plot(ha_files, 'HA')
plt.savefig('HA_multi-object_validation.pdf')

# HE
multi_object_validation_plot(he_files, 'HE')
plt.savefig('HE_multi-object_validation.pdf')