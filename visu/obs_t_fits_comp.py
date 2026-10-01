from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from tqdm import tqdm
from datetime import datetime

# Walk directories and follow symlinks
base = Path("/cosmos99/nirps/apero-data/nirps_he_online/objects/HD195094")
datasets = [str(p) for p in base.rglob('*t.fits')]

data_list = []

# Step 1: Collect metadata and flux metrics across all files
print("Loading headers and measuring fluxes...")
for filename in tqdm(datasets[::2]):
    with fits.open(filename) as hdul:
        hdr = hdul[0].header
        wave = hdul['WaveA'].data
        flux = hdul['FluxA'].data

    date_str = hdr['DATE']
    # Convert date string to datetime object for proper sorting/plotting
    try:
        date_dt = datetime.fromisoformat(date_str)
    except ValueError:
        date_dt = datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%S.%f')

    mean_airmass = np.nanmean([hdr['HIERARCH ESO TEL AIRM START'], hdr['HIERARCH ESO TEL AIRM END']])
    seeing = float(hdr['HIERARCH ESO INS2 AOS ATM SEEING'])

    # Define approximate wavelength boundaries for J-band (1100 - 1400 nm) and H-band (1500 - 1800 nm)
    j_mask = (wave >= 1100) & (wave <= 1400)
    h_mask = (wave >= 1500) & (wave <= 1800)

    median_flux_j = np.nanmedian(flux[j_mask]) if np.any(j_mask) else np.nan
    median_flux_h = np.nanmedian(flux[h_mask]) if np.any(h_mask) else np.nan

    data_list.append({
        'filename': filename,
        'wave': wave,
        'flux': flux,
        'date_str': date_str,
        'date': date_dt,
        'airmass': mean_airmass,
        'seeing': seeing,
        'j_flux': median_flux_j,
        'h_flux': median_flux_h
    })

# Step 2: Select spectra based ONLY on date (earliest, median, latest)
dates = [d['date'] for d in data_list]
sorted_date_idx = np.argsort(dates)

earliest_idx = sorted_date_idx[0]
median_idx = sorted_date_idx[len(sorted_date_idx) // 2]
latest_idx = sorted_date_idx[-1]

date_key_indices = [
    (earliest_idx, "Earliest Date"),
    (median_idx, "Median Date"),
    (latest_idx, "Latest Date")
]

# Distinct marker per highlighted date, reused across all bottom plots
date_markers = {
    "Earliest Date": "^",
    "Median Date": "s",
    "Latest Date": "*"
}

# Step 3: Set up grid layout
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, height_ratios=[1.5, 1])

ax_top = fig.add_subplot(gs[0, :])     # Full width top plot
ax_bottom1 = fig.add_subplot(gs[1, 0]) # Seeing plot
ax_bottom2 = fig.add_subplot(gs[1, 1]) # Airmass plot
ax_bottom3 = fig.add_subplot(gs[1, 2]) # Date plot

# Plot selected date spectra on top
for idx, tag in date_key_indices:
    d = data_list[idx]
    label = f"{tag}: {d['date_str'][:10]} | Air: {d['airmass']:.2f} | See: {d['seeing']:.2f}"
    ax_top.plot(d['wave'].ravel(), d['flux'].ravel(), label=label, alpha=0.7, lw=0.8)

ax_top.set_title("Spectra: Earliest, Median, and Latest Observation Dates")
ax_top.set_xlabel("Wavelength [nm]")
ax_top.set_ylabel("Flux")
ax_top.legend(loc='upper right', fontsize=8)

# Prepare series data for trend plots
seeings = [d['seeing'] for d in data_list]
airmasses = [d['airmass'] for d in data_list]
j_fluxes = [d['j_flux'] for d in data_list]
h_fluxes = [d['h_flux'] for d in data_list]

# Overlay the earliest/median/latest points on a bottom plot using
# distinct markers so they can be spotted against the full scatter
def highlight_dates(ax, xs):
    for idx, tag in date_key_indices:
        marker = date_markers[tag]
        ax.scatter(xs[idx], j_fluxes[idx], marker=marker, s=120,
                   facecolors='none', edgecolors='k', linewidths=1.3,
                   zorder=5, label=tag)
        ax.scatter(xs[idx], h_fluxes[idx], marker=marker, s=120,
                   facecolors='none', edgecolors='k', linewidths=1.3,
                   zorder=5)

# Bottom Plot 1: Flux vs Seeing (Scatter)
ax_bottom1.scatter(seeings, j_fluxes, alpha=0.7, s=20, label='J band')
ax_bottom1.scatter(seeings, h_fluxes, alpha=0.7, s=20, label='H band')
highlight_dates(ax_bottom1, seeings)
ax_bottom1.set_xlabel("Seeing")
ax_bottom1.set_ylabel("Median Flux")
ax_bottom1.set_title("Median Flux vs Seeing")
ax_bottom1.legend(fontsize=7)

# Bottom Plot 2: Flux vs Airmass (Scatter)
ax_bottom2.scatter(airmasses, j_fluxes, alpha=0.7, s=20, label='J band')
ax_bottom2.scatter(airmasses, h_fluxes, alpha=0.7, s=20, label='H band')
highlight_dates(ax_bottom2, airmasses)
ax_bottom2.set_xlabel("Airmass")
ax_bottom2.set_ylabel("Median Flux")
ax_bottom2.set_title("Median Flux vs Airmass")
ax_bottom2.legend(fontsize=7)

# Bottom Plot 3: Flux vs Date (Scatter)
ax_bottom3.scatter(dates, j_fluxes, alpha=0.7, s=20, label='J band')
ax_bottom3.scatter(dates, h_fluxes, alpha=0.7, s=20, label='H band')
highlight_dates(ax_bottom3, dates)
ax_bottom3.set_xlabel("Date")
ax_bottom3.set_ylabel("Median Flux")
ax_bottom3.set_title("Median Flux vs Date")
ax_bottom3.tick_params(axis='x', rotation=30)
ax_bottom3.legend(fontsize=7)

plt.tight_layout()
plt.show()