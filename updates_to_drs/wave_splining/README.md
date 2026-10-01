# Wave-grid spline utilities

## `apply_ohpca.py`

Spline OH PCA vectors from an old wavelength grid onto a new one.

### What it does

- Loads `wave_old`, `wave_new`, and `pca_old` from FITS files.
- Checks that `pca_old.shape[1] == wave_old.size`.
- Reshapes each PCA component to `(n_orders, n_pix_old)`.
- Spline-interpolates each order onto `wave_new`.
- Flattens back to `(n_components, wave_new.size)`.
- Writes output to `pca_new`.
- If `pca_new` already exists, makes a timestamped backup in the script
  directory before writing.

### Quick run

```bash
python /scratch2/spirou/drs-bin/apero-utils/updates_to_drs/wave_splining/apply_ohpca.py
```

### Dry run (shape checks + interpolation only)

```bash
python /scratch2/spirou/drs-bin/apero-utils/updates_to_drs/wave_splining/apply_ohpca.py --dry-run
```

### Synthetic self-test

```bash
python /scratch2/spirou/drs-bin/apero-utils/updates_to_drs/wave_splining/apply_ohpca.py --self-test
```

### Interactive comparison plots

Show old vs new `(wavelength, PCA)` curves for selected components:

```bash
python /scratch2/spirou/drs-bin/apero-utils/updates_to_drs/wave_splining/apply_ohpca.py --plot --plot-components 0,1,5
```

This opens matplotlib windows (`plt.show()`) and does not save figures.
Each figure includes all spectral orders, with a top panel showing old/new
curves and a bottom panel showing `(new interpolated to old grid) - old`.

## `apply_ohsky_model.py`

Spline all 2-D sky-model extensions (including `SCI_SKY`, `CAL_SKY`,
`WAVE`, `REG_ID`, `WEIGHTS`, `GRADIENT`) onto a target `wave_new` grid.

- Uses `WAVE` extension in the input model as the old wavelength grid.
- Uses instrument defaults for target grids when `--wave-new` is omitted:
  `DEFAULT_WAVE_NEW_HA` and `DEFAULT_WAVE_NEW_HE`.
- No external `wave_old` input is required.
- Updates each matching 2-D image extension from old pixel sampling
  (e.g. 4088) to `wave_new.shape[1]` (e.g. 8176).
- Preserves integer dtype for discrete extensions like `REG_ID`.
- If writing in place, creates a timestamped backup in the script directory.

### Example: NIRPS-HA model

```bash
python /scratch2/spirou/drs-bin/apero-utils/updates_to_drs/wave_splining/apply_ohsky_model.py \
--model-in /scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/apero/apero-assets/nirps_ha/reset/telludb/sky_model_ha.fits \
--wave-new /path/to/new_wave.fits \
--dry-run --plot
```

### Example: NIRPS-HE model

```bash
python /scratch2/spirou/drs-bin/apero-utils/updates_to_drs/wave_splining/apply_ohsky_model.py \
--model-in /scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/apero/apero-assets/nirps_he/reset/telludb/sky_model_he.fits \
--wave-new /path/to/new_wave.fits \
--dry-run --plot
```

### Synthetic self-test

```bash
python /scratch2/spirou/drs-bin/apero-utils/updates_to_drs/wave_splining/apply_ohsky_model.py --self-test
```

- Supports table-format PCA files with shape `(n_pix_flat, n_col)` where
  column `0` is the old wavelength grid and columns `1..N` are PCA values.
- In that format, column `0` is set to the new wave grid directly and only
  columns `1..N` are spline-interpolated.
