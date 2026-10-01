# APERO Astrometrics Export Helper

This folder includes `export_astrometrics_yamls.py`, a small utility that
copies the astrometric yaml assets used by
`apero.core.drs_astrometrics` into a local directory.

The script uses `apero.dev.get_base_params(...)` to bootstrap temporary APERO
configuration, so an existing APERO profile is **not** required.
If the astrometric assets are not already present under `PATH.ASSETS`, the
script attempts to refresh/download APERO assets via
`apero.tools.module.setup.drs_assets`.

## Requirements

- Python environment with `apero-core` and `apero-drs` importable

## Quick run

```bash
PYTHONPATH="/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-core:/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs" \
python /scratch2/spirou/drs-bin/apero-utils/general/apero_astrometrics/export_astrometrics_yamls.py
```

By default, files are copied into `./apero_astremetrics` (note the directory
`./apero_astrometrics`.

The export root also gets a user README at `README.md` describing:

- what `verified/`, `pending/`, and `rejected/` mean,
- and key yaml fields including `APERO_NAME` and `SIMBAD_NAME`.

## Useful options

```bash
PYTHONPATH="/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-core:/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs" \
python /scratch2/spirou/drs-bin/apero-utils/general/apero_astrometrics/export_astrometrics_yamls.py \
    --instrument SPIROU \
    --output-dir /tmp/apero_astrometrics \
    --force-asset-download \
    --clean-output
```

