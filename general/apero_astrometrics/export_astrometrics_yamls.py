#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Export APERO astrometric yaml files into a local directory.

This script uses ``apero.dev`` to bootstrap a temporary configuration,
so no existing APERO profile is required. It then finds the astrometric
asset directory used by ``apero.core.drs_astrometrics`` and copies every
``*.yaml`` entry into a local output directory.

If assets are missing, the script can refresh them via APERO setup tools.
You can also force a non-local tarfile path to exercise the download logic.

Examples
--------
Copy yaml files into the default ``apero_astrometrics`` directory::

    python export_astrometrics_yamls.py

Choose an instrument and destination::

    python export_astrometrics_yamls.py --instrument NIRPS_HE \
        --output-dir /tmp/apero_astrometrics --clean-output

Force APERO to retrieve assets without reading a local assets tarfile::

    python export_astrometrics_yamls.py --force-asset-download
"""

import argparse
import pathlib
import shutil
import sys
from typing import Any, Optional, Sequence, Tuple

from apero.core import drs_astrometrics
from apero.dev import get_base_params


# =============================================================================
# Define functions
# =============================================================================
def get_args(argv: Sequence[str]) -> argparse.Namespace:
    """Parse command-line arguments.

    Parameters
    ----------
    argv : Sequence[str]
        Command-line arguments excluding the executable name.

    Returns
    -------
    argparse.Namespace
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description=('Copy APERO astrometric yaml assets into a local '
                     'directory without requiring an APERO profile.')
    )
    parser.add_argument(
        '--instrument',
        default='SPIROU',
        help=('APERO instrument used to bootstrap dev parameters '
              '(default: SPIROU).')
    )
    parser.add_argument(
        '--output-dir',
        default='apero_astrometrics',
        help=('Destination directory for copied yaml files '
              '(default: ./apero_astrometrics).')
    )
    parser.add_argument(
        '--clean-output',
        action='store_true',
        help='Remove the output directory before copying files.'
    )
    parser.add_argument(
        '--force-asset-download',
        action='store_true',
        help=('Force APERO assets update to avoid re-using a local assets '
              'tar file, so the download branch is exercised.')
    )
    return parser.parse_args(list(argv))


def _resolve_astrometrics_dir(
    instrument: str,
    force_asset_download: bool,
) -> pathlib.Path:
    """Return the astrometric asset directory for an instrument.

    Parameters
    ----------
    instrument : str
        APERO instrument name passed to ``apero.dev.get_base_params``.

    Returns
    -------
    pathlib.Path
        Absolute path to ``PATH.ASSETS/astrometrics``.

    Raises
    ------
    FileNotFoundError
        Raised when the astrometric assets directory does not exist.
    """
    # Use apero.dev so we do not require a pre-existing APERO profile.
    params = get_base_params(instrument)
    assets_root = pathlib.Path(str(params['PATH.ASSETS'])).resolve()
    astrom_dir = _ensure_astrometric_assets(
        params=params,
        assets_root=assets_root,
        force_asset_download=force_asset_download,
    )
    if not astrom_dir.is_dir():
        emsg = 'Astrometric assets directory not found: {0}'
        raise FileNotFoundError(emsg.format(astrom_dir))
    return astrom_dir


def _ensure_astrometric_assets(
    params: Any,
    assets_root: pathlib.Path,
    force_asset_download: bool,
) -> pathlib.Path:
    """Ensure astrometric assets are available on disk.

    Parameters
    ----------
    params : Any
        APERO parameter dictionary from ``apero.dev.get_base_params``.
    assets_root : pathlib.Path
        Root assets directory (``PATH.ASSETS``).

    Returns
    -------
    pathlib.Path
        Path to the astrometrics asset directory.

    Raises
    ------
    RuntimeError
        Raised when assets cannot be retrieved.
    """
    from apero.tools.module.setup import drs_assets

    astrom_dir = assets_root / drs_astrometrics.ASTROM_SUBDIR
    if astrom_dir.is_dir() and (not force_asset_download):
        yaml_files = drs_astrometrics.iter_yaml_files(str(astrom_dir))
        if len(yaml_files) > 0:
            return astrom_dir

    # If checks fail or data are missing, ask APERO setup helper to refresh.
    try:
        update_needed = bool(drs_assets.check_local_assets(params))
    except Exception:
        update_needed = True

    backup_tar = None
    if force_asset_download:
        backup_tar = _disable_local_assets_tar(params, assets_root)

    if force_asset_download or update_needed or not astrom_dir.is_dir():
        try:
            drs_assets.update_local_assets(params)
        except Exception as exc:
            _restore_assets_tar(backup_tar)
            emsg = 'Failed to retrieve APERO assets: {0}'
            raise RuntimeError(emsg.format(exc)) from exc
        _finalize_assets_tar_backup(backup_tar)

    yaml_files = drs_astrometrics.iter_yaml_files(str(astrom_dir))
    if len(yaml_files) == 0:
        emsg = 'No astrometric yaml files found in {0}'
        raise RuntimeError(emsg.format(astrom_dir))
    return astrom_dir


def _disable_local_assets_tar(
    params: Any,
    assets_root: pathlib.Path,
) -> Optional[Tuple[pathlib.Path, pathlib.Path]]:
    """Temporarily move local assets tar so APERO must use download logic.

    Parameters
    ----------
    params : Any
        APERO parameter dictionary.
    assets_root : pathlib.Path
        Root assets directory (``PATH.ASSETS``).

    Returns
    -------
    tuple or None
        ``(tar_path, backup_path)`` when a local tar existed and was moved,
        otherwise ``None``.
    """
    _ = assets_root
    tar_path = _get_assets_tar_path(params=params)
    if tar_path is None or (not tar_path.exists()):
        return None
    backup_path = pathlib.Path(str(tar_path) + '.skip_local')
    if backup_path.exists():
        backup_path.unlink()
    tar_path.rename(backup_path)
    return tar_path, backup_path


def _restore_assets_tar(
    backup_tar: Optional[Tuple[pathlib.Path, pathlib.Path]],
) -> None:
    """Restore a temporarily moved assets tar file after update failure."""
    if backup_tar is None:
        return
    tar_path, backup_path = backup_tar
    if backup_path.exists() and (not tar_path.exists()):
        backup_path.rename(tar_path)


def _finalize_assets_tar_backup(
    backup_tar: Optional[Tuple[pathlib.Path, pathlib.Path]],
) -> None:
    """Remove backup tar after a successful forced refresh."""
    if backup_tar is None:
        return
    _tar_path, backup_path = backup_tar
    if backup_path.exists():
        backup_path.unlink()


def _get_assets_tar_path(
    params: Any,
) -> Optional[pathlib.Path]:
    """Return the local APERO assets tar path from checksum metadata."""
    from apero.base import base as apero_base
    from aperocore.base import base
    from apero.utils import drs_data

    try:
        assets_rel = str(params['IPATH.RESET_ASSETS'])
        cdata_rel = str(params['IPATH.CDATA'])
        assets_abs = pathlib.Path(
            drs_data.construct_path(params, '', assets_rel))
        cdata_abs = pathlib.Path(drs_data.construct_path(params, '', cdata_rel))
        checksum_path = cdata_abs / apero_base.CHECKSUM_FILE
        yaml_dict = base.load_yaml(str(checksum_path))
        tar_name = str(yaml_dict['setup']['tarfile'])
    except Exception:
        return None
    return assets_abs / tar_name


def _write_export_readme(output_dir: pathlib.Path) -> None:
    """Write a user-facing README describing exported astrometric data."""
    lines = [
        '# APERO Astrometrics YAML Export',
        '',
        'This directory contains APERO astrometric object entries.',
        'Each object is stored as one `.yaml` file.',
        '',
        '## Sub-directories',
        '',
        '- `verified/`: entries validated and ready for production use.',
        '- `pending/`: entries waiting for human validation.',
        '- `rejected/`: entries explicitly rejected or excluded.',
        '',
        'Some installs may also include legacy `.yaml` files directly in the',
        'root of this directory.',
        '',
        '## YAML structure (key points)',
        '',
        '- `APERO_NAME`: canonical APERO object name used as primary key.',
        '- `SIMBAD_NAME`: SIMBAD-resolved object name when available.',
        '',
        'Typical files include many other fields (coordinates, proper motion,',
        'parallax, aliases, and metadata).',
        '',
    ]
    readme_path = output_dir / 'README.md'
    readme_path.write_text('\n'.join(lines), encoding='utf-8')


def export_yaml_files(
    instrument: str,
    output_dir: pathlib.Path,
    clean_output: bool,
    force_asset_download: bool,
) -> Tuple[int, pathlib.Path]:
    """Copy all astrometric yaml files into ``output_dir``.

    Parameters
    ----------
    instrument : str
        APERO instrument used for developer bootstrap.
    output_dir : pathlib.Path
        Directory where yaml files are copied.
    clean_output : bool
        If True, delete ``output_dir`` first.

    Returns
    -------
    tuple
        ``(count, source_dir)`` where ``count`` is the number of copied
        yaml files and ``source_dir`` is the resolved source directory.
    """
    source_dir = _resolve_astrometrics_dir(
        instrument=instrument,
        force_asset_download=force_asset_download,
    )

    # Optionally start from a clean destination tree.
    if clean_output and output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Use APERO helper so status sub-dirs are handled exactly like pipeline.
    yaml_files = drs_astrometrics.iter_yaml_files(str(source_dir))
    for yaml_file in yaml_files:
        source_file = pathlib.Path(yaml_file)
        try:
            relative = source_file.relative_to(source_dir)
        except ValueError:
            # Fallback for unexpected paths; keep filename only.
            relative = pathlib.Path(source_file.name)
        target_file = output_dir / relative
        target_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_file, target_file)

    # Add user documentation to the export root.
    _write_export_readme(output_dir)

    return len(yaml_files), source_dir


def main(argv: Sequence[str]) -> int:
    """Run the astrometric yaml export command.

    Parameters
    ----------
    argv : Sequence[str]
        Command-line arguments excluding the executable name.

    Returns
    -------
    int
        Process exit code (0 on success, non-zero on failure).
    """
    args = get_args(argv)
    output_dir = pathlib.Path(args.output_dir).expanduser().resolve()

    try:
        nfiles, source_dir = export_yaml_files(
            instrument=str(args.instrument),
            output_dir=output_dir,
            clean_output=bool(args.clean_output),
            force_asset_download=bool(args.force_asset_download),
        )
    except Exception as exc:
        print('Error: {0}'.format(exc), file=sys.stderr)
        return 1

    print('Astrometric source: {0}'.format(source_dir))
    print('Output directory:   {0}'.format(output_dir))
    print('Copied yaml files:  {0}'.format(nfiles))
    return 0


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':
    main(sys.argv[1:])

# =============================================================================
# End of code
# =============================================================================

