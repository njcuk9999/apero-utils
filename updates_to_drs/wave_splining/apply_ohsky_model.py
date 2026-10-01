#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Spline NIRPS sky-model extensions onto a new wavelength grid.

This script updates 2-D sky-model image extensions from their original
pixel sampling (typically 4088) to the target sampling defined by a
``wave_new`` file (typically 8176 pixels).

The model FITS file already contains its old wavelength grid in extension
``WAVE``, so no external ``wave_old`` input is required.
"""

import argparse
import datetime
import os
import shutil
import tempfile

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy import interpolate

# =============================================================================
# Define variables
# =============================================================================
DEFAULT_MODEL_HA = (
    '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/apero/'
    'apero-assets/nirps_ha/reset/telludb/sky_model_ha.fits'
)
DEFAULT_MODEL_HE = (
    '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/apero/'
    'apero-assets/nirps_he/reset/telludb/sky_model_he.fits'
)
DEFAULT_WAVE_NEW_HA = (
    '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/apero/'
    'apero-assets/nirps_ha/calib/static_wave_ref_A.fits'
)
DEFAULT_WAVE_NEW_HE = (
    '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/apero/'
    'apero-assets/nirps_he/calib/static_wave_ref_A.fits'
)


# =============================================================================
# Define functions
# =============================================================================
def parse_name_list(value: str) -> list[str]:
    """Parse a comma-separated list of extension names.

    :param value: string like ``'SCI_SKY,CAL_SKY'``
    :return: uppercase extension-name list
    """
    tokens = [token.strip().upper() for token in value.split(',')]
    return [token for token in tokens if token]


def resolve_wave_new_path(model_in_path: str, wave_new_arg: str | None) -> str:
    """Resolve wave_new path from explicit argument or instrument defaults.

    :param model_in_path: input model file path
    :param wave_new_arg: optional explicit wave_new file path
    :return: chosen wave_new file path
    :raises ValueError: if instrument cannot be inferred
    """
    if wave_new_arg is not None:
        return wave_new_arg

    model_path = model_in_path.lower()
    if 'nirps_ha' in model_path or 'sky_model_ha' in model_path:
        return DEFAULT_WAVE_NEW_HA
    if 'nirps_he' in model_path or 'sky_model_he' in model_path:
        return DEFAULT_WAVE_NEW_HE

    emsg = (
        'Cannot infer wave_new path from model filename. '
        'Use --wave-new explicitly.'
    )
    raise ValueError(emsg)


def first_data_hdu_index(hdul: fits.HDUList) -> int:
    """Return index of first HDU containing image data.

    :param hdul: opened FITS HDU list
    :return: index of first HDU with non-None data
    :raises ValueError: if no HDU has data
    """
    for index, hdu in enumerate(hdul):
        if hdu.data is not None:
            return index
    raise ValueError('No image data found in FITS file.')


def load_wave_new(path: str) -> np.ndarray:
    """Load new wavelength grid as 2-D float array.

    :param path: wave_new FITS file path
    :return: 2-D wavelength grid
    :raises ValueError: if data is not 2-D
    """
    with fits.open(path) as hdul:
        hdu_index = first_data_hdu_index(hdul)
        wave = np.array(hdul[hdu_index].data, dtype=float)
    if wave.ndim != 2:
        raise ValueError('wave_new must be 2-D.')
    return wave


def backup_file(source_path: str, backup_dir: str) -> str:
    """Create timestamped backup of ``source_path`` in ``backup_dir``.

    :param source_path: file to back up
    :param backup_dir: destination backup directory
    :return: backup file path
    """
    os.makedirs(backup_dir, exist_ok=True)
    stamp = datetime.datetime.now(datetime.UTC).strftime('%Y%m%dT%H%M%SZ')
    base_name = os.path.basename(source_path)
    backup_name = f'{base_name}.backup_{stamp}'
    backup_path = os.path.join(backup_dir, backup_name)
    shutil.copy2(source_path, backup_path)
    return backup_path


def get_wave_extension_index(hdul: fits.HDUList) -> int:
    """Return HDU index for extension named ``WAVE``.

    :param hdul: opened FITS HDU list
    :return: extension index
    :raises ValueError: if ``WAVE`` extension is missing
    """
    for index, hdu in enumerate(hdul):
        if hdu.name == 'WAVE' and hdu.data is not None:
            return index
    raise ValueError('Input model FITS must contain image extension WAVE.')


def orient_orders_pixels(
    array_2d: np.ndarray,
    n_orders: int,
    n_pix: int,
) -> tuple[np.ndarray, bool]:
    """Orient 2-D array as ``(n_orders, n_pix)``.

    :param array_2d: input 2-D array
    :param n_orders: expected order count
    :param n_pix: expected pixel count
    :return: tuple(oriented array, transpose_used)
    :raises ValueError: if shape does not match expected dimensions
    """
    if array_2d.shape == (n_orders, n_pix):
        return np.array(array_2d), False
    if array_2d.shape == (n_pix, n_orders):
        return np.array(array_2d).T, True
    emsg = (
        'Array shape incompatible with expected order/pixel dimensions: '
        f'{array_2d.shape} vs ({n_orders}, {n_pix}) or ({n_pix}, {n_orders})'
    )
    raise ValueError(emsg)


def interpolate_one_order(
    wave_old_order: np.ndarray,
    data_old_order: np.ndarray,
    wave_new_order: np.ndarray,
    is_discrete: bool,
) -> np.ndarray:
    """Interpolate one order from old wavelength grid to new wavelength grid.

    :param wave_old_order: old wavelength values for one order
    :param data_old_order: old data values for one order
    :param wave_new_order: new wavelength values for one order
    :param is_discrete: use nearest-neighbor if True, spline otherwise
    :return: interpolated values for one order on new wave grid
    """
    finite = np.isfinite(wave_old_order) & np.isfinite(data_old_order)
    x_valid = wave_old_order[finite]
    y_valid = data_old_order[finite]

    if x_valid.size < 2:
        return np.full_like(wave_new_order, np.nan, dtype=float)

    sort_index = np.argsort(x_valid)
    x_sorted = x_valid[sort_index]
    y_sorted = y_valid[sort_index]

    x_unique, unique_index = np.unique(x_sorted, return_index=True)
    y_unique = y_sorted[unique_index]

    if x_unique.size < 2:
        return np.full_like(wave_new_order, np.nan, dtype=float)

    if is_discrete:
        interp = interpolate.interp1d(
            x_unique,
            y_unique,
            kind='nearest',
            bounds_error=False,
            fill_value=(y_unique[0], y_unique[-1]),
            assume_sorted=True,
        )
        return np.array(interp(wave_new_order), dtype=float)

    degree = min(3, x_unique.size - 1)
    spline = interpolate.InterpolatedUnivariateSpline(
        x_unique,
        y_unique,
        k=degree,
        ext=3,
    )
    return np.array(spline(wave_new_order), dtype=float)


def interpolate_extension(
    data_old: np.ndarray,
    wave_old: np.ndarray,
    wave_new: np.ndarray,
    is_discrete: bool,
) -> np.ndarray:
    """Interpolate a full extension from old to new wave grid.

    :param data_old: old extension data, shape ``(n_orders, n_pix_old)``
    :param wave_old: old wave data, shape ``(n_orders, n_pix_old)``
    :param wave_new: new wave data, shape ``(n_orders, n_pix_new)``
    :param is_discrete: nearest-neighbor interpolation flag
    :return: interpolated extension, shape ``(n_orders, n_pix_new)``
    """
    n_orders = wave_old.shape[0]
    n_pix_new = wave_new.shape[1]
    out = np.full((n_orders, n_pix_new), np.nan, dtype=float)

    for order_num in range(n_orders):
        out[order_num] = interpolate_one_order(
            wave_old[order_num],
            data_old[order_num],
            wave_new[order_num],
            is_discrete,
        )

    return out


def plot_extension_comparison(
    ext_name: str,
    wave_old: np.ndarray,
    wave_new: np.ndarray,
    data_old: np.ndarray,
    data_new: np.ndarray,
) -> None:
    """Plot all-order comparison and difference for one extension.

    :param ext_name: extension name for figure title
    :param wave_old: old wave array, shape ``(n_orders, n_pix_old)``
    :param wave_new: new wave array, shape ``(n_orders, n_pix_new)``
    :param data_old: old extension array, shape ``(n_orders, n_pix_old)``
    :param data_new: new extension array, shape ``(n_orders, n_pix_new)``
    :return: None
    """
    n_orders = wave_old.shape[0]
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=False)
    ax_top = axes[0]
    ax_bot = axes[1]

    for order_num in range(n_orders):
        w_old = wave_old[order_num]
        w_new = wave_new[order_num]
        y_old = data_old[order_num]
        y_new = data_new[order_num]

        ax_top.plot(
            w_old,
            y_old,
            linestyle='--',
            color='0.50',
            alpha=0.35,
            linewidth=0.8,
        )
        ax_top.plot(
            w_new,
            y_new,
            linestyle='-',
            color='tab:blue',
            alpha=0.25,
            linewidth=0.8,
        )

        order_sort = np.argsort(w_new)
        w_new_sort = w_new[order_sort]
        y_new_sort = y_new[order_sort]
        w_new_unique, unique_index = np.unique(w_new_sort, return_index=True)
        y_new_unique = y_new_sort[unique_index]

        if w_new_unique.size > 1:
            y_new_on_old = np.interp(
                w_old,
                w_new_unique,
                y_new_unique,
                left=np.nan,
                right=np.nan,
            )
            delta = y_new_on_old - y_old
            ax_bot.plot(
                w_old,
                delta,
                linestyle='-',
                color='tab:red',
                alpha=0.25,
                linewidth=0.8,
            )

    ax_top.set_title(f'{ext_name}: all orders, old vs new wave grids')
    ax_top.set_ylabel('Value')
    ax_top.grid(True, alpha=0.2)
    ax_top.plot(
        [],
        [],
        linestyle='--',
        color='0.40',
        label='old on old wave grid',
    )
    ax_top.plot(
        [],
        [],
        linestyle='-',
        color='tab:blue',
        label='new on new wave grid',
    )
    ax_top.legend(loc='upper right')

    ax_bot.axhline(0.0, color='k', alpha=0.25, linewidth=0.8)
    ax_bot.set_title('Difference: new(old-grid interp) - old')
    ax_bot.set_xlabel('Wavelength')
    ax_bot.set_ylabel('Delta value')
    ax_bot.grid(True, alpha=0.2)
    fig.tight_layout()


def run_conversion(
    model_in_path: str,
    model_out_path: str,
    wave_new_path: str,
    dry_run: bool,
    make_plot: bool,
    plot_extensions: list[str],
    backup_dir: str | None,
) -> None:
    """Execute sky-model wave-grid conversion for all matching extensions.

    :param model_in_path: input model FITS path
    :param model_out_path: output model FITS path
    :param wave_new_path: target wave-grid FITS path
    :param dry_run: if True, skip backup and write
    :param make_plot: if True, display plots
    :param plot_extensions: extension names to plot
    :param backup_dir: backup destination directory (None -> script dir)
    :return: None
    """
    print(f'Using wave_new:       {wave_new_path}')
    wave_new = load_wave_new(wave_new_path)

    with fits.open(model_in_path) as hdul:
        hdul_out = fits.HDUList([hdu.copy() for hdu in hdul])

    wave_ext = get_wave_extension_index(hdul_out)
    wave_old_raw = np.array(hdul_out[wave_ext].data, dtype=float)

    if wave_old_raw.ndim != 2:
        raise ValueError('WAVE extension must be 2-D.')

    old_dims = wave_old_raw.shape
    n_orders = min(old_dims)
    n_pix_old = max(old_dims)
    wave_old, wave_transposed = orient_orders_pixels(
        wave_old_raw,
        n_orders,
        n_pix_old,
    )

    if wave_new.shape[0] != n_orders:
        emsg = (
            'wave_new order count does not match model WAVE order count: '
            f'{wave_new.shape[0]} != {n_orders}'
        )
        raise ValueError(emsg)

    print(f'model wave_old shape: {wave_old.shape}')
    print(f'wave_new shape:       {wave_new.shape}')

    plot_payload = []
    for ext_index, hdu in enumerate(hdul_out):
        if hdu.data is None:
            continue

        data = np.array(hdu.data)
        if data.ndim != 2:
            continue

        if n_pix_old not in data.shape:
            continue

        data_old_oriented, used_transpose = orient_orders_pixels(
            data,
            n_orders,
            n_pix_old,
        )

        is_discrete = np.issubdtype(data.dtype, np.integer)
        if hdu.name == 'WAVE':
            data_new_oriented = np.array(wave_new, dtype=float)
        else:
            data_new_oriented = interpolate_extension(
                data_old_oriented,
                wave_old,
                wave_new,
                is_discrete,
            )

        if is_discrete:
            data_new_oriented = np.rint(data_new_oriented)
            data_new_oriented = data_new_oriented.astype(data.dtype)

        if used_transpose:
            data_new = data_new_oriented.T
        else:
            data_new = data_new_oriented

        hdul_out[ext_index].data = data_new
        print(
            f'Updated {hdu.name}: {data.shape} -> {data_new.shape} '
            f'(discrete={is_discrete})'
        )

        if make_plot and hdu.name in plot_extensions:
            plot_payload.append(
                (
                    hdu.name,
                    np.array(data_old_oriented, dtype=float),
                    np.array(data_new_oriented, dtype=float),
                )
            )

    if make_plot:
        for ext_name, data_old_plot, data_new_plot in plot_payload:
            plot_extension_comparison(
                ext_name,
                wave_old,
                wave_new,
                data_old_plot,
                data_new_plot,
            )
        if plot_payload:
            plt.show()

    if dry_run:
        print('Dry run enabled: no backup or write performed.')
        return

    if model_in_path == model_out_path and os.path.exists(model_out_path):
        if backup_dir is None:
            backup_dir_use = os.path.dirname(os.path.abspath(__file__))
        else:
            backup_dir_use = backup_dir
        backup_path = backup_file(model_out_path, backup_dir_use)
        print(f'Backup written to: {backup_path}')

    hdul_out.writeto(model_out_path, overwrite=True)
    print(f'Wrote updated model to: {model_out_path}')


def run_self_test() -> None:
    """Run synthetic self-test with small model and wave arrays.

    :return: None
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        model_in_path = os.path.join(tmp_dir, 'sky_model_in.fits')
        model_out_path = os.path.join(tmp_dir, 'sky_model_out.fits')
        wave_new_path = os.path.join(tmp_dir, 'wave_new.fits')

        n_orders = 3
        n_pix_old = 4
        n_pix_new = 7

        wave_old = np.array(
            [
                [100.0, 101.0, 102.0, 103.0],
                [200.0, 201.0, 202.0, 203.0],
                [300.0, 301.0, 302.0, 303.0],
            ]
        )
        wave_new = np.array(
            [
                np.linspace(100.0, 103.0, n_pix_new),
                np.linspace(200.0, 203.0, n_pix_new),
                np.linspace(300.0, 303.0, n_pix_new),
            ]
        )

        sci_sky = wave_old * 0.1
        cal_sky = wave_old * 0.2
        reg_id = np.array(
            [
                [1, 1, 2, 2],
                [2, 2, 3, 3],
                [3, 3, 4, 4],
            ],
            dtype=np.int64,
        )
        weights = np.ones_like(wave_old)
        gradient = np.gradient(sci_sky, axis=1)

        hdul = fits.HDUList()
        hdul.append(fits.PrimaryHDU())
        hdul.append(fits.ImageHDU(data=sci_sky.T, name='SCI_SKY'))
        hdul.append(fits.ImageHDU(data=cal_sky.T, name='CAL_SKY'))
        hdul.append(fits.ImageHDU(data=wave_old.T, name='WAVE'))
        hdul.append(fits.ImageHDU(data=reg_id.T, name='REG_ID'))
        hdul.append(fits.ImageHDU(data=weights.T, name='WEIGHTS'))
        hdul.append(fits.ImageHDU(data=gradient.T, name='GRADIENT'))
        hdul.writeto(model_in_path, overwrite=True)

        fits.writeto(wave_new_path, wave_new, overwrite=True)

        run_conversion(
            model_in_path,
            model_out_path,
            wave_new_path,
            dry_run=False,
            make_plot=False,
            plot_extensions=['SCI_SKY'],
            backup_dir=tmp_dir,
        )

        with fits.open(model_out_path) as hdul_out:
            sci_shape = hdul_out['SCI_SKY'].data.shape
            cal_shape = hdul_out['CAL_SKY'].data.shape
            wave_shape = hdul_out['WAVE'].data.shape
            rid_shape = hdul_out['REG_ID'].data.shape

            expected = (n_pix_new, n_orders)
            if sci_shape != expected:
                raise ValueError(f'SCI_SKY shape {sci_shape} != {expected}')
            if cal_shape != expected:
                raise ValueError(f'CAL_SKY shape {cal_shape} != {expected}')
            if wave_shape != expected:
                raise ValueError(f'WAVE shape {wave_shape} != {expected}')
            if rid_shape != expected:
                raise ValueError(f'REG_ID shape {rid_shape} != {expected}')

            if not np.issubdtype(hdul_out['REG_ID'].data.dtype, np.integer):
                raise ValueError('REG_ID dtype should remain integer.')

        print('Self-test completed successfully.')


def build_argparser() -> argparse.ArgumentParser:
    """Build command-line parser.

    :return: configured parser
    """
    parser = argparse.ArgumentParser(
        description='Spline sky-model FITS extensions onto a new wave grid.'
    )
    parser.add_argument(
        '--model-in',
        default=DEFAULT_MODEL_HA,
        help='Input sky-model FITS file (e.g. HA or HE model).',
    )
    parser.add_argument(
        '--model-out',
        default=None,
        help='Output model FITS. Default: same as --model-in.',
    )
    parser.add_argument(
        '--wave-new',
        default=None,
        help=(
            'Target wave FITS file defining output pixel size. If omitted, '
            'script uses DEFAULT_WAVE_NEW_HA/HE based on --model-in path.'
        ),
    )
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument(
        '--plot-extensions',
        default='SCI_SKY,CAL_SKY',
        help='Comma-separated extension names for comparison plots.',
    )
    parser.add_argument('--self-test', action='store_true')
    return parser


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':
    args = build_argparser().parse_args()

    if args.self_test:
        run_self_test()
    else:
        model_out = args.model_out
        if model_out is None:
            model_out = args.model_in
        wave_new_path = resolve_wave_new_path(args.model_in, args.wave_new)
        run_conversion(
            model_in_path=args.model_in,
            model_out_path=model_out,
            wave_new_path=wave_new_path,
            dry_run=args.dry_run,
            make_plot=args.plot,
            plot_extensions=parse_name_list(args.plot_extensions),
            backup_dir=None,
        )

# =============================================================================
# End of code
# =============================================================================

