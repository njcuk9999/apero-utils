#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Spline OH PCA vectors from an old wave grid onto a new wave grid.

This script reads an old wavelength solution, a new wavelength solution,
and an old OH PCA matrix. The old PCA vectors are reshaped to match the
old wave grid, spline-interpolated order-by-order onto the new wave grid,
flattened back to a 2-D PCA matrix, and written to an output FITS file.
"""

import argparse
import datetime
import os
import shutil
import tempfile

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import interpolate

# =============================================================================
# Define variables
# =============================================================================
DEFAULT_WAVE_OLD = (
    '/scratch2/spirou/drs-data/spirou_xxs_07/calib/'
    '2F3798BAE7a_pp_e2dsff_AB_wavesol_ref_AB.fits'
)
DEFAULT_WAVE_NEW = (
    '/scratch2/spirou/drs-data/spirou_xxs_08/calib/REF/'
    '2F3798BAE7a_pp_e2dsff_AB_wavesol_ref_AB.fits'
)
DEFAULT_PCA_OLD = (
    '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/apero/'
    'apero-assets/spirou/telluric/sky_PCs.fits'
)
DEFAULT_PCA_NEW = (
    '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/apero/'
    'apero-assets/spirou/telluric/sky_PCs.fits'
)


# =============================================================================
# Define functions
# =============================================================================
def first_data_hdu_index(hdul: fits.HDUList) -> int:
    """Return index of first HDU with image data.

    :param hdul: opened FITS HDU list
    :return: index of first HDU that contains a numpy array
    :raises ValueError: if no HDU has image data
    """
    for index, hdu in enumerate(hdul):
        if hdu.data is not None:
            return index
    raise ValueError('No image data found in FITS file.')


def load_fits_array(path: str) -> np.ndarray:
    """Read and return the first image-data array from a FITS file.

    :param path: absolute or relative FITS file path
    :return: numpy array copy of the first image HDU data
    :raises ValueError: if file does not contain image data
    """
    with fits.open(path) as hdul:
        hdu_index = first_data_hdu_index(hdul)
        return np.array(hdul[hdu_index].data)


def backup_file(source_path: str, backup_dir: str) -> str:
    """Create a timestamped backup of a file in ``backup_dir``.

    :param source_path: file to back up
    :param backup_dir: directory where backup is written
    :return: full path to the created backup file
    """
    os.makedirs(backup_dir, exist_ok=True)
    stamp = datetime.datetime.now(datetime.UTC).strftime('%Y%m%dT%H%M%SZ')
    base_name = os.path.basename(source_path)
    backup_name = f'{base_name}.backup_{stamp}'
    backup_path = os.path.join(backup_dir, backup_name)
    shutil.copy2(source_path, backup_path)
    return backup_path


def interpolate_one_order(
    wave_old_order: np.ndarray,
    pca_old_order: np.ndarray,
    wave_new_order: np.ndarray,
) -> np.ndarray:
    """Spline-interpolate one PCA order from old to new wavelength points.

    :param wave_old_order: old wavelength vector for one order
    :param pca_old_order: old PCA values for one order
    :param wave_new_order: target wavelength vector for one order
    :return: interpolated PCA values on ``wave_new_order``
    """
    # Keep only finite sample pairs before fitting a spline.
    finite = np.isfinite(wave_old_order) & np.isfinite(pca_old_order)
    x_valid = wave_old_order[finite]
    y_valid = pca_old_order[finite]

    # If too little data is available, return NaNs for this full order.
    if x_valid.size < 2:
        return np.full_like(wave_new_order, np.nan, dtype=float)

    # Enforce monotonic x for spline fitting.
    sort_index = np.argsort(x_valid)
    x_sorted = x_valid[sort_index]
    y_sorted = y_valid[sort_index]

    # Remove duplicate wavelengths to avoid singular spline fits.
    x_unique, unique_index = np.unique(x_sorted, return_index=True)
    y_unique = y_sorted[unique_index]

    # With <4 points, lower the spline degree to keep the fit valid.
    if x_unique.size < 2:
        return np.full_like(wave_new_order, np.nan, dtype=float)
    degree = min(3, x_unique.size - 1)

    # Use boundary-value extrapolation if new grid steps outside old range.
    spline = interpolate.InterpolatedUnivariateSpline(
        x_unique,
        y_unique,
        k=degree,
        ext=3,
    )
    return spline(wave_new_order)


def spline_pca(
    wave_old: np.ndarray,
    wave_new: np.ndarray,
    pca_old: np.ndarray,
) -> np.ndarray:
    """Interpolate old flattened PCA vectors from ``wave_old`` to ``wave_new``.

    Accepts both PCA layouts: ``(n_comp, n_orders * n_pix_old)`` and
    ``(n_orders * n_pix_old, n_comp)``. Output keeps the same orientation as
    the input PCA matrix.

    Special case: for layout ``(n_orders * n_pix_old, n_col)``, if
    ``n_col > 1`` the first column is treated as the old wave grid for
    interpolation and the output first column is set directly to
    ``wave_new.flatten()``.

    :param wave_old: old wavelength grid, shape ``(n_orders, n_pix_old)``
    :param wave_new: new wavelength grid, shape ``(n_orders, n_pix_new)``
    :param pca_old: old PCA matrix in one of the accepted layouts
    :return: new PCA matrix matching input orientation
    :raises ValueError: for incompatible wave/PCA dimensions
    """
    if wave_old.ndim != 2:
        raise ValueError('wave_old must be 2-D.')
    if wave_new.ndim != 2:
        raise ValueError('wave_new must be 2-D.')
    if pca_old.ndim != 2:
        raise ValueError('pca_old must be 2-D.')

    n_orders_old, n_pix_old = wave_old.shape
    n_orders_new, n_pix_new = wave_new.shape

    if n_orders_old != n_orders_new:
        emsg = (
            'wave_old and wave_new must have the same number of orders: '
            f'{n_orders_old} != {n_orders_new}'
        )
        raise ValueError(emsg)

    expected_size = n_orders_old * n_pix_old

    # Handle table-like format where rows are flattened pixels and the
    # first column carries the old wavelength grid.
    if pca_old.shape[0] == expected_size and pca_old.shape[1] > 1:
        n_columns = pca_old.shape[1]
        wave_old_from_pca = pca_old[:, 0].reshape(n_orders_old, n_pix_old)
        out = np.full((n_orders_new * n_pix_new, n_columns), np.nan)
        out[:, 0] = wave_new.reshape(-1)

        for column_num in range(1, n_columns):
            old_col = pca_old[:, column_num].reshape(n_orders_old, n_pix_old)
            new_col = np.full((n_orders_new, n_pix_new), np.nan)
            for order_index in range(n_orders_old):
                new_col[order_index] = interpolate_one_order(
                    wave_old_from_pca[order_index],
                    old_col[order_index],
                    wave_new[order_index],
                )
            out[:, column_num] = new_col.reshape(-1)
        return out

    if pca_old.shape[1] == expected_size:
        use_transpose = False
        pca_old_use = pca_old
    elif pca_old.shape[0] == expected_size:
        use_transpose = True
        pca_old_use = pca_old.T
    else:
        emsg = (
            'pca_old must have one axis equal to wave_old flattened size: '
            f'{pca_old.shape} with expected flattened size {expected_size}'
        )
        raise ValueError(emsg)

    n_components = pca_old_use.shape[0]
    pca_old_cube = pca_old_use.reshape(n_components, n_orders_old, n_pix_old)
    pca_new_cube = np.full(
        (n_components, n_orders_new, n_pix_new),
        np.nan,
        dtype=float,
    )

    # Interpolate each component/order independently to preserve structure.
    for comp_index in range(n_components):
        for order_index in range(n_orders_old):
            pca_new_cube[comp_index, order_index] = interpolate_one_order(
                wave_old[order_index],
                pca_old_cube[comp_index, order_index],
                wave_new[order_index],
            )

    pca_new_comp_first = pca_new_cube.reshape(
        n_components,
        n_orders_new * n_pix_new,
    )
    if use_transpose:
        return pca_new_comp_first.T
    return pca_new_comp_first


def write_pca_like_template(
    output_path: str,
    template_path: str,
    new_pca: np.ndarray,
) -> None:
    """Write ``new_pca`` into first image HDU of a template FITS file.

    :param output_path: destination FITS file path
    :param template_path: FITS file used as header/extension template
    :param new_pca: PCA matrix to write
    :return: None
    """
    with fits.open(template_path) as hdul:
        hdul_out = fits.HDUList([hdu.copy() for hdu in hdul])

    hdu_index = first_data_hdu_index(hdul_out)
    hdul_out[hdu_index].data = np.array(new_pca)
    hdul_out.writeto(output_path, overwrite=True)


def parse_index_list(value: str) -> list[int]:
    """Parse a comma-separated index list.

    :param value: string like ``'0,1,2'``
    :return: list of non-negative integer indices
    :raises ValueError: if parsing fails or an index is negative
    """
    tokens = [token.strip() for token in value.split(',') if token.strip()]
    if not tokens:
        return []

    out = []
    for token in tokens:
        index = int(token)
        if index < 0:
            raise ValueError('Indices must be non-negative.')
        out.append(index)
    return out


def plot_pca_comparison(
    wave_old: np.ndarray,
    wave_new: np.ndarray,
    pca_old: np.ndarray,
    pca_new: np.ndarray,
    component_indices: list[int],
) -> None:
    """Plot old/new PCA curves and their difference for components.

    Accepts both PCA layouts: ``(n_comp, wave_size)`` and
    ``(wave_size, n_comp)``.

    :param wave_old: old wavelength grid, shape ``(n_orders, n_pix_old)``
    :param wave_new: new wavelength grid, shape ``(n_orders, n_pix_new)``
    :param pca_old: old flattened PCA
    :param pca_new: new flattened PCA
    :param component_indices: component IDs to show
    :return: None
    """
    n_orders_old, n_pix_old = wave_old.shape
    n_orders_new, n_pix_new = wave_new.shape
    if n_orders_old != n_orders_new:
        raise ValueError('Cannot plot: order counts differ between grids.')

    old_size = n_orders_old * n_pix_old
    new_size = n_orders_new * n_pix_new
    has_wave_column = False

    if pca_old.shape[1] == old_size:
        pca_old_use = pca_old
    elif pca_old.shape[0] == old_size:
        pca_old_use = pca_old.T
        # Table-format input: row 0 is wavelength column, not a PCA vector.
        has_wave_column = True
    else:
        emsg = 'Cannot plot: pca_old shape incompatible with wave_old.'
        raise ValueError(emsg)

    if pca_new.shape[1] == new_size:
        pca_new_use = pca_new
    elif pca_new.shape[0] == new_size:
        pca_new_use = pca_new.T
    else:
        emsg = 'Cannot plot: pca_new shape incompatible with wave_new.'
        raise ValueError(emsg)

    n_rows = pca_old_use.shape[0]
    old_cube = pca_old_use.reshape(n_rows, n_orders_old, n_pix_old)
    new_cube = pca_new_use.reshape(n_rows, n_orders_new, n_pix_new)

    orders = range(n_orders_old)

    for comp in component_indices:
        if has_wave_column:
            max_comp = n_rows - 2
            if comp > max_comp:
                emsg = (
                    f'Skipping component {comp}: '
                    f'max PCA component index is {max_comp}'
                )
                print(emsg)
                continue
            data_row = comp + 1
        else:
            max_comp = n_rows - 1
            if comp > max_comp:
                emsg = (
                    f'Skipping component {comp}: '
                    f'max PCA component index is {max_comp}'
                )
                print(emsg)
                continue
            data_row = comp

        if comp < 0:
            emsg = f'Skipping component {comp}: index must be >= 0.'
            print(emsg)
            continue

        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=False)
        ax_top = axes[0]
        ax_bot = axes[1]

        for order_num in orders:
            wave_old_ord = wave_old[order_num]
            wave_new_ord = wave_new[order_num]
            pca_old_ord = old_cube[data_row, order_num]
            pca_new_ord = new_cube[data_row, order_num]

            ax_top.plot(
                wave_old_ord,
                pca_old_ord,
                linestyle='--',
                color='0.50',
                alpha=0.35,
                linewidth=0.8,
            )
            ax_top.plot(
                wave_new_ord,
                pca_new_ord,
                linestyle='-',
                color='tab:blue',
                alpha=0.25,
                linewidth=0.8,
            )

            # Compare on old grid by interpolating the new curve to old x.
            order_sort = np.argsort(wave_new_ord)
            wave_new_sort = wave_new_ord[order_sort]
            pca_new_sort = pca_new_ord[order_sort]
            wave_new_unique, unique_idx = np.unique(
                wave_new_sort,
                return_index=True,
            )
            pca_new_unique = pca_new_sort[unique_idx]
            if wave_new_unique.size > 1:
                pca_new_on_old = np.interp(
                    wave_old_ord,
                    wave_new_unique,
                    pca_new_unique,
                    left=np.nan,
                    right=np.nan,
                )
                pca_diff = pca_new_on_old - pca_old_ord
                ax_bot.plot(
                    wave_old_ord,
                    pca_diff,
                    linestyle='-',
                    color='tab:red',
                    alpha=0.25,
                    linewidth=0.8,
                )

        ax_top.set_title(
            f'PCA component {comp}: all orders, old vs new wavelength grids'
        )
        ax_top.set_ylabel('PCA value')
        ax_top.grid(True, alpha=0.2)
        ax_top.plot(
            [],
            [],
            linestyle='--',
            color='0.40',
            label='old PCA on old wave grid',
        )
        ax_top.plot(
            [],
            [],
            linestyle='-',
            color='tab:blue',
            label='new PCA on new wave grid',
        )
        ax_top.legend(loc='upper right')

        ax_bot.axhline(0.0, color='k', alpha=0.25, linewidth=0.8)
        ax_bot.set_title('Difference: new(old-grid interp) - old')
        ax_bot.set_xlabel('Wavelength')
        ax_bot.set_ylabel('Delta PCA')
        ax_bot.grid(True, alpha=0.2)
        fig.tight_layout()

    if component_indices:
        plt.show()


def run_conversion(
    wave_old_path: str,
    wave_new_path: str,
    pca_old_path: str,
    pca_new_path: str,
    dry_run: bool,
    make_plot: bool,
    plot_components: list[int],
) -> None:
    """Execute the full wave-grid PCA conversion pipeline.

    :param wave_old_path: old wave FITS path
    :param wave_new_path: new wave FITS path
    :param pca_old_path: source PCA FITS path
    :param pca_new_path: destination PCA FITS path
    :param dry_run: if True, do not backup or write files
    :param make_plot: if True, show PCA-vs-wavelength comparison plots
    :param plot_components: PCA component IDs to plot
    :return: None
    """
    wave_old = load_fits_array(wave_old_path)
    wave_new = load_fits_array(wave_new_path)
    pca_old = load_fits_array(pca_old_path)

    print(f'wave_old shape: {wave_old.shape}')
    print(f'wave_new shape: {wave_new.shape}')
    print(f'pca_old shape:  {pca_old.shape}')

    pca_new = spline_pca(wave_old, wave_new, pca_old)
    print(f'pca_new shape:  {pca_new.shape}')

    if make_plot:
        plot_pca_comparison(
            wave_old,
            wave_new,
            pca_old,
            pca_new,
            plot_components,
        )

    if dry_run:
        print('Dry run enabled: no backup or output write performed.')
        return

    script_dir = os.path.dirname(os.path.abspath(__file__))
    if os.path.exists(pca_new_path):
        backup_path = backup_file(pca_new_path, script_dir)
        print(f'Backup written to: {backup_path}')

    # Prefer preserving output-file metadata when output already exists.
    if os.path.exists(pca_new_path):
        template_path = pca_new_path
    else:
        template_path = pca_old_path
    write_pca_like_template(pca_new_path, template_path, pca_new)
    print(f'Wrote interpolated PCA to: {pca_new_path}')


def run_self_test() -> None:
    """Run a tiny synthetic test to validate shape handling and interpolation.

    :return: None
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        wave_old_path = os.path.join(tmp_dir, 'wave_old.fits')
        wave_new_path = os.path.join(tmp_dir, 'wave_new.fits')
        pca_old_path = os.path.join(tmp_dir, 'pca_old.fits')
        pca_new_path = os.path.join(tmp_dir, 'pca_new.fits')

        wave_old = np.array([[1.0, 2.0, 3.0], [10.0, 11.0, 12.0]])
        wave_new = np.array([[1.0, 1.5, 2.5, 3.0], [10.0, 10.5, 11.5, 12.0]])

        comp0 = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        comp1 = np.array([[2.0, 4.0, 6.0], [8.0, 10.0, 12.0]])
        pca_old = np.vstack([comp0.ravel(), comp1.ravel()])

        fits.writeto(wave_old_path, wave_old, overwrite=True)
        fits.writeto(wave_new_path, wave_new, overwrite=True)
        fits.writeto(pca_old_path, pca_old, overwrite=True)
        fits.writeto(pca_new_path, pca_old, overwrite=True)

        wave_old_read = load_fits_array(wave_old_path)
        wave_new_read = load_fits_array(wave_new_path)
        pca_old_read = load_fits_array(pca_old_path)
        pca_new = spline_pca(wave_old_read, wave_new_read, pca_old_read)
        write_pca_like_template(pca_new_path, pca_new_path, pca_new)

        out = load_fits_array(pca_new_path)
        expected_shape = (2, 2 * 4)
        if out.shape != expected_shape:
            emsg = (
                'Self-test failed: output shape '
                f'{out.shape} != {expected_shape}'
            )
            raise ValueError(emsg)
        print('Self-test completed successfully.')


def build_argparser() -> argparse.ArgumentParser:
    """Create and return command-line parser.

    :return: configured ``argparse.ArgumentParser`` instance
    """
    parser = argparse.ArgumentParser(
        description='Spline OH PCA vectors from old wave grid to new wave grid.'
    )
    parser.add_argument('--wave-old', default=DEFAULT_WAVE_OLD)
    parser.add_argument('--wave-new', default=DEFAULT_WAVE_NEW)
    parser.add_argument('--pca-old', default=DEFAULT_PCA_OLD)
    parser.add_argument('--pca-new', default=DEFAULT_PCA_NEW)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument(
        '--plot-components',
        default='0,1,2',
        help='Comma-separated PCA component indices (zero-based).',
    )
    return parser


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':
    args = build_argparser().parse_args()

    if args.self_test:
        run_self_test()
    else:
        plot_components = parse_index_list(args.plot_components)
        run_conversion(
            args.wave_old,
            args.wave_new,
            args.pca_old,
            args.pca_new,
            args.dry_run,
            args.plot,
            plot_components,
        )

# =============================================================================
# End of code
# =============================================================================
