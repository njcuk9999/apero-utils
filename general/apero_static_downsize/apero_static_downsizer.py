import glob
from astropy.io import fits
import numpy as np
import os


IMAGE_HDU_TYPES = (fits.ImageHDU, fits.CompImageHDU)

# =============================================================================
# Define variables
# =============================================================================
# Define path to downsize files
PATH = '/scratch2/spip/misc/spip_March2026/neil_engineering/'
# define the files to try and downsize
FILE_PATTERN = 'static_*.fits'
# define the minimum file size in MB to consider for downcasting
MIN_SIZE_MB = 50


# =============================================================================
# Define functions
# =============================================================================
def batch_downcast(path: str, pattern: str = "*.fits",
                   min_size_mb: float = 50):
    """
    Batch downcast image FITS HDUs in a directory to 32-bit float if they are
    above a certain size threshold, while leaving table HDUs unchanged.

    :param path: string,
    :param pattern:
    :param min_size_mb:
    :return:
    """
    # Convert MB to Bytes
    min_size_bytes = min_size_mb * 1024 * 1024

    # Find all fits files
    files = glob.glob(os.path.join(path, pattern))
    # Print the number of files found
    print(f"Found {len(files)} files matching pattern '{pattern}' in '{path}'")

    for filename in files:
        # Skip gzipped files or directories
        if filename.endswith('.gz') or not os.path.isfile(filename):
            continue

        # Check file size
        file_size = os.path.getsize(filename)
        if file_size <= min_size_bytes:
            # Optional: print(f"Skipping {filename} (too small: {file_size / (1024**2):.1f} MB)")
            continue

        print(f"Processing: {filename} ({file_size / (1024 ** 2):.1f} MB)")
        try:
            with fits.open(filename) as hdul:
                modified = False
                for hdu in hdul:
                    if isinstance(hdu, IMAGE_HDU_TYPES) and hdu.data is not None:
                        # Convert to float32 (BITPIX -32)
                        hdu.data = hdu.data.astype(np.float32)
                        modified = True

                if modified:
                    output_name = filename.replace(".fits", "_32bit.fits")
                    hdul.writeto(output_name, overwrite=True)
                    print(f"  -> Saved as {output_name}")

        except Exception as e:
            print(f"  !! Failed to process {filename}: {e}")


if __name__ == "__main__":
    # Run batch downcasting
    batch_downcast(PATH, FILE_PATTERN, MIN_SIZE_MB)