import argparse
import os
from typing import Any

from astropy.io import fits


def update_header(filename: str, key: str, value: Any, ext: int = 0):
    """
    Update the header of {filename} with {key} and {value} for extension {ext}
    :param filename:
    :param key:
    :param value:
    :param ext:
    :return:
    """
    # deal with file not on disk
    if not os.path.exists(filename):
        print('File does not exist: {}'.format(filename))
        return
    # Open the file in update mode
    with fits.open(filename, mode='update') as hdul:
        try:
            # Access the extension you want
            hdr = hdul[ext].header
            # Deal with bad key
            if key not in hdr:
                uinput = input('Key {} not found in header, add? [Y]es or [N]o:'
                               '\t'.format(key))
                if not ('y' in uinput.lower()):
                    print('Key "{}" not found in header of {}'.format(key, filename))
                    return
            else:
                # Modify the header key
                uargs = [key, hdr[key], value]
                uinput = input('Update key {} = {} --> {}? '
                               '[Y]es or [N]o:\t'.format(*uargs))

            # try to evaluate values
            try:
                value = eval(value)
            except:
                value = str(value).strip()

            if 'y' in uinput.lower():
                hdr[key] = value
            else:
                return
        except Exception as e:
            eargs = [filename, type(e), str(e)]
            print('Error updating file: {}\n\t{}: {}'.format(*eargs))


if __name__ == "__main__":
    # Get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str,
                        help='Fits file to modify')
    parser.add_argument('--key', type=str, required=True,
                        help='Header key to change')
    parser.add_argument('--value', type=str, required=True,
                        help='Header value to change')
    parser.add_argument('--ext', type=int, default=0,
                        help='Extension header is in')
    # parse arguments
    args = parser.parse_args()
    # update header (open + update)
    update_header(args.filename, args.key, args.value, args.ext)




