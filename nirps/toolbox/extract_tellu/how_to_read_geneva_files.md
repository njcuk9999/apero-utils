# How to read Geneva drs files

## Downloading files

There is 2 ways to download the Geneva drs files. 

The first way is to copy them directly from `/cosmos99/nirps/geneva-data/DRS-3.0.0/reduced/...` as you would do with APERO data (See [How do I get data](https://github.com/njcuk9999/apero-drs/wiki/How-do-I-get-data)).

The second way is to download them from DACE using the python API. The full documentation of this method can be found [here](https://dace-query.readthedocs.io/en/latest/).
If it is your first time using the `dace-query` package, you will need to do the Authentication step to have access to NIRPS data.

## Opening fits files

After installing the astropy package in your python environment, start your file with
`from astropy.io import fits`. To open the file in Python, simply write `hdul = fits.open(FILENAME)` with FILENAME being the name of the fits file.

## s1d file structure
`hdul[0].header` is the main header.

`hdul[1].data` is the data of the file. This data is stored in a `FITS_rec` object. `FITS_rec` is a subclass of [numpy's record arrays](https://numpy.org/doc/stable/user/basics.rec.html#record-arrays).
So every method that works for numpy record array will work on the data.

The Geneva s1d files contains 5 attributes:
- (vacuum) wavelength
- wavelength_air
- flux
- error
- quality

All of these can be access by using the dot operator (e.g. `hdul[1].data.wavelength` for the wavelength array).

## s2d file structure

`hdul[0].header` is the main header.

The data is contained in 7 matrices. They are the following:
- `hdul[1].data` is the science data (SCIDATA)
- `hdul[2].data` is the error data (ERRDATA)
- `hdul[3].data` is the quality data (QUALDATA)
- `hdul[4].data` is the vacuum wavelength (WAVEDATA_VAC_BARY)
- `hdul[5].data` is the air wavelength (WAVEDATA_AIR_BARY)
- `hdul[6].data` has extension name DLLDATA_VAC_BARY, not sure what it means
- `hdul[7].data` has extension name DLLDATA_AIR_BARY, not sure what it means

These matrices are of dimension 71 (number of orders) by 4084 (number of pixels per order).

## Differences with APERO

In the s2d files, the Geneva drs has less orders than APERO. Some orders between photometric bands, where the tellurics are the most present, are not in the Geneva s2d files.
This creates a difference in the numbering of the orders between Geneva and APERO. I haven't found a robust way to equate APERO order to Geneva order yet.
If someone knows a robust way to do this, please write it in this guide.

The wavelengths in all Geneva drs files are already BERV corrected.
`hdul[0].header['HIERARCH ESO QC BERV']` contains the BERV, in km/s, used in the correction.