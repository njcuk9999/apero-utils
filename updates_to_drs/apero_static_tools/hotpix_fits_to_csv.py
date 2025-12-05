from astropy.table import Table
import os

PATH = '/data/spirou/drs-data/July2025/assets/engineering/'
HOTPIX_FILE = 'static_hotpix_spirou.fits'
outfile = HOTPIX_FILE.replace('.fits', '.csv')

data = Table.read(os.path.join(PATH, HOTPIX_FILE))

data.write(os.path.join(PATH, outfile), format='csv', overwrite=True)

