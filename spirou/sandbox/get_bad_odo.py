import requests
from astropy.table import Table
# Bad odometers

URL_BASE = ('https://docs.google.com/spreadsheets/d/'
            '{}/gviz/tq?tqx=out:csv&sheet={}')
SHEET_ID = '1gvMp1nHmEcKCUpxsTxkx-5m115mLuQIGHhxJCyVoZCM'
WORKSHEET = 0
BAD_ODO_URL = URL_BASE.format(SHEET_ID, WORKSHEET)
def get_bad_odo(url):
    # fetch data
    data = requests.get(url)
    tbl = Table.read(data.text, format='ascii')
    # Convert types
    tbl['ODOMETER'] = tbl['ODOMETER'].astype(str)
    tbl['PP'] = tbl['PP'] == 'TRUE'
    tbl['RV'] = tbl['RV'] == 'TRUE'
    return tbl['ODOMETER']

# Remove bad odometers
tbl_bad_odo = get_bad_odo(BAD_ODO_URL)
