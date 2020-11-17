"""
Standalone script to load all bad odometers in an astropy table.

This can also be used for any other google sheet by changing the id and
the tab.
"""
import requests
from astropy.table import Table

URL_BASE = ('https://docs.google.com/spreadsheets/d/'
            '{}/gviz/tq?tqx=out:csv&sheet={}')
SHEET_ID = '1gvMp1nHmEcKCUpxsTxkx-5m115mLuQIGHhxJCyVoZCM'
WORKSHEET = 0

# fetch data
url = URL_BASE.format(SHEET_ID, WORKSHEET)
data = requests.get(url)
tbl = Table.read(data.text, format='ascii')

# Convert boolean
tbl['PP'] = tbl['PP'] == 'TRUE'
tbl['RV'] = tbl['RV'] == 'TRUE'
