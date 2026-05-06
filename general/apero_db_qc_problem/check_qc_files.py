from astropy.io import fits
from astropy.table import Table
import os


PID = 'PID-00017429490357374720-S4HF'

DB_KIND = 'mysql+pymysql'
DB_HOST = 'cosmos.astro.umontreal.ca'
DB_USER = 'spirou'
DB_DATABASE = 'spirou'

FINDEX_TABLENAME = 'findex_spirou_offline_db'
LOG_TABLENAME = 'log_spirou_offline_db'

QUERY = """
SELECT BLOCK_KIND, OBS_DIR, FILENAME, PASSED_ALL_QC 
FROM {FINDEX_TABLENAME} fdb
LEFT JOIN (
    SELECT PID, MAX(PASSED_ALL_QC) AS PASSED_ALL_QC
    FROM {LOG_TABLENAME}
    GROUP BY PID
) ldb
    ON fdb.KW_PID = ldb.PID
WHERE fdb.KW_PID = '{PID}'
"""

BLOCK_KIND: dict = {
    'raw':   '/cosmos99/spirou/apero-data/spirou_offline/raw/',
    'tmp':   '/cosmos99/spirou/apero-data/spirou_offline/tmp/',
    'calib': '/cosmos99/spirou/apero-data/spirou_offline/calib/',
    'red':   '/cosmos99/spirou/apero-data/spirou_offline/red/',
    'tellu': '/cosmos99/spirou/apero-data/spirou_offline/tellu/',
    'out':   '/cosmos99/spirou/apero-data/spirou_offline/out/',
    'lbl':   '/cosmos99/spirou/apero-data/spirou_offline/lbl/',
}


def get_files_from_db(password):

    from sqlalchemy import create_engine
    import pandas as pd

    # Create a database engine (replace with your actual database URI)
    engine = create_engine(f'{DB_KIND}://{DB_USER}:'
                           f'{password}@{DB_HOST}/{DB_DATABASE}')

    # Execute the query and fetch results into a DataFrame
    df = pd.read_sql_query(QUERY.format(FINDEX_TABLENAME=FINDEX_TABLENAME,
                                        LOG_TABLENAME=LOG_TABLENAME,
                                        PID=PID), engine)

    # construct filenames
    filenames = []
    db_qc = []

    for row in range(len(df)):

        block_kind = BLOCK_KIND[df['BLOCK_KIND'].iloc[row]]

        filenames.append(os.path.join(block_kind, df['OBS_DIR'].iloc[row],
                                      df['FILENAME'].iloc[row]))
        db_qc.append(df['PASSED_ALL_QC'].iloc[row])

    return filenames, db_qc


def check_header(filename):
    hdr = fits.getheader(filename)
    if 'QCC_ALL' in hdr:
        print('QCC_ALL', '=', hdr['QCC_ALL'])
    else:
        print('No QCC_ALL in header')


def check_param_table(filename):
    table = Table.read(filename, hdu='PARAM_TABLE')
    for i in range(len(table)):
        if 'PASSED_ALL_QC' in table[i]['NAME']:
            print(table[i]['NAME'], '=', table[i]['VALUE'])



if __name__ == '__main__':

    # ask for db password
    password = input(f'Enter db password for {DB_USER}@{DB_HOST}:\t')
    # get all files with the same PID
    files, db_qcs = get_files_from_db(password)
    # loop around files
    for filename in files:
        print('\n')
        print('='*50)
        print(filename)
        print('='*50)

        print('\nFrom DB:')
        print('PASSED_ALL_QC =', db_qcs[files.index(filename)])

        print('\nFrom header:')
        check_header(filename)

        print('\nFrom PARAM_TABLE:')
        check_param_table(filename)

