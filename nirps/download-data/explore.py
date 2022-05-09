"""
Useful commands to explore the ESO TAP interface
"""
import yaml

# Constants
ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
AUTH_FILE = "auth.yaml"
RAMPS_DIR = "test-harps"
# RAMPS_DIR = "test-nirps"
READS_DIR = "test-reads"
CALIB_TYPE = "raw2raw"
INSTRUMENT = "NIRPS"  # Instrument to get data from
DATESTR = "2022"  # Get data after that date (inclusively)
# less HARPS data and more varied files. Use for testing
# INSTRUMENT = "HARPS"  # Instrument to get data from
# DATESTR = "2021"  # Get data after that date (inclusively)

# %%
with open(AUTH_FILE, "r") as authfile:
    auth_info = yaml.safe_load(authfile)


list_tables_query = """SELECT table_name, description
FROM TAP_SCHEMA.tables
ORDER BY table_name"""

raw_cols_query = """SELECT column_name, datatype, arraysize, xtype, unit, ucd, description
FROM TAP_SCHEMA.columns
WHERE table_name='dbo.raw'
"""

red_cols_query = """SELECT column_name, datatype, arraysize, xtype, unit, ucd, description
FROM TAP_SCHEMA.columns
WHERE table_name='ivoa.ObsCore'
"""

raw_query = f"""SELECT * FROM dbo.raw
WHERE prog_id='{auth_info['program_id']}'
AND date_obs > '{DATESTR}'
AND instrument='{INSTRUMENT}'"""

# %%
raw_query = f"""SELECT * FROM dbo.raw
WHERE prog_id='{auth_info['program_id']}'
AND date_obs >= '2022-05-05'
AND date_obs <= '2022-05-05T23:59:59'"""

# %%
red_query = f"""SELECT * FROM ivoa.ObsCore
WHERE proposal_id='{auth_info['program_id']}'
AND obs_release_date > '{DATESTR}'
AND instrument_name='{INSTRUMENT}'"""
