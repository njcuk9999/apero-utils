# %%
import eso_programmatic as esop
import pyvo as vo
import tqdm
import yaml

# Constants
ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
TOKEN_AUTHENTICATION_URL = "https://www.eso.org/sso/oidc/token"
AUTH_FILE = "auth.yaml"
RAMPS_DIR = "test-harps"
# RAMPS_DIR = "test-nirps"
READS_DIR = "test-reads"
# TODO: Should this be both raw2raw and raw2master ?
CALIB_TYPE = "raw2raw"
INSTRUMENT = "NIRPS"  # Instrument to get data from
DATESTR = "2022"  # Get data after that date (inclusively)
# less HARPS data and more varied files. Use for testing
# INSTRUMENT = "HARPS"  # Instrument to get data from
# DATESTR = "2021"  # Get data after that date (inclusively)

# %%
with open(AUTH_FILE, "r") as authfile:
    auth_info = yaml.safe_load(authfile)

# %%
token = esop.getToken(auth_info["user"], auth_info["password"])

session = esop.createSession(token)

tap = vo.dal.TAPService(ESO_TAP_OBS, session=session)


# %%
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

# TODO: Choose columns
raw_query = f"""SELECT * FROM dbo.raw
WHERE prog_id='{auth_info['program_id']}'
AND date_obs > '{DATESTR}'
AND instrument='{INSTRUMENT}'"""

# TODO: Choose columns
red_query = f"""SELECT * FROM ivoa.ObsCore
WHERE proposal_id='{auth_info['program_id']}'
AND obs_release_date > '{DATESTR}'
AND instrument_name='{INSTRUMENT}'"""

# TODO: Implement async query option
files_tbl = tap.search(query=raw_query).to_table()

# TODO: Download raw files
# TODO: Handle ramp and reads separately
reads_mask = files_tbl["dp_cat"] == "TECHNICAL"
reads_tbl = files_tbl[reads_mask]
ramps_tbl = files_tbl[~reads_mask]

# %%
# TODO: Handle dir creation
# TODO: Confirm file names
# TODO: Move raw files (sci and calibs) to nightly directories
# TODO: Resolve duplicates
# TODO: Should this loop be only on science files, or CALIB as well?
# NOTE: (Bold in eso calselector ref linked below) Under certain circumstances, including a mixture of science and non-science (e.g. acquisition images) files in the same request results in CalSelector behaving erratically. Users are strongly encouraged to include only science files in their query (see below for the details on how to proceed) and let CalSelector associate the other relevant files.
# WARN: Remove next line after debug
ramps_tbl = ramps_tbl[:3]
for ramp in tqdm.tqdm(ramps_tbl):
    access_url = ramp["access_url"]
    status, filename = esop.downloadURL(access_url, session=session, dirname=RAMPS_DIR)

# %%
# TODO: Do we need to store file assoications somewhere (for APERO, for e.g.)
# TODO: Do we care about complete/certified file and cascade info?
# TODO: Do we want this loop or should all calib files be found with the "main" loop above
# For example, do we overwrite a file if was certified and now is? Do we wait for it to be?
# Do we just flag if downloaded and still not certified, etc.
# Ref: http://archive.eso.org/cms/application_support/calselectorInfo.html
for ramp in ramps_tbl:

    datalink_url = ramp["datalink_url"]
    datalink = vo.dal.adhoc.DatalinkResults.from_result_url(datalink_url, session=session)

    esop.printTableTransposedByTheRecord(datalink.to_table())

    semantics = f"http://archive.eso.org/rdf/datalink/eso#calSelector_{CALIB_TYPE}"

    calib_list_url = next(datalink.bysemantics(semantics)).access_url

    # The calselector URL points to a list of files. Get that list of files
    # Ref: http://archive.eso.org/programmatic/rdf/datalink/eso/
    associated_calib_files = vo.dal.adhoc.DatalinkResults.from_result_url(
        calib_list_url, session=session
    )

    # TODO: Do we do this or we also want sibling_raw?
    # create and use a mask to get only the #calibration entries,
    # given that other entries, like #this or ...#sibiling_raw, could be present:
    calibrator_mask = associated_calib_files['semantics'] == '#calibration'
    calib_urls = associated_calib_files.to_table()[calibrator_mask]['access_url', 'eso_category']

    esop.printTableTransposedByTheRecord(associated_calib_files.to_table())

    # TODO: Display warnings here using caslector info output
    for url, _category in calib_urls:
        status, filename = esop.downloadURL(url, session=session, dirname=RAMPS_DIR)

# %%
# TODO: Do these have associated calibs?
# WARN: Remove next line after debug
reads_tbl = reads_tbl[:2]
for read in reads_tbl:
    access_url = read["access_url"]
    status, filename = esop.downloadURL(access_url, session=session, dirname=READS_DIR)

# %%
# TODO: Download associated calibs
# TODO: Which type of calibs?

# TODO: Handle uncompressing files... Probably best to download in one place, uncompress, and then dispatch? Or uncompress at download time and move the file then

# TODO: Download reduced data
# TODO: Which types of calib? How do we handle last modified date vs original date (re-reductions, etc.)
