# %%
import yaml
import eso_programmatic as esop
import pyvo as vo


# Constants
ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
TOKEN_AUTHENTICATION_URL = "https://www.eso.org/sso/oidc/token"
AUTH_FILE = "auth.yaml"
INSTRUMENT = "NIRPS"
DATESTR = "2022%"

# %%
with open(AUTH_FILE, 'r') as authfile:
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
WHERE proposal_id='{auth_info['program_id']}'
AND date_obs LIKE '{DATESTR}'
AND instrument='{INSTRUMENT}'"""

# TODO: Choose columns
red_query = f"""SELECT * FROM ivoa.ObsCore
WHERE proposal_id='{auth_info['program_id']}'
AND date_obs LIKE '{DATESTR}'
AND instrument='{INSTRUMENT}'"""
