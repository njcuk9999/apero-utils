import yaml
import eso_programmatic as esop
import pyvo as vo


# Constants
ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
TOKEN_AUTHENTICATION_URL = "https://www.eso.org/sso/oidc/token"
AUTH_FILE = "auth.yaml"

with open(AUTH_FILE, 'r') as authfile:
    auth_info = yaml.safe_load(authfile)

token = esop.getToken(auth_info["username"], auth_info["password"])

# %%
session = esop.createSession(token)


tap = vo.dal.TAPService(ESO_TAP_OBS, session=session)
