"""
file:        eso_programmatic.py
description: ESO python library defining methods to programmatically access the ESO science archive
status:      alpha release
version:     0.3
date:        2021-07-29
type:        not yet a package, just only an evolving library of methods
contact:     https://support.eso.org/
Github repo: https://github.com/almicol/eso_authentication_and_authorisation
LICENSE: GPL3

This is a modified version of the eso_programmatic module.
Some edits are just formatting. Others aim to make the tool a bit more flexible,
for example regarding file path formats.

Original file by @almicol/ESO
Modified by @vandalt (Thomas Vandal)
"""
import cgi
import json
import os
import re

import requests

__name__ = "eso_programmatic"
__version__ = "0.3"

# Definition of ESO services' endpoints
ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
TOKEN_AUTHENTICATION_URL = "https://www.eso.org/sso/oidc/token"


def getToken(username, password):
    """
    Token based authentication to ESO: provide username and password to receive back a JSON Web Token.

    :param username: ESO username
    :type username: str
    :param password: ESO archive password
    :type password: str
    :returns: token
    :rtype: str, optional
    """
    if username is None or password is None:
        return None

    token = None
    try:
        response = requests.get(
            TOKEN_AUTHENTICATION_URL,
            params={
                "response_type": "id_token token",
                "grant_type": "password",
                "client_id": "clientid",
                "username": username,
                "password": password,
            },
        )
        token_response = json.loads(response.content)
        token = token_response["id_token"] + "=="
    except NameError as e:
        print(e)
    except Exception:
        print(
            "*** AUTHENTICATION ERROR: Invalid credentials provided for username %s"
            % (username)
        )

    return token


def createSession(token=None):
    """
    Create requests session, with Auth token if provided

    :param token: ESO token (lasts 8 hrs), defaults to None
    :type token: str, optional
    """
    session = requests.Session()
    if token:
        session.headers["Authorization"] = "Bearer " + token
    return session


def downloadURL(file_url, dirname=".", filename=None, session=None):
    """
    Method to download a file, either anonymously (no session or session not "tokenized"),
    or authenticated (if session with token is provided).

    :param file_url: URL of the file to download on the archive
    :type file_url: str
    :param dirname: Directory where file will be saved, defaults to "."
    :type dirname: str, optional
    :param filename: name of the file on disk, defaults to None
    :type filename: str, optional
    :param session: rerquest session to use for auth, defaults to None
    :type session: requests.Session, optional
    :returns: http status, and filepath on disk (if successful)
    :rtype: Tuple[int, str]
    """

    dirname = dirname or "."
    if not os.access(dirname, os.W_OK):
        raise PermissionError(f"Directory {dirname} is not writable")

    if session is not None:
        response = session.get(file_url, stream=True)
    else:
        # no session -> no authentication
        response = requests.get(file_url, stream=True)

    filename = filename or get_default_filename(file_url, response=response)
    filename = filename.replace(":", "_")

    # define the file path where the file is going to be stored
    filepath = os.path.join(dirname, filename)

    if response.status_code == 200:
        with open(filepath, "wb") as f:
            for chunk in response.iter_content(chunk_size=50000):
                f.write(chunk)

    return (response.status_code, filepath)


def get_default_filename(file_url, response=None):
    # If not provided, define the filename from the response header
    if response is not None:
        # Content disposstion gives attachment; filename={filename} so can extract
        contentdisposition = response.headers.get("Content-Disposition")
        if contentdisposition is not None:
            try:
                value, params = cgi.parse_header(contentdisposition)
                filename = params["filename"]
            except KeyError:
                pass

    # if the response header does not provide a name, derive a name from the URL
    if filename is None:
        # last chance: get anything after the last '/'
        filename = file_url.rsplit("/", maxsplit=1)[-1]

    return filename


# Let's define some methods to nicely print reference files information from the calselector service:

# calselectorInfo(description): [internal] parsing a calselector description


def calselectorInfo(description):
    """Parse the main calSelector description, and fetch: category, complete, certified, mode, and messages."""

    category = ""
    complete = ""
    certified = ""
    mode = ""
    messages = ""

    m = re.search('category="([^"]+)"', description)
    if m:
        category = m.group(1)
    m = re.search('complete="([^"]+)"', description)
    if m:
        complete = m.group(1).lower()
    m = re.search('certified="([^"]+)"', description)
    if m:
        certified = m.group(1).lower()
    m = re.search('mode="([^"]+)"', description)
    if m:
        mode = m.group(1).lower()
    m = re.search('messages="([^"]+)"', description)
    if m:
        messages = m.group(1)

    return category, complete, certified, mode, messages


# - printCalselectorInfo(description, mode_requested): method that returns
#   possible alerts and warnings on the obtained calibration cascade,
#   while printing most relevant info.
def printCalselectorInfo(description, mode_requested):
    """Print the most relevant params contained in the main calselector description."""

    category, complete, certified, mode_executed, messages = calselectorInfo(
        description
    )

    alert = ""
    if complete != "true":
        alert = "ALERT: incomplete calibration cascade"

    mode_warning = ""
    if mode_executed != mode_requested:
        mode_warning = "WARNING: requested mode (%s) could not be executed" % (
            mode_requested
        )

    certified_warning = ""
    if certified != "true":
        certified_warning = 'WARNING: certified="%s"' % (certified)

    print("    calibration info:")
    print("    ------------------------------------")
    print("    science category=%s" % (category))
    print("    cascade complete=%s" % (complete))
    print("    cascade messages=%s" % (messages))
    print("    cascade certified=%s" % (certified))
    print("    cascade executed mode=%s" % (mode_executed))
    print("    full description: %s" % (description))

    return alert, mode_warning, certified_warning


# - printTableTransposedByTheRecord: Utility method to print a table transposed, one record at the time
def printTableTransposedByTheRecord(table):
    """Utility method to print a table transposed, one record at the time"""
    prompt = "    "
    rec_sep = "-" * 105
    print("=" * 115)
    for row in table:
        for col in row.columns:
            print("{0}{1: <14} = {2}".format(prompt, col, row[col]))
        print("{0}{1}".format(prompt, rec_sep))
