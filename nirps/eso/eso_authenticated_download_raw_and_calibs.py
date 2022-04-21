#!/usr/bin/env python
# ------------------------------------------------------------------------------------
# Name: eso_authenticated_download_raw_and_calibs.py
# Version: 2020-09-04
# Author: A.Micol, Archive Science Group, ESO
# Purpose: Python 3 example on how to authenticate and download raw science frames
#          and, for each of them, the associated calibration files.
# !!!!!!!!
# !CAVEAT! This is just an example showing the business logic.
# !CAVEAT! There is basically:
# !CAVEAT!  - no error handling in place.
# !CAVEAT!  - no optimisation regarding files downloaded
# !CAVEAT!    (if N different raw science frames require the same calibration files,
# !CAVEAT!     these will be downloaded N times)
# !CAVEAT!  - no check if enough disk space is available for download.
# !CAVEAT!  - etc.
# !CAVEAT! Use at your own risk.
# !!!!!!!!
#
# Business logic:
# - AUTHENTICATE WITH ESO, SO THAT I CAN DOWNLOAD BOTH THE PUBLIC FILES AND THE PROPRIETARY FILES YOU HAVE ACCESS TO
# - FIND SOME SCIENCE RAW FILES
# - FOR EACH SCIENCE RAW FILE:
#     1.- RETRIEVE IT
#     2.- GET ITS LIST OF CALIBRATION FILES (without saving it into a file)
#     3.- PRINT CASCADE INFORMATION AND MAIN DESCRIPTION
#     4.- DOWNLOAD EACH #calibration FILE FROM THE LIST
#     5.- PRINT ANY ALERT OR WARNING ENCOUNTERED BY THE PROCESS THAT GENERATES THE CALIBRATION CASCADE
#
# Limitation: This script does not download siblings, nor the association trees.
#
# Documentation: http://archive.eso.org/cms/application_support/calselectorInfo.html
#
# Contact: In case of questions, please send an email to: usd-help@eso.org
#          with the following subject: programmatic access (eso_authenticated_download_raw_and_calibs.py)
# -------
# Executing the script you might get a WARNING from astropy:
# WARNING: W35: None:6:6: W35: 'value' attribute required for INFO elements [astropy.io.votable.tree]
# That is due to a astropy bug, that does not affect the execution of the script, and which should be
# solved soon, see: https://github.com/astropy/astropy/issues/9646
# ------------------------------------------------------------------------------------
import os
import sys
import math
import pyvo
from pyvo.dal import tap
import requests
import cgi
import re
import json
import getpass

# Decide what you want to download:
mode_requested = "raw2master"  # other choice: raw2raw

# Set a User Agent (modify as you like, but please let intact the python version used for our usage statistics):
thisscriptname = os.path.basename(__file__)
headers = {}
headers = {'User-Agent': '%s (ESO script drc %s)' % (requests.utils.default_headers()['User-Agent'], thisscriptname)}


# If instead of the requests package, urllib is used:
# headers={'User-Agent': '%s (ESO script drc %s)'%(urllib.request.URLopener.version, thisscriptname)}

def getToken(username, password):
    """Token based authentication to ESO: provide username and password to receive back a JSON Web Token."""
    if username == None or password == None:
        return None
    token_url = "https://www.eso.org/sso/oidc/token"
    token = None
    try:
        response = requests.get(token_url,
                                params={"response_type": "id_token token", "grant_type": "password",
                                        "client_id": "clientid",
                                        "username": username, "password": password})
        token_response = json.loads(response.content)
        token = token_response['id_token'] + '=='
    except NameError as e:
        print(e)
    except:
        print("*** AUTHENTICATION ERROR: Invalid credentials provided for username %s" % (username))

    return token


def download_asset(file_url, filename=None, token=None):
    headers = None
    if token != None:
        headers = {"Authorization": "Bearer " + token}
        response = requests.get(file_url, headers=headers)
    else:
        # Trying to download anonymously
        response = requests.get(file_url, stream=True, headers=headers)

    if filename == None:
        contentdisposition = response.headers.get('Content-Disposition')
        if contentdisposition != None:
            value, params = cgi.parse_header(contentdisposition)
            filename = params["filename"]

        if filename == None:
            # last chance: get anything after the last '/'
            filename = url[url.rindex('/') + 1:]

    if response.status_code == 200:
        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=50000):
                f.write(chunk)

    return (response.status_code, filename)


def calselector_info(description):
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


def print_calselector_info(description, mode_requested):
    """Print the most relevant params contained in the main calselector description."""

    category, complete, certified, mode_executed, messages = calselector_info(description)

    alert = ""
    if complete != "true":
        alert = "ALERT: incomplete calibration cascade"

    mode_warning = ""
    if mode_executed != mode_requested:
        mode_warning = "WARNING: requested mode (%s) could not be executed" % (mode_requested)

    certified_warning = ""
    if certified != "true":
        certified_warning = "WARNING: certified=\"%s\"" % (certified)

    print("    calibration info:")
    print("    ------------------------------------")
    print("    science category=%s" % (category))
    print("    cascade complete=%s" % (complete))
    print("    cascade messages=%s" % (messages))
    print("    cascade certified=%s" % (certified))
    print("    cascade executed mode=%s" % (mode_executed))
    print("    full description: %s" % (this_description))

    return alert, mode_warning, certified_warning


### MAIN ###

# --- Defining the run from which files will be downloaded:
run_id = '099.B-0118(A)'

print()
print("Looking for SCIENCE frames belonging to the %s observing run." % (run_id))
print()

# --- Authenticate with your ESO User Portal credentials: a token is received

print("Authenticating:")
username = input("Type your ESO username: ")
password = getpass.getpass(prompt="%s user's password: " % (username), stream=None)

# With successful authentication you get a valid token,
# which needs to be added to the HTTP header of your GET request (see later),
token = getToken(username, password)
if token == None:
    print("Could not authenticate. Continuing as anonymous")
else:
    print("Authentication successful")
print()

# --- instantiate the ESO TAP service for raw and processed data:

ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
tapobs = tap.TAPService(ESO_TAP_OBS)

# --- Query TAP for science raw
print("Querying the ESO TAP service at %s" % (ESO_TAP_OBS))

query = """SELECT top 4 dp_id
FROM dbo.raw
WHERE prog_id='%s'
AND dp_cat='SCIENCE'
""" % (run_id)

print("")
print(query)
print("")

rawframes = tapobs.search(query=query)
print(rawframes.to_table())
print("")

# --- Downloading those science raw frames and the calibration files necessary to calibrate them

nfiles = 0
calib_association_trees = []

print("Downloading %d files and their calibration cascades (trying with %s mode)." % (
len(rawframes.to_table()), mode_requested))
print("Note: Even if present in the cascade, siblings are not downloaded by this script.")
print("")

# FOR EACH SCIENCE RAW FRAME...
# -----------------------------
for raw in rawframes:
    nfiles += 1
    mode_different_than_requested = ""

    #  ... 1.- RETRIEVE THE SCIENCE RAW FRAME:
    # ---------------------------------------
    sci_url = "https://dataportal.eso.org/dataportal_new/file/%s" % raw["dp_id"].decode()
    status, sci_filename = download_asset(sci_url, token=token)
    print("SCIENCE: %4d/%d dp_id: %s downloaded" % (nfiles, len(rawframes), sci_filename))

    #  ... 2.- GET ITS LIST OF CALIBRATION FILES (without saving it into a file):
    # ---------------------------------------
    calselector_url = "http://archive.eso.org/calselector/v1/associations?dp_id=%s&mode=%s&responseformat=votable" % (
    raw["dp_id"].decode(), mode_requested)
    datalink = pyvo.dal.adhoc.DatalinkResults.from_result_url(calselector_url)

    #  ... 3.- PRINT CASCADE INFORMATION AND MAIN DESCRIPTION
    # ---------------------------------------
    this_description = next(datalink.bysemantics('#this')).description
    alert, mode_warning, certified_warning = print_calselector_info(this_description, mode_requested)

    # create and use a mask to get only the #calibration entries:
    calibrators = datalink['semantics'] == b'#calibration'
    calib_urls = datalink.to_table()[calibrators]['access_url', 'eso_category']

    #  ... 4.- DOWNLOAD EACH #calibration FILE FROM THE LIST
    # ---------------------------------------
    i_calib = 0
    for url, category in calib_urls:
        i_calib += 1
        # Calibration files are accessible anonymously
        status, filename = download_asset(url.decode())
        if status == 200:
            print(
                "    CALIB: %4d/%d dp_id: %s (%s) downloaded" % (i_calib, len(calib_urls), filename, category.decode()))
        else:
            print("    CALIB: %4d/%d dp_id: %s (%s) NOT DOWNLOADED (http status:%d)" % (
            i_calib, len(calib_urls), filename, category.decode(), status))

    #  ... 5.- PRINT ANY ALERT OR WARNING ENCOUNTERED BY THE PROCESS THAT GENERATES THE CALIBRATION CASCADE
    # ---------------------------------------
    if alert != "":
        print("    %s" % (alert))
    if mode_warning != "":
        print("    %s" % (mode_warning))
    if certified_warning != "":
        print("    %s" % (certified_warning))

    print("------------------------------------------------------------------------------------------------")

