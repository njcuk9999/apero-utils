#### file:        eso_programmatic.py
#### description: ESO python library defining methods to programmatically access the ESO science archive
#### status:      alpha release
#### version:     0.3
#### date:        2021-07-29
#### type:        not yet a package, just only an evolving library of methods
#### contact:     https://support.eso.org/
#### ---------------------------------------------------------------------------------------------------

__name__ = "eso_programmatic"
__version__ = "0.3"

##   - Definition of ESO services' endpoints

ESO_TAP_OBS              = "http://archive.eso.org/tap_obs"
TOKEN_AUTHENTICATION_URL = "https://www.eso.org/sso/oidc/token"

##   - getToken: Method to authenticate with ESO credentials;
##               Returns: a token (or None)

import requests
import json
def getToken(username, password):
    """Token based authentication to ESO: provide username and password to receive back a JSON Web Token."""
    if username==None or password==None:
        return None

    token = None
    try:
        response = requests.get(TOKEN_AUTHENTICATION_URL,
                            params={"response_type": "id_token token",
                                    "grant_type":    "password",
                                    "client_id":     "clientid",
                                    "username":      username,
                                    "password":      password})
        token_response = json.loads(response.content)
        token = token_response['id_token']+'=='
    except NameError as e:
        print(e)
    except:
        print("*** AUTHENTICATION ERROR: Invalid credentials provided for username %s" %(username))
    
    return token

##   - createSession(token=None):

def createSession(token=None):
    session = requests.Session()
    if token:
        session.headers['Authorization'] = "Bearer " + token
    return session

##   - downloadURL(file_url[, dirname, filename, session]): Method to download a file given its URL,
# either anonymously or with a token.
# Returns: http status, filepath on disk (if successul)
import cgi
import os
import sys
def downloadURL(file_url, dirname='.', filename=None, session=None):
    """Method to download a file, either anonymously (no session or session not "tokenized"), or authenticated (if session with token is provided).
       It returns: http status, and filepath on disk (if successful)"""

    if dirname != None:
        if not os.access(dirname, os.W_OK):
            print("ERROR: Provided directory (%s) is not writable" % (dirname))
            sys.exit(1)
      
    if session!=None:
        response = session.get(file_url, stream=True)
    else:
        # no session -> no authentication
        response = requests.get(file_url, stream=True)

    # If not provided, define the filename from the response header
    if filename == None:
        contentdisposition = response.headers.get('Content-Disposition')
        if contentdisposition != None:
            value, params = cgi.parse_header(contentdisposition)
            filename = params["filename"]

        # if the response header does not provide a name, derive a name from the URL
        if filename == None:
            # last chance: get anything after the last '/'
            filename = file_url[file_url.rindex('/')+1:]

    # define the file path where the file is going to be stored
    if dirname == None:
        filepath = filename
    else:
        filepath = dirname + '/' + filename

    if response.status_code == 200:
        with open(filepath, 'wb') as f:
            for chunk in response.iter_content(chunk_size=50000):
                f.write(chunk)

    return (response.status_code, filepath)


# Let's define some methods to nicely print reference files information from the calselector service:

# calselectorInfo(description): [internal] parsing a calselector description 
import re
def calselectorInfo(description):
    """Parse the main calSelector description, and fetch: category, complete, certified, mode, and messages."""

    category=""
    complete=""
    certified=""
    mode=""
    messages=""

    m = re.search('category="([^"]+)"', description)
    if m:
        category=m.group(1)
    m = re.search('complete="([^"]+)"', description)
    if m:
        complete=m.group(1).lower()
    m = re.search('certified="([^"]+)"', description)
    if m:
        certified=m.group(1).lower()
    m = re.search('mode="([^"]+)"', description)
    if m:
        mode=m.group(1).lower()
    m = re.search('messages="([^"]+)"', description)
    if m:
        messages=m.group(1)

    return category, complete, certified, mode, messages

##   - printCalselectorInfo(description, mode_requested): method that returns
##     possible alerts and warnings on the obtained calibration cascade,
##     while printing most relevant info.

def printCalselectorInfo(description, mode_requested):
    """Print the most relevant params contained in the main calselector description."""

    category, complete, certified, mode_executed, messages = calselectorInfo(description)

    alert=""
    if complete!= "true":
        alert = "ALERT: incomplete calibration cascade"

    mode_warning=""
    if mode_executed != mode_requested:
        mode_warning = "WARNING: requested mode (%s) could not be executed" % (mode_requested)

    certified_warning=""
    if certified != "true":
        certified_warning = "WARNING: certified=\"%s\"" %(certified)

    print("    calibration info:")
    print("    ------------------------------------")
    print("    science category=%s" % (category))
    print("    cascade complete=%s" % (complete))
    print("    cascade messages=%s" % (messages))
    print("    cascade certified=%s" % (certified))
    print("    cascade executed mode=%s" % (mode_executed))
    print("    full description: %s" % (description))

    return alert, mode_warning, certified_warning

##   - printTableTransposedByTheRecord: Utility method to print a table transposed, one record at the time    

def printTableTransposedByTheRecord(table):
    """Utility method to print a table transposed, one record at the time"""
    prompt='    '
    rec_sep='-' * 105
    print('=' * 115)
    for row in table:
        for col in row.columns:
            print("{0}{1: <14} = {2}".format(prompt, col, row[col]) )
        print("{0}{1}".format(prompt,rec_sep))


