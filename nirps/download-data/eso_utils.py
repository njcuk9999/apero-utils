"""
This is a modified version of the eso_programmatic.py module. Some edits are
just formatting. Others aim to make the tool a bit more flexible, for example
regarding file path formats. Also adding functionality that overlaps with
astroquery

Original file by @almicol/ESO
Modified by @vandalt (Thomas Vandal)
"""
import cgi
import datetime
import io
import json
import os
import re
import shutil
import subprocess
import warnings
from typing import Dict, List, Optional, Union

import astropy.utils.data
import pytz
import pyvo
import requests
import tqdm
import yaml
from astropy.table import Column, Table, unique, vstack
from astropy.time.core import Time
from requests import Session

TOKEN_AUTHENTICATION_URL = "https://www.eso.org/sso/oidc/token"
ESO_TAP_OBS = "http://archive.eso.org/tap_obs"


class EsoSession(Session):
    def __init__(self, auth_file):

        super(EsoSession, self).__init__()

        self.auth_file = auth_file

        self.set_token()

    def set_token(self):

        with open(self.auth_file, "r") as authfile:
            auth_info = yaml.safe_load(authfile)

        eso_token_lifetime = 8
        self.expiration_time = datetime.datetime.now() + datetime.timedelta(
            hours=eso_token_lifetime
        )
        self.token = get_token(auth_info["user"], auth_info["password"])
        self.headers["Authorization"] = "Bearer " + self.token

    def pre_request_check(self):
        # 15 minutes buffer for token validity
        if (
            self.token is None
            or datetime.datetime.now() + datetime.timedelta(minutes=15)
            >= self.expiration_time
        ):
            self.set_token()

    def send(self, request, **kwargs):
        self.pre_request_check()
        return super(EsoSession, self).send(request, **kwargs)


def get_night_from_date(date: Union[Time, Column], fmt="isot"):
    if hasattr(date, "dtype"):
        if date.dtype == "object":
            date = date.astype(str)

    # date_obs is in UTC time
    time_utc = Time(date, format=fmt, scale="utc")
    # Convert to chile time. Pytz is aware of date at which UTC offset changes
    datetime_cl = time_utc.to_datetime(pytz.timezone("America/Santiago"))
    night_list = []
    start_night = datetime.time(15, 00, 00)
    for dtcl in datetime_cl:
        # Drop timezone-aware info (value is still chile time)
        dtcl = dtcl.replace(tzinfo=None)

        if dtcl.time() >= start_night:
            night = dtcl.date()
        else:
            night = dtcl.date() - datetime.timedelta(days=1)

        night_list.append(night.isoformat())

    return night_list


def get_eso_tbl(
    data_type: str = "raw",
    program_id: Optional[str] = None,
    instrument: Optional[str] = None,
    start: Optional[str] = None,
    end: Optional[str] = None,
    cols: Optional[Union[List[str], str]] = None,
    dp_cat: Optional[Union[str, List[str]]] = None,
    session: Optional[Session] = None,
    query_async: bool = False,
    to_table: bool = True,
    show_query: bool = False,
    sun: bool = False,
):

    tap = pyvo.dal.TAPService(ESO_TAP_OBS, session=session)

    if cols is not None:
        if isinstance(cols, str):
            colstr = cols
        if isinstance(cols, list):
            colstr = ", ".join(cols)
    else:
        colstr = "*"

    if data_type == "raw":
        tbl_str = "dbo.raw"
        prog_col = "prog_id"
        date_col = "date_obs"
        inst_col = "instrument"
    elif data_type == "red":
        tbl_str = "ivoa.ObsCore"
        prog_col = "proposal_id"
        date_col = "obs_release_date"
        inst_col = "instrument_name"

    cond_list = []
    if program_id is not None:
        cond_list.append(f"{prog_col}='{program_id}'")
    if instrument is not None:
        cond_list.append(f"{inst_col}='{instrument}'")
    if start is not None:
        cond_list.append(f"{date_col}>='{start}'")
    if end is not None:
        cond_list.append(f"{date_col}>='{end}'")
    if dp_cat is not None:
        if data_type == "red":
            warnings.warn("dp_cat is not a column for reduced data. Ignoring it.")
        else:
            if isinstance(dp_cat, (tuple, list)):
                if dp_cat[0].startswith("~"):
                    op = " NOT IN "
                    if any([not d.startswith("~") for d in dp_cat]):
                        raise ValueError("Cannot mix include and exclude int dp_cat")
                    dp_cat = [f"'{d[1:]}'" for d in dp_cat]
                else:
                    op = " IN "
                    if any([d.startswith("~") for d in dp_cat]):
                        raise ValueError("Cannot mix include and exclude int dp_cat")
                    dp_cat = [f"'{d}'" for d in dp_cat]
                dp_cat = "(" + ", ".join(dp_cat) + ")"
            else:
                if dp_cat.startswith("~"):
                    op = "<>"
                    dp_cat = f"'{dp_cat[1:]}'"
                else:
                    dp_cat = f"'{dp_cat}'"
                    op = "="
            cond_list.append(f"dp_cat{op}{dp_cat}")
    # TODO: Unsafe: could remove words with sun in it. Need to handle these cases:
    # - Sun, SUN or sun alone
    # - All the above with , before and/or after
    if not sun:
        cond_list.append("object NOT LIKE '%sun%'")

    if len(cond_list) > 0:
        cond_str = " WHERE " + " AND ".join(cond_list)
    else:
        cond_str = ""

    query_str = f"SELECT {colstr} FROM {tbl_str}{cond_str}"

    if show_query:
        print(f"'{query_str}'")

    # TODO: Add async query option
    if not query_async:
        result = tap.search(query=query_str)
    else:
        raise NotImplementedError("Async queries not implemented yet for ESO table")

    return result.to_table() if to_table else result


def get_calibs_iter(
    data_table, calib_mode, session=None, query_async=False, to_table: bool = True
):

    tap = pyvo.dal.TAPService(ESO_TAP_OBS, session=session)
    cal_semantics = f"http://archive.eso.org/rdf/datalink/eso#calSelector_{calib_mode}"

    # NOTE: This is very slow and takes > 1 min for the 28 files I tested with...
    # NOTE: And if we want to make sure all calibs are there can't really skip already-downloaded science files
    all_associated_calibs = []
    # TODO: Make loop faster by generating URL by hand instead of 2nd from_result
    for row in tqdm.tqdm(data_table):

        datalink_url = row["datalink_url"]
        datalink = pyvo.dal.adhoc.DatalinkResults.from_result_url(
            datalink_url, session=session
        ).to_table()
        datalink = datalink[datalink["semantics"] == cal_semantics]

        for dl_row in datalink:
            calib_list_url = dl_row["access_url"]

            # The calselector URL points to a list of files. Get that list of files
            # Ref: http://archive.eso.org/programmatic/rdf/datalink/eso/
            associated_calibs = pyvo.dal.adhoc.DatalinkResults.from_result_url(
                calib_list_url, session=session
            ).to_table()

            sibling_mask = associated_calibs["semantics"] == "#sibling_raw"
            calib_mask = associated_calibs["semantics"] == "#calibration"
            mask = sibling_mask | calib_mask
            associated_calibs = associated_calibs[mask]

            if len(associated_calibs) > 0:
                all_associated_calibs.append(associated_calibs)

    if len(all_associated_calibs) == 0:
        return Table([])

    all_associated_calibs = unique(vstack(all_associated_calibs))

    # Get calib files in same TAP query format than regular file
    calib_ids = []
    for cf in all_associated_calibs:
        calib_ids.append(cf["access_url"].split("ivo://eso.org/ID?")[1].split("&")[0])
    calib_ids = [f"'{cid}'" for cid in calib_ids]
    cal_dp_ids = "(" + ", ".join(calib_ids) + ")"
    calib_query = f"SELECT * FROM dbo.raw WHERE dp_id IN {cal_dp_ids}"
    if not query_async:
        result = tap.search(query=calib_query)
    else:
        raise NotImplementedError("Async queries not implemented yet for ESO table")

    return result.to_table() if to_table else result


def get_calib_table(
    data_table: Table,
    calib_mode: str,
    id_col: str = "dp_id",
    session: Optional[Session] = None,
    query_async: bool = False,
    to_table: bool = True,
):
    if len(data_table) == 0:
        return Table([])

    tap = pyvo.dal.TAPService(ESO_TAP_OBS, session=session)

    ids_str = ",".join(data_table[id_col])

    cal_semantics = f"http://archive.eso.org/rdf/datalink/eso#calSelector_{calib_mode}"

    # Use format of first row, but replace ID with all IDs in table
    template_datalink_url, template_id = data_table[0]["datalink_url", id_col]
    datalink_url = template_datalink_url.replace(template_id, ids_str)

    # This returns a table with datalink info for each table. Keep calibs only
    datalink = pyvo.dal.adhoc.DatalinkResults.from_result_url(
        datalink_url, session=session
    ).to_table()
    datalink = datalink[datalink["semantics"] == cal_semantics]

    if len(datalink) == 0:
        return Table([])

    # Each row in datalink is calselector URL with link that uses dp id
    # Replace dp id with lis tof IDs to query all calib files at once
    url_all_calibs = re.sub(
        f"{id_col}=.*&", f"{id_col}={ids_str}&", datalink[0]["access_url"]
    )

    # NOTE: Do we want only #calibration or also #sibling_raw?
    # WARN: Maybe because big query, but sometimes get error. Re-running a 2nd time seems to work
    all_associated_calibs = pyvo.dal.adhoc.DatalinkResults.from_result_url(
        url_all_calibs, session=session
    ).to_table()
    sibling_mask = all_associated_calibs["semantics"] == "#sibling_raw"
    calib_mask = all_associated_calibs["semantics"] == "#calibration"
    mask = sibling_mask | calib_mask
    all_associated_calibs = all_associated_calibs[mask]

    # Get calib files in same TAP query format than regular file
    calib_ids = []
    for cf in all_associated_calibs:
        calib_ids.append(cf["access_url"].split("ivo://eso.org/ID?")[1].split("&")[0])
    calib_ids = [f"'{cid}'" for cid in calib_ids]
    cal_dp_ids = "(" + ", ".join(calib_ids) + ")"
    calib_query = f"SELECT * FROM dbo.raw WHERE dp_id IN {cal_dp_ids}"
    if not query_async:
        result = tap.search(query=calib_query)
    else:
        raise NotImplementedError("Async queries not implemented yet for ESO table")

    return result.to_table() if to_table else result


def get_token(username, password):
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


def create_eso_session(token=None):
    """
    Create requests session, with Auth token if provided

    :param token: ESO token (lasts 8 hrs), defaults to None
    :type token: str, optional
    """
    session = requests.Session()
    if token:
        session.headers["Authorization"] = "Bearer " + token
    return session


def download_file(
    file_url: str,
    local_filepath: str,
    verbose: bool = False,
    show_progress: bool = True,
    session: Optional[Session] = None,
    timeout: Optional[int] = None,
    continuation: bool = True,
    cache: bool = False,
    method: str = "GET",
    **kwargs,
):
    """
    Simplified version of astroquery.BaseQuery._download_file()

    All extra kwargs are passed to the http requests

    From astroquery: Licensed under a 3-clause BSD style license - see astroquery/LICENSE.rst

    Parameters
    ----------
    url : string
    local_filepath : string
    timeout : int
    continuation : bool
        If the file has already been partially downloaded *and* the server
        supports HTTP "range" requests, the download will be continued
        where it left off. Skipped if file already fully downloaded
    cache : bool
        Provided that:
            - continuation is False OR
            - HTTP "ranges" are not supported OR
        If the file on is found on disk, the download with be skipped if:
            - File sizes match when content length is available from HTTP request
            - If no file length is avaible from the server
    method : "GET" or "POST"
    """

    if session is not None:
        response = session.request(
            method, file_url, timeout=timeout, stream=True, **kwargs
        )
    else:
        # no session -> no authentication
        response = requests.request(
            method, file_url, timeout=timeout, stream=True, **kwargs
        )

    # Get length of content to be downloaded
    response.raise_for_status()
    if "content-length" in response.headers:
        length = int(response.headers["content-length"])
        if length == 0:
            warnings.warn("URL {0} has length=0".format(file_url), RuntimeWarning)
    else:
        length = None

    if (
        os.path.exists(local_filepath)
        and ("Accept-Ranges" in response.headers)
        and continuation
    ):
        open_mode = "ab"

        existing_file_length = os.stat(local_filepath).st_size
        if length is not None and existing_file_length >= length:
            # all done!
            print(
                "Found cached file {0} with expected size {1}.".format(
                    local_filepath, existing_file_length
                )
            )
            response.close()
            return
        elif existing_file_length == 0:
            open_mode = "wb"
        else:
            print(
                "Continuing download of file {0}, with {1} bytes to "
                "go ({2}%)".format(
                    local_filepath,
                    length - existing_file_length,
                    (length - existing_file_length) / length * 100,
                )
            )

            # bytes are indexed from 0:
            # https://en.wikipedia.org/wiki/List_of_HTTP_header_fields#range-request-header
            end = "{0}".format(length - 1) if length is not None else ""
            session.headers["Range"] = "bytes={0}-{1}".format(existing_file_length, end)

            response = session.request(
                method, file_url, timeout=timeout, stream=True, **kwargs
            )
            response.raise_for_status()
            del session.headers["Range"]

    elif cache and os.path.exists(local_filepath):
        if length is not None:
            statinfo = os.stat(local_filepath)
            if statinfo.st_size != length:
                warnings.warn(
                    f"Found cached file {local_filepath} with size {statinfo.st_size} "
                    f"that is different from expected size {length}",
                    RuntimeWarning,
                )
                open_mode = "wb"
            else:
                print(
                    "Found cached file {0} with expected size {1}.".format(
                        local_filepath, statinfo.st_size
                    )
                )
                response.close()
                return
        else:
            print("Found cached file {0}.".format(local_filepath))
            response.close()
            return
    else:
        open_mode = "wb"

    blocksize = astropy.utils.data.conf.download_block_size

    if verbose:
        print(
            f"Downloading URL {file_url} to {local_filepath} with size {length} "
            f"by blocks of {blocksize}"
        )

    bytes_read = 0

    if show_progress:
        progress_stream = None  # Astropy default
    else:
        progress_stream = io.StringIO()

    print(f"Downloading URL {file_url} to {local_filepath} ...")
    with tqdm.tqdm(total=length) as pb:
        with open(local_filepath, open_mode) as f:
            for block in response.iter_content(blocksize):
                f.write(block)
                bytes_read += len(block)
                if length is not None:
                    pb.update(bytes_read if bytes_read <= length else length)
                else:
                    pb.update(bytes_read)
    response.close()
    return response


def dispatch_files(path_mappings: Optional[Dict[str, str]] = None, unzip: bool = True):

    if unzip:
        unzip_cmd = "uncompress"
        if not shutil.which(unzip_cmd):
            raise OSError(f"Command {unzip_cmd} not found... Cannot extract file.")

    for tmp_path, final_path in path_mappings.items():

        dest_dir = os.path.dirname(final_path)
        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)

        if unzip:
            subprocess.run([unzip_cmd, tmp_path])
            tmp_path = tmp_path.rstrip(".Z")
            final_path = final_path.rstrip(".Z")

        shutil.move(tmp_path, final_path)


def get_default_filename(file_url, response=None):
    # If not provided, define the filename from the response header
    if response is not None:
        # Content disposstion gives attachment; filename={filename} so can extract
        contentdisposition = response.headers.get("Content-Disposition")
        if contentdisposition is not None:
            try:
                _, params = cgi.parse_header(contentdisposition)
                filename = params["filename"]
            except KeyError:
                pass

    # if the response header does not provide a name, derive a name from the URL
    if filename is None:
        # last chance: get anything after the last '/'
        filename = file_url.rsplit("/", maxsplit=1)[-1]

    return filename


def get_calselector_info(description):
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

    category, complete, certified, mode_executed, messages = get_calselector_info(
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


def print_table_transposed_by_record(table):
    """Utility method to print a table transposed, one record at the time"""
    prompt = "    "
    rec_sep = "-" * 105
    print("=" * 115)
    for row in table:
        for col in row.columns:
            print("{0}{1: <14} = {2}".format(prompt, col, row[col]))
        print("{0}{1}".format(prompt, rec_sep))
