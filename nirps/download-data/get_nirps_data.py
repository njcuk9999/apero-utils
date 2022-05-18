import argparse
import os
import re
from typing import Optional

import eso_utils as eut
import yaml
from astropy.table import unique, vstack

# Constants
# RAMPS_DIR = "raw-ramps"
# TMP_DIR = "raw-dump"
# READS_DIR = "test-reads"
# CALIB_TYPE = "raw2raw"
# INSTRUMENT = "NIRPS"  # Instrument to get data from
# DATESTR = "2022"  # Get data after that date (inclusively)
# less HARPS data and more varied files. Use for testing
# INSTRUMENT = "HARPS"  # Instrument to get data from
# DATESTR = "2021"  # Get data after that date (inclusively)


def get_nirps_data(
    program_id,
    instrument: str = "NIRPS",
    data_category: str = "ramps",
    calib_mode: Optional[str] = None,
    tmp_dir: str = "tmp",
    auth_file: str = "auth.yaml",
    destination: Optional[str] = None,
    unzip: bool = True,
    start: Optional[str] = None,
    test_run: bool = False,
    per_night: bool = True,
    overwrite: bool = False,
    continue_tmp: bool = False,
    cache_tmp: bool = False,
):

    session = eut.EsoSession(auth_file)

    if data_category == "ramps":
        dp_cat = "~TECHNICAL"
    elif data_category == "reads":
        dp_cat = "TECHNICAL"
    else:
        raise ValueError("Unknown data category.")

    files_tbl = eut.get_eso_tbl(
        program_id=program_id,
        instrument=instrument,
        start=start,
        dp_cat=dp_cat,
        session=session,
        sun=False,
    )

    if calib_mode is not None:

        if calib_mode not in ["raw2raw", "raw2master"]:
            raise ValueError("calib_mode must be raw2raw or raw2master")

        # NOTE: This often crashes, probably because many IDs in one query.
        # TODO: Would be good to have an alternative that is not iteration w/ 2 requests
        try:
            calibs_tbl = eut.get_calib_table(files_tbl, calib_mode, session=session)
        except ValueError:
            print("Calib table query failed... Trying a second time")

            try:
                calibs_tbl = eut.get_calib_table(files_tbl, calib_mode, session=session)
            except ValueError:
                print("Calib table query failed twice. Trying the iterative way.")

                # This fallback function gets calibs one-by-one.
                # It is slow because there are 2 HTTP requests per file
                calibs_tbl = eut.get_calibs_iter(files_tbl, calib_mode, session=session)

        if len(calibs_tbl) > 0:
            files_tbl = unique(
                vstack([files_tbl, calibs_tbl]), keys=["access_url", "dp_id", "dp_cat"]
            )

    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)

    files_tbl["local_night"] = eut.get_night_from_date(files_tbl["date_obs"])

    downloaded_files = dict()
    for row in files_tbl:
        access_url = row["access_url"]
        id_for_path = re.sub("\\.|\\:", "_", row["dp_id"])
        local_path = os.path.join(tmp_dir, id_for_path + ".fits.Z")
        if destination is not None:
            if per_night:
                dest_dir = os.path.join(destination, row["local_night"])
            else:
                dest_dir = destination
            final_destination = os.path.join(dest_dir, id_for_path + ".fits")
            if not unzip:
                final_destination = final_destination + ".Z"

        if (
            destination is None
            or (destination is not None and not os.path.isfile(final_destination))
            or overwrite
        ):
            if destination is not None:
                downloaded_files[local_path] = final_destination
            if not test_run:
                eut.download_file(
                    access_url,
                    local_path,
                    session=session,
                    continuation=continue_tmp,
                    cache=cache_tmp,
                )
            else:
                if destination is not None:
                    print(local_path, "->", final_destination)
                else:
                    print(local_path)
        elif destination is not None:
            print(
                f"File {final_destination} exists and overwrite=False. Skipping download."
            )

    if not test_run and destination is not None:
        eut.dispatch_files(downloaded_files, unzip=unzip)

    return files_tbl


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get NIRPS data from ESO Archive.")
    parser.add_argument(
        "category",
        type=str,
        help="Type of data to download",
        choices=["ramps", "reads"],
    )
    parser.add_argument(
        "-m",
        "--calib-mode",
        dest="mode",
        type=str,
        default=None,
        help="Calibration mode from ESO",
        choices=["raw2raw", "raw2master", "None"]
    )
    parser.add_argument(
        "-s",
        "--start",
        type=str,
        help="Start date",
        default="2022-05-05",
    )
    parser.add_argument(
        "--test",
        dest="test_run",
        action="store_true",
        help="Do not download data on disk, just get file info.",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        default="config.yaml",
        help="Config file with info about dir structure and obs program.",
    )
    args = parser.parse_args()

    config_file = args.config
    with open(config_file, "r") as authfile:
        config_info = yaml.safe_load(authfile)

    if args.category == "ramps":
        data_category = "ramps"
        tmp_dir = config_info["raw_tmp"]
        destination = config_info["raw_dir"]
        if args.mode is None:
            calib_mode = "raw2raw"
        elif args.mode == "None":
            calib_mode = None
        else:
            calib_mode = args.mode
        per_night = True
        unzip = True

    if args.category == "reads":
        data_category = "reads"
        tmp_dir = config_info["reads_tmp"]
        destination = config_info["reads_dir"]
        if args.mode is None:
            calib_mode = None
        elif args.mode == "None":
            calib_mode = None
        else:
            calib_mode = args.mode
        per_night = False
        unzip = False

    nirps_files = get_nirps_data(
        config_info["program_id"],
        instrument="NIRPS",
        data_category=data_category,
        calib_mode=calib_mode,
        tmp_dir=tmp_dir,
        destination=destination,
        per_night=per_night,
        start=args.start,
        test_run=args.test_run,
        cache_tmp=True,
        continue_tmp=True,
        overwrite=False,
        unzip=unzip,
    )
