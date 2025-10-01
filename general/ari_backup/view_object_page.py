#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-09-23 at 12:07

@author: cook
"""
import os
import shutil
import tempfile
from sphinx.application import Sphinx


# =============================================================================
# Define variables
# =============================================================================
# define the paths for each apero profile
APERO_PROFILES = dict()
# APERO_PROFILES['nirps_he_online'] = '/cosmos99/nirps/apero-data/nirps_he_online/other/ari/nirps_he_online_udem/object_pages/'
# APERO_PROFILES['nirps_ha_online'] = '/cosmos99/nirps/apero-data/nirps_ha_online/other/ari/nirps_ha_online_udem/object_pages/'
APERO_PROFILES['spirou_offline'] = '/cosmos99/spirou/apero-data/spirou_offline/other/ari/spirou_offline_udem/object_pages/'
# -----------------------------------------------------------------------------
MINI_CONF = """project = 'TempDocs'
extensions = []
master_doc = 'index'
"""

OUTDIR = '/cosmos99/nirps/home/ARI-MANUAL'

# =============================================================================
# Start of functions
# =============================================================================
def build_with_sphinx(rst_dir, out_dir):
    # make a temporary project directory
    tmpdir = tempfile.mkdtemp()
    srcdir = os.path.join(tmpdir, "source")
    os.makedirs(srcdir)

    # copy rst files into source/
    shutil.copytree(rst_dir, srcdir, dirs_exist_ok=True)

    # minimal conf.py
    with open(os.path.join(srcdir, "conf.py"), "w") as f:
        f.write(MINI_CONF)

    # create index.rst if missing
    index = os.path.join(srcdir, "index.rst")
    if not os.path.exists(index):
        # include all rst files (except index itself)
        files = [f[:-4] for f in os.listdir(srcdir)
                 if f.endswith(".rst") and f != "index.rst"]
        toc = "\n   ".join(files)
        with open(index, "w") as f:
            f.write("TempDocs\n=========\n\n.. toctree::\n   :maxdepth: 2\n\n   " + toc)

    # prepare output folders
    doctreedir = os.path.join(tmpdir, "doctrees")
    os.makedirs(doctreedir)

    # run sphinx
    app = Sphinx(
        srcdir=srcdir,
        confdir=srcdir,
        outdir=out_dir,
        doctreedir=doctreedir,
        buildername="html",
    )
    app.build(force_all=True)

    # clean up temporary directory
    shutil.rmtree(tmpdir)

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":

    # ----------------------------------------------------------------------
    for apero_profile in APERO_PROFILES:
        print(f"Apero profile: {apero_profile}")

        # get list of object directories
        object_names = [d for d in os.listdir(APERO_PROFILES[apero_profile])
                        if os.path.isdir(os.path.join(APERO_PROFILES[apero_profile], d))]

        for object_name in object_names:

            # -------------------------------------------------------------------------
            # construct the file path
            in_path = os.path.join(APERO_PROFILES[apero_profile], object_name)
            out_path = os.path.join(OUTDIR, apero_profile, object_name)
            os.makedirs(out_path, exist_ok=True)

            # check if the file exists
            if not os.path.exists(in_path):
                print(f"Directory {in_path} does not exist.")
                continue

            # Build the page using sphinx
            build_with_sphinx(in_path, out_path)

# =============================================================================
# End of code
# =============================================================================
