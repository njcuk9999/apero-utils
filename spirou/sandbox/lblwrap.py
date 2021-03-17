"""
Thin wrapper for the lbl codes with convienient options.
"""
import argparse

import numpy as np
from astropy.table import Table

from lbl import lbl
from compilbl import compilbl
from zeropointcorr import zpcorr
from diagnosticplots import pdfplots

parser = argparse.ArgumentParser(
    description="Thin wrapper to analyse list of objects with lbl."
)
parser.add_argument(
    "filename",
    type=str,
    help="File containing csv table with an object and a reference.",
)
parser.add_argument(
    "-l",
    "--skip-lbl",
    dest="do_lbl",
    action="store_false",
    help="Skip the lbl (first step) calculations.",
)
parser.add_argument(
    "-c",
    "--skip-compil",
    dest="do_compil",
    action="store_false",
    help="Skip the compilbl (second step) calculations.",
)
parser.add_argument(
    "-z",
    "--zp-file",
    dest="zp_file",
    default=None,
    help="File with zero-point corrections.",
)
parser.add_argument(
    "-p",
    "--plots",
    dest="do_plots",
    action="store_true",
    help="Produce plots of lbl output (see diagnosticplots.py).",
)
parser.add_argument(
    "-f",
    "--force",
    dest="force",
    action="store_true",
    help="Force recalculations if files exist (applied to each step).",
)
clargs = parser.parse_args()

tbl = Table.read(clargs.listfile, format="csv")

random_order = np.argsort(np.random.random(len(tbl)))
tbl = tbl[random_order]

for i in range(len(tbl)):
    try:
        if clargs.do_lbl:
            lbl(tbl["OBJECT"][i], obj_template=tbl["TEMPLATE"][i], force=clargs.force)
        if clargs.do_compil:
            compilbl(
                tbl["OBJECT"][i], obj_template=tbl["TEMPLATE"][i], force=clargs.force
            )
        if clargs.zp_file is not None:
            # NOTE: this is a (very) slightly modified zpcorr that takes full filename
            #  the file is in the apero-utils repo
            zpcorr(
                clargs.zp_file,
                tbl["OBJECT"][i],
                obj_template=tbl["TEMPLATE"][i],
                force=clargs.force,
            )

        # GL* are named "GJ "* in exoplanet archive
        if tbl["OBJECT"][i].startswith("GL"):
            exoarchive_name = tbl["OBJECT"][i].replace("GL", "GJ ")
        elif tbl["OBJECT"][i].startswith("GJ"):
            exoarchive_name = tbl["OBJECT"][i].replace("GJ", "GJ ")
        else:
            exoarchive_name = tbl["OBJECT"][i]
        if clargs.do_plots:
            pdfplots(
                tbl["OBJECT"][i],
                obj_template=tbl["TEMPLATE"][i],
                exoarchive_name=exoarchive_name,
                doplot_debug=False,
            )
    except Exception as e:
        print(e)
