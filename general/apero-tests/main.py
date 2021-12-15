"""
Main script to run the apero tests

@author: charles
"""
import os
from datetime import datetime
from pathlib import Path

from apero.core import constants
from jinja2 import Environment, FileSystemLoader, select_autoescape

import apero_tests.utils as ut
from apero_tests.drs_test import OUTDIR, TEMPLATEDIR, DrsTest
from apero_tests.spirou.test_definitions import tests

# TODO: Make independent of instrument
params = constants.load("SPIROU")
setup = os.environ["DRS_UCONFIG"]  # setup
instrument = params["INSTRUMENT"]  # instrument


def summary(test: DrsTest):
    """
    Create summary info for a test.

    :param test: Test to summarize
    :type test: DrsTest
    """

    f = open(os.path.join(OUTDIR, test, test + ".html"), "r")
    html = f.read()

    npassed = html.count("Lime")
    ncond = html.count("Yellow")
    nfailed = html.count("Red")

    ntotal = npassed + ncond + nfailed

    if npassed == ntotal:
        color = "Lime"
    elif ncond >= 1 and nfailed == 0:
        color = "Yellow"
    else:
        color = "Red"

    return ntotal, npassed, ncond, nfailed, color


date_ini = datetime.now()  # initial date and time

# =============================================================================
# Run the tests
# =============================================================================
n = len(tests)  # number of tests
for i, test in enumerate(tests):

    # Create dir only if it don't exist (ignore error auto if it's there)
    p = Path(OUTDIR, test.test_id)
    p.mkdir(parents=True, exist_ok=True)

    print("Test {0}/{1}".format(i + 1, n))
    print("Testing {0}.py".format(test.recipe_name))

    test.run_test()

    print()

print("All tests done")

# =============================================================================
# Write HTML summary
# =============================================================================
# build HTML table element
summary_list = []
for i, test in enumerate(tests):
    test_path = Path(OUTDIR, test.test_id, test.test_id + ".html")
    if test_path.is_file():
        ntotal, npassed, ncond, nfailed, color = summary(test.test_id)
        inspect = os.path.join(*test_path.parts[-2:])  # rel. path for inspect
    else:
        ntotal = npassed = ncond = nfailed = color = inspect = None
    test_dict = {
        "name": test.name,
        "color": color,
        "ntotal": ntotal,
        "passed": npassed,
        "ncond": ncond,
        "nfailed": nfailed,
        "inspect": inspect,
    }
    summary_list.append(test_dict)

date_final = datetime.now()  # final date

delta_date = date_final - date_ini

date_final = date_final.strftime("%Y-%m-%d %H:%M:%S")

# This dictionary is used to fill jinja template with real values
html_dict = {
    "setup": setup,
    "instrument": instrument,
    "date_final": date_final,
    "delta_date": delta_date,
    "summary_list": summary_list,
}

# build main .html summary
env = Environment(
    loader=FileSystemLoader(TEMPLATEDIR),
    autoescape=select_autoescape(["html", "xml"]),
)
env.tests["series"] = ut.is_series

summary_file = "summary.html"

template = env.get_template(summary_file)

html_text = template.render(html_dict)

with open(os.path.join(OUTDIR, summary_file), "w") as f:
    f.write(html_text)
