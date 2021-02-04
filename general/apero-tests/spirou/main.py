"""
Main script to run apero tests

@author: charles
"""
import os
from pathlib import Path
from datetime import datetime

from jinja2 import Environment, FileSystemLoader, select_autoescape
from apero.core import constants
 
from tests import TEMPLATEDIR, OUTDIR
from tests.factory import get_test, TESTDICT
from tests.utils import summary

# =============================================================================
# Define constants
# =============================================================================
params = constants.load('SPIROU')

setup = os.environ['DRS_UCONFIG']  # setup
instrument = params['INSTRUMENT']  # instrument
date_ini = datetime.now()          # initial date

test_list_long = [
        'Preprocessing Recipe Test #1',
        'Dark Master Recipe Test #1',
        'Bad Pixel Correction Recipe Test #1',
        'Localisation Recipe Test #1',
        'Shape Master Recipe Test #1',
        'Shape (per night) Recipe Test #1',
        'Flat/Blaze Correction Recipe Test #1',
        'Thermal Correction Recipe Test #1',
        'Master Leak Correction Recipe Test #1',
        'Master Wavelength Solution Recipe Test #1',
        'Wavelength Solution (per night) Test #1',
        'Leak (per night) Correction Test #1',
        'Extraction Recipe Test #1',
        'Extraction Recipe Test #2',
        'Extraction Recipe Test #3',
        'Make Telluric Recipe Test #1',
        'Fit Telluric Recipe Test #1',
        'Make Template Recipe Test #1',
        'CCF Recipe Test #1'
        ]


test_list_short = list(TESTDICT)

# =============================================================================
# Run the tests
# =============================================================================
n = len(test_list_short)  # number of tests
for i, testid in enumerate(test_list_short):

    # Generate test object
    test = get_test(testid)

    # Create dir only if not exist (ignore error auto if it's there)
    p = Path(OUTDIR, testid)
    p.mkdir(parents=True, exist_ok=True)

    print('test {0}/{1}'.format(i+1, n))
    print('running {0}.py'.format(test_list_short[i]))

    if test is not None:
        test.runtest()
    else:
        print('test {} not implemented'.format(testid))

    print()

print('all tests done')

# =============================================================================
# Write html summary
# =============================================================================
# build table element
summary_list = []
for i in range(n):
    test_path = Path(OUTDIR, test_list_short[i], test_list_short[i]+'.html')
    if test_path.is_file():
        ntotal, npassed, ncond, nfailed, color = summary(test_list_short[i])
        inspect = '/'.join(test_path.parts[-2:])  # rel. path for inspect
    else:
        ntotal = npassed = ncond = nfailed = color = inspect = None
    test_dict = {
            'name': test_list_long[i],
            'color': color,
            'ntotal': ntotal,
            'passed': npassed,
            'ncond': ncond,
            'nfailed': nfailed,
            'inspect': inspect
            }
    summary_list.append(test_dict)

date_final = datetime.now()  # final date

delta_date = date_final - date_ini

date_final = date_final.strftime("%Y-%m-%d %H:%M:%S")

html_dict = {
        'setup': setup,
        'instrument': instrument,
        'date_final': date_final,
        'delta_date': delta_date,
        'summary_list': summary_list,
        }

# build main .html doc
env = Environment(
    loader=FileSystemLoader(TEMPLATEDIR),
    autoescape=select_autoescape(['html', 'xml'])
)

summary_file = 'summary.html'

template = env.get_template(summary_file)

html_text = template.render(html_dict)

with open(os.path.join(OUTDIR, summary_file), 'w') as f:
    f.write(html_text)
