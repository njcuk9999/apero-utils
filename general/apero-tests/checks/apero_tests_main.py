"""
Main script to run apero tests

@author: charles
"""
import os
from pathlib import Path
from datetime import datetime

from jinja2 import Environment, FileSystemLoader, select_autoescape

from apero_test_factory import TESTDICT, get_test
from apero_tests_func import summary
from apero.core import constants

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
        'Thermal Correction Test #1',
        'Master Leak Correction Recipe Test #1',
        'Leak (per night) Correction Test #1',
        'Master Wavelength Solution Recipe Test #1',
        'Wavelength Solution (per night) Test #1',
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
    p = Path(os.path.join('..', 'out', testid))
    p.mkdir(exist_ok=True)

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
# TODO: Move this to jinja2 template if not too complicated
html_table = []
for i in range(n):
    if os.path.isfile("{0}/{0}.html".format(test_list_short[i])):
        html_str, color = summary(test_list_short[i])
        html_table.append("""
        <tr>
          <td>{1}</td>
          <td>{2}</td>
          <td bgcolor={3}><a href='{0}/{0}.html'>Inspect</a></td>
        </tr>
        """.format(test_list_short[i], test_list_long[i], html_str, color))
    else:
        html_table.append("""
        <tr>
          <td>{0}</td>
          <td></td>
          <td></td>
        </tr>
        """.format(test_list_long[i]))
html_table = "".join(html_table)


date_final = datetime.now()  # final date

delta_date = date_final - date_ini

date_final = date_final.strftime("%Y-%m-%d %H:%M:%S")

html_dict = {
        'setup': setup,
        'instrument': instrument,
        'date_final': date_final,
        'delta_date': delta_date,
        'html_table': html_table,
        }

# build main .html doc
# TODO: Use package or define path in a main/init file?
env = Environment(
    loader=FileSystemLoader('../templates'),
    autoescape=select_autoescape(['html', 'xml'])
)

template = env.get_template('summary.html')

html_text = template.render(html_dict)

with open(os.path.join('../out', 'summary.html'), 'w') as f:
    f.write(html_text)
