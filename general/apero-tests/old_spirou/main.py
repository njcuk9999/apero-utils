"""
Main script to run the apero tests

@author: charles
"""
import os
from pathlib import Path
from datetime import datetime

#from jinja2 import Environment, FileSystemLoader, select_autoescape
#from apero.core import constants
 
from test_definitions import tests
from tests import TEMPLATEDIR, OUTDIR
#from tests.factory import get_test, TESTDICT
#from tests.utils import summary

date_ini = datetime.now() # initial date

# =============================================================================
# Run the tests
# =============================================================================
n = len(tests)  # number of tests
for i, test in enumerate(tests):

    # Create dir only if it don't exist (ignore error auto if it's there)
    p = Path(OUTDIR, test.name)
    p.mkdir(parents=True, exist_ok=True)

    print('Test {0}/{1}'.format(i+1, n))
    print('Running {0}.py'.format(test.test_id))

    if test.recipe is not None:
        test.runtest()
    else:
        print('Test {} not implemented'.format(test.test_id))

    print()

print('All tests done')

# =============================================================================
# Write html summary
# =============================================================================
# build table element
#summary_list = []
#for i in range(n):
#    test_path = Path(OUTDIR, test_list_short[i], test_list_short[i]+'.html')
#    if test_path.is_file():
#        ntotal, npassed, ncond, nfailed, color = summary(test_list_short[i])
#        inspect = '/'.join(test_path.parts[-2:])  # rel. path for inspect
#    else:
#        ntotal = npassed = ncond = nfailed = color = inspect = None
#    test_dict = {
#            'name': test_list_long[i],
#            'color': color,
#            'ntotal': ntotal,
#            'passed': npassed,
#            'ncond': ncond,
#            'nfailed': nfailed,
#            'inspect': inspect
#            }
#    summary_list.append(test_dict)

#date_final = datetime.now()  # final date

#delta_date = date_final - date_ini

#date_final = date_final.strftime("%Y-%m-%d %H:%M:%S")

#html_dict = {
#        'setup': setup,
#        'instrument': instrument,
#        'date_final': date_final,
#        'delta_date': delta_date,
#        'summary_list': summary_list,
#        }

# build main .html doc
#env = Environment(
#    loader=FileSystemLoader(TEMPLATEDIR),
#    autoescape=select_autoescape(['html', 'xml'])
#)

#summary_file = 'summary.html'

#template = env.get_template(summary_file)

#html_text = template.render(html_dict)

#with open(os.path.join(OUTDIR, summary_file), 'w') as f:
#    f.write(html_text)