"""
Init file with definitions and imports in parent package namespace.

@author: vandalt
"""
import os

# Parent directory of the tests module
# NOTE: this should be if the tests are merged with APERO or distributed
# as a package
PARENTDIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
TEMPLATEDIR = os.path.join(PARENTDIR, 'templates')
OUTDIR = os.path.join(PARENTDIR, 'out')
