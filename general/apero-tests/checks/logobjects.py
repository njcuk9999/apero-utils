"""
Log and index objects from the DRS
"""
from astropy.io import fits

class Index:
    """
    DRS index output
    """

    def __init__(self, path: str):
        """__init__.

        :param path: path to index.fits file
        :type path: str
        """

        tbl = fits.getdata(path)
        self.tbl = tbl
        self.len = len(tbl)
        self.filename = tbl['FILENAME']
        self.nights = tbl['NIGHTNAME']
        self.object = tbl['KW_OBJNAME']


class Log:

    def __init__(self, path: str):
        """__init__.

        :param path: Path to log.fits file
        :type path: str
        """

        tbl = fits.getdata(path)
        self.tbl = tbl
        self.len = len(tbl)
        self.recipe = tbl['RECIPE']
        self.QC = tbl['PASSED_ALL_QC']
        self.ENDED = tbl['ENDED']

        self.indexQCtrue = tbl['PASSED_ALL_QC']
        self.indexQCfalse = ~tbl['PASSED_ALL_QC']
        self.indexENDEDtrue = tbl['ENDED']
        self.indexENDEDfalse = ~tbl['ENDED']

        self.nights = tbl['DIRECTORY']
        self.runstr = tbl['RUNSTRING']
        self.args = tbl['ARGS']

        self.QCstr = tbl['QC_STRING']
        self.QCvalue = tbl['QC_VALUES']
        self.ERRORS = tbl['ERRORS']
        self.LOGFILE = tbl['LOGFILE']
