"""
Utility functions to manage object ID sheet.

@author: vandalt
"""
import numpy as np


def hdr_to_dict(adict, kwd, ahdr):
    """
    Add value from fits header to list in dict with same keyword.

    Args:
        adict   (dict): dictionnary of list where values are added
        kwd      (str): key corresponding to value
        ahdr  (Header): fits header to fetch value from
    Returns:
        Nothing returned, adict is edited directly

    Notes:
        - If key is not in header, nan is added.
        - ahdr could in fact be any dict-like object.
    """
    try:
        adict[kwd].append(ahdr[kwd])
    except KeyError:
        adict[kwd].append(np.nan)
