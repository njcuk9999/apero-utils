#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-10-07

@author: cook
"""
import numpy as np


# =============================================================================
# Define variables
# =============================================================================


# =============================================================================
# Define functions
# =============================================================================
def get_files(times, N):
    # get a time-sorted vector. It will be trimmed
    # within the while loop to reject points that
    # are closest in times. The resulting 'times2'
    # is very uniformly distributed with time.
    times2 = times[np.argsort(times)]

    while len(times2) > N:
        dt1 = np.abs(times2 - np.roll(times2, 1))
        dt2 = np.abs(times2 - np.roll(times2, -1))

        dt1[dt2 < dt1] = dt2[dt2 < dt1]

        times2 = np.delete(times2, np.argmin(dt1))

    # find the values that have been kept and construct
    # a mask
    mask = np.array([times[i] in times2 for i in range(len(times))])

    return mask


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # print hello world
    print('Hello World')

# =============================================================================
# End of code
# =============================================================================
