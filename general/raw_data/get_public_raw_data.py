#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-11-10 at 13:55

@author: cook
"""
import glob
import os
import warnings
from typing import List, Optional, Tuple

import numpy as np
from astropy.io import fits
from astroquery.cadc import Cadc
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
DOWNLOAD = False


# =============================================================================
# Define functions
# =============================================================================
def get_dates(case: Optional[str] = None,
              path: Optional[str] = None) -> List[Tuple[float, float]]:
    """
    Get the dates to push into get_data()
    Only possible if data is on disk


    :param case: optional str, if set use predefined paths
    :param path: optional str, if set use this path to search for files
                 must have directory structure

                 path/night1/*.fits
                 path/night2/*.fits
                 path/night3/*.fits

    :return: list of tuples, start and end time (MJD) for each night in path
             directory list
    """
    if path is not None:
        pass
    elif case == 'minidata1':
        path = '/scratch2/spirou/drs-data/common/minidata1'
    else:
        path = '/scratch2/spirou/drs-data/common/minidata2'
    # get a list of directories
    directories = os.listdir(path)
    # storage for dates
    dates = []
    # loop around directories
    for directory in directories:
        # if we are not dealing with a directory continue
        if not os.path.isdir(directory):
            continue
        # print progress
        print(f'Processing {directory}')
        # get file list in sub-directory
        files = glob.glob(os.path.join(path, directory, '*.fits'))
        # start the min and max times as extremely high/low respectively
        mindate = np.inf
        maxdate = -np.inf
        # loop around files
        for filename in tqdm(files):
            # get the header of this file
            hdr = fits.getheader(filename)
            # if we have a minimum date change mindate
            if hdr['MJDATE'] < mindate:
                mindate = float(hdr['MJDATE']) - 1 / 24.0
            # if we have a maximum date change maxdate
            if hdr['MJDATE'] > maxdate:
                maxdate = float(hdr['MJDATE']) + 1 / 24.0
        # append the best found min and max time to dates
        dates.append((mindate, maxdate))
    # return dates for use in get_data()
    return dates


def get_data(case: Optional[str] = None,
             dates: Optional[List[Tuple[float, float]]] = None
             ) -> List[str]:
    """
    Get public CFHT data for specific dates

    :param case: optional str, if 'minidata1' or 'minidata2' we have the dates
                 already
    :param dates: optional list, can use the MJDATE to get dates, dates is in
                  form:

                  [(start0, end0), (start1, end1), ... (start_n, end_n)]

                  where tuples are individual nights of data
                  start0 is an MJD time before first observation
                  end0 is an MJD time after last observation

    :return: list of strings, the download urls for all available files
    """
    # this is the query to get the table
    # noinspection SqlDialectInspection
    query = """
    SELECT 
    	 Observation.sequenceNumber AS odometer,

    	 Plane.productID AS product_id,
    	 COORD1(CENTROID(Plane.position_bounds)) AS RA,
    	 COORD2(CENTROID(Plane.position_bounds)) AS Dec,
    	 Observation.target_name AS Target_Name,
    	 Plane.time_bounds_lower AS Start_Date,
    	 isDownloadable(Plane.publisherID) AS DOWNLOADABLE,
    	 Plane.publisherID AS publisherID
    FROM caom2.Plane AS Plane 
    	JOIN caom2.Observation AS Observation ON Plane.obsID = Observation.obsID 
    WHERE  ( Observation.instrument_name = 'SPIRou' 
    AND Observation.collection = 'CFHT' 
    AND  ( Plane.quality_flag IS NULL OR Plane.quality_flag != 'junk' )

    AND (Plane.productID LIKE '%a' OR Plane.productID LIKE '%o' 
    	 OR Plane.productID LIKE '%f' OR Plane.productID LIKE '%d'
    	 OR Plane.productID LIKE '%c')
    AND (Plane.time_bounds_lower BETWEEN {MINDATE} and {MAXDATE})
    )
    """

    # these are the MJDATE for minidata one and minidata2
    #  (start, end) -/+ an hour
    # these were made using the get_dates() function
    if dates is not None and case is None:
        pass
    elif case == 'minidata1':
        dates = [(58983.000760133335, 58983.780495566665),
                 (59030.00074343334, 59030.76773646667),
                 (59059.00074323334, 59059.75131846667),
                 (59060.000742133336, 59060.768704066664),
                 (59062.00076043334, 59062.77709166666),
                 (59063.00074183333, 59063.734092666666),
                 (59064.00075893333, 59064.77337426667),
                 (59065.00075683334, 59065.77617046666),
                 (59071.00074233334, 59071.762133966666),
                 (59092.00074163334, 59092.77317936666),
                 (59115.00075953334, 59115.77995036666),
                 (59129.00074253334, 59129.78254386666),
                 (59154.00074043334, 59154.78506546666)]
    else:
        dates = [(58983.000760133335, 58983.780495566665),
                 (59030.00074343334, 59030.76773646667),
                 (59059.00074323334, 59059.75131846667),
                 (59060.000742133336, 59060.768704066664),
                 (59062.00076043334, 59062.77709166666),
                 (59063.00074183333, 59063.734092666666),
                 (59064.00075893333, 59064.77337426667),
                 (59065.00075683334, 59065.77617046666),
                 (59071.00074233334, 59071.762133966666),
                 (59092.00074163334, 59092.77317936666),
                 (59115.00075953334, 59115.77995036666),
                 (59129.00074253334, 59129.78254386666),
                 (59154.00074043334, 59154.78506546666)]

    # get the urls for download
    data_urls = []
    # loop around the dates
    for it, date in enumerate(dates):
        print('Processing day {0} of {1}')
        # get the query keyword arguments
        qkwargs = dict(MINDATE=date[0], MAXDATE=date[1])
        # this may take some time
        cadc = Cadc()
        table = cadc.exec_sync(query.format(**qkwargs))
        # create a mask
        mask = np.ones(len(table), dtype=bool)
        # filter rows to the data we can download
        # loop around rows in table
        for row in range(len(table)):
            # only keep those we can download
            if len(table['DOWNLOADABLE'][row]) == 0:
                mask[row] = False
        # cut down the table
        get_table = table[mask]
        # break up the data into chunks to get urls
        chunks = 10
        for chunk in tqdm(range(len(get_table) // chunks)):
            start = chunk
            end = chunk + chunks
            with warnings.catch_warnings(record=True) as _:
                data_urls += cadc.get_data_urls(get_table[start:end])

    # return urls
    return data_urls


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":

    # -------------------------------------------------------------------------
    # get the dates (if we already have mini data)
    # get_dates('minidata1')
    # -------------------------------------------------------------------------
    # get data urls for a mini run
    urls = get_data('minidata1')
    # -------------------------------------------------------------------------
    # download data with wget or something else (untested)
    if DOWNLOAD:
        import wget
        for url in urls:
            wget.download(url)

# =============================================================================
# End of code
# =============================================================================
