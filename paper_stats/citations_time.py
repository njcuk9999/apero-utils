#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-05-13 at 14:05

@author: cook
"""
import glob
import os
from collections import Counter
from datetime import datetime
from time import sleep
from typing import Any, Dict, List

import ads
import yaml
import matplotlib.pyplot as plt
from astropy.time import Time
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
if 'ADS_TOKEN' in os.environ:
    ads.config.token = os.environ['ADS_TOKEN']
else:
    question = ('Get token from here: '
                'https://ui.adsabs.harvard.edu/user/settings/token'
                '\n\tEnter token:\t')

    ads.config.token = input(question)
    # set environment
    os.environ['ADS_TOKEN'] = ads.config.token

# -----------------------------------------------------------------------------

BIBCODES = dict()
BIBCODES['APERO2022'] = '2022PASP..134k4509C'
BIBCODES['LBL2022'] = '2022AJ....164...84A'


# =============================================================================
# Define functions
# =============================================================================
def pubdate_to_decimal_year(pubdate_str):
    try:
        # Handle YYYY-MM-00 by replacing with mid-month
        if pubdate_str.endswith("-00"):
            pubdate_str = pubdate_str[:-3] + "-15"

        date = datetime.strptime(pubdate_str, "%Y-%m-%d")
    except Exception:
        return None
    year_start = datetime(date.year, 1, 1)
    year_end = datetime(date.year + 1, 1, 1)
    return round(date.year + (date - year_start).total_seconds() /
                 (year_end - year_start).total_seconds(), 3)


def get_citations(target_bibcode):
    # Get list of bibcodes that cite this paper
    citations = list(ads.SearchQuery(
        q=f'citations(bibcode:{target_bibcode})',
        fl=['bibcode', 'pubdate'],  # fields you want
        rows=10000  # increase if needed
    ))

    return citations

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # loop around papers and get references
    citations: Dict[str, List[float]] = dict()

    for bibcode in BIBCODES:
        # Fetch the paper
        clist = get_citations(BIBCODES[bibcode])

        if len(clist) > 0:
            # loop around references and get the year of each paper
            citations[bibcode] = []
            for c_obj in tqdm(clist):

                pubdate = getattr(c_obj, 'pubdate', None)

                if pubdate is None:
                    dec_year = c_obj.year + 0.5
                else:
                    dec_year = pubdate_to_decimal_year(pubdate)


                if dec_year is None:
                    dec_year = int(c_obj.year) + 0.5

                citations[bibcode].append(dec_year)

        else:
            print("Paper not found.")


    plt.close()
    fig, frame = plt.subplots()

    for bibcode in BIBCODES:
        # Convert decimal years to (year, month)
        year_month = [(int(y), int((y % 1) * 12) + 1) for y in citations[bibcode]]

        # Count occurrences per (year, month)
        counts = Counter(year_month)

        # Sort and generate cumulative counts
        sorted_keys = sorted(counts)
        months = []
        ccounts = []
        decimal = []
        total = 0

        for ym in sorted_keys:
            total += counts[ym]
            months.append(f"{ym[0]}-{ym[1]:02d}")
            ccounts.append(total)
            decimal.append(ym[0] + ym[1]/12)

        dates = Time(decimal, format='decimalyear')

        # Plot
        frame.plot(dates.plot_date, ccounts, label=bibcode,
                   marker='o')

        frame.xaxis.axis_date()
        frame.tick_params(axis='x', rotation=45)
        frame.set_xlabel("Month")
        frame.set_ylabel("Cumulative Citations")
        frame.set_title("Cumulative Citations per Month")

    plt.legend(loc=0)
    plt.tight_layout()
    plt.grid(True)
    plt.show(block=True)

# =============================================================================
# End of code
# =============================================================================
