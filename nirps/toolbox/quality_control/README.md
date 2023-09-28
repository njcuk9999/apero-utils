# Quality control data working group 

author: Neil Cook
last update: 2023-07-31

Please edit this file when you add or change files in this directory.


## Index

- qc_over_time.py
    
        Utilities to plot apero QC checks over time for reduced files. Contains functions to:
            - plot all QCCs for a given file type
            - plot a QCC by name for a given file type and/or return a list of files that failed the QCC
            - plot a QCC by number (eg. '001') for a given file type.
        Alternatively, the script can be run directly to plot all QCCs for all file types specified in the
        files_to_check list.
        The destination path for figures can be modified at the beginning of the script.
