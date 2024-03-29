# APERO Object ID Management
The [APERO Object ID Google Sheet](https://docs.google.com/spreadsheets/d/1jwlux8AJjBMMVrbg6LszJIpFJrk6alhbT5HA7BiAHD8/edit#gid=0) stores information about objects analyzed with APERO. These
files facilitate the management of the sheet with two executable scripts and various utility
functions.

See the [Wiki](https://github.com/njcuk9999/apero-utils/wiki) for information on [APERO object identification](https://github.com/njcuk9999/apero-utils/wiki/APERO-Object-Identification) and [using Google sheets with Python](https://github.com/njcuk9999/apero-utils/wiki/Working-with-Google-Sheets).

## Files
- `cleandfits` is a short bash script to pre-format the dfits output as astropy table tab-separated table (stripping spaces around columns)
  - Simply run it with the files as arguments (e.g. `cleandfits *.txt`). Files will be backed up by adding the `.bak format`
  - Not necesseary, the code should handle dfits outputs. This is just provided as a utility script for convenience.
- `full_update.py`: Performs a full update of the sheet using local data. This performs the following tasks:
  - Fetch data from `.fits` headers or use a pre-constructed local `.csv` file.
  - Include new objects in the main sheet.
  - Add aliases for existing objects.
  - Automatically flag files that require attention in the [Need Sorting](https://docs.google.com/spreadsheets/d/1jwlux8AJjBMMVrbg6LszJIpFJrk6alhbT5HA7BiAHD8/edit#gid=1396997839) tab.
  - Update information in the [Info](https://docs.google.com/spreadsheets/d/1jwlux8AJjBMMVrbg6LszJIpFJrk6alhbT5HA7BiAHD8/edit#gid=0) and [Maintenance](https://docs.google.com/spreadsheets/d/1jwlux8AJjBMMVrbg6LszJIpFJrk6alhbT5HA7BiAHD8/edit#gid=247392185) tabs.
- `remote_update.py`: Updates the information in the remote sheet without adding new objects.
  - Main purpose: automatically move files marked CHECKED from Need Sorting to Info
  - Also updates the information in Info and Maintenance
- `astro_object`: Used to cross-match positions with Gaia.
- `maintenance.py`: Main building blocks of the two scripts.
- `utils.py`: Utility functions for smaller tasks, used by maintenance and the scripts.
