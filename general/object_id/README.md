# APERO Object Identification

These scripts facilitate the management of the object identification database for APERO.

### Verifying Object ID
The _check_and_update.py_ script is used to verify the Gaia ID of an object and to obtain the corresponding 2MASS ID. The steps are:
- String manipulation to cross-match the object name and the aliases from SIMBAD. This step finds most of the objects
- Cross-checking with TIC to verify the Gaia ID of TOIs, which are often not explicitly named in SIMBAD. Only objects not found previously are checked.
- Searching for the object name in SIMBAD to see if the Gaia ID matches the one we have. Only object still not found are checked.

To obtain the 2MASS ID, we:
- Use the Gaia catalog to cross-match with 2MASS. This finds most objects.
- We then search the SIMBAD aliases in the database for a 2MASS ID, only for objects not found above.

### Updating RV and Teff

To do.
