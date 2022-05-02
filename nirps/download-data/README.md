# NIRPS Data Download

Scripts to get NIRPS data from the eso archvie.

- `eso-tools`: scripts and notebooks from [ESO Science Archive Programmatic and Tools Access](http://archive.eso.org/programmatic/#SCRIPT).

**Note about ESO interface**: `astroquery.eso` uses the http requests which will
be depcreated in 2022 according to [this
comment](https://github.com/astropy/astroquery/issues/1818#issuecomment-1031372421).
There also seem to be some bugs present with the http interface that are not
present in the TAP interface. We use the TAP interface directly for now. We can
switch to `astroquery.eso` once they will have added support for TAP.
