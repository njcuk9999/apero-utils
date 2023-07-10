# APERO Checks

version: 0.0.12
date: 2023-07-06



# 1. Tasks to check NIRPS data at various points:

Task list here: [APERO NIRPS UdeM Checks Allocation](https://docs.google.com/spreadsheets/d/1s116aabnMH0zJ5YbXrGBWIVP17lmYXl6tYzLteoSEAc/edit#gid=0)

## Task descriptions

### 1.  raw data + red checks

- Go to here: [APERO NIRPS UdeM data monitor](https://docs.google.com/spreadsheets/d/1zvU_XFA1ZOJE111qZKiav7v6ptYveWpDkMHjdhiN06M/edit#gid=2059069236)
- Check for any "False" parameters on "RAW-NIRPS-HA" and "RED-NIRPS-HA"
- Check for any "False" parameters on "RAW-NIRPS-HE" and "RED-NIRPS-HE"
- If there are any "False" on raw data checks run [APERO raw data checks](https://github.com/njcuk9999/apero-drs/wiki/apero-tools#apero-raw-data-checks) for that column name. Email someone! Follow up!
- If there are any "False" on red data checsk run [APERO red checks](https://github.com/njcuk9999/apero-drs/wiki/apero-tools#apero-red-checks)
- Note you can run a single observation directory again if you think it might not be up-to-date
- Fill out the [NIRPS Data Check form](https://docs.google.com/forms/d/e/1FAIpQLScSeonrbNdVO7E_mXCW-Wkv0TgsM4pFqWda1LFkXYxoytXlDg/viewform) tick only what you have checked


###  2. online/offline reduction

- Run the [manual trigger](https://github.com/njcuk9999/apero-drs/wiki/apero-tools#apero-manual-trigger)

### 3. ARI

- Run the [APERO reduction interface] (https://github.com/njcuk9999/apero-drs/wiki/apero-tools#apero-reduction-interface-ari)
- This should be done after APERO reduction
- This should be done again after LBL 

### 4. LBL

- Run the [LBL](https://github.com/njcuk9999/apero-drs/wiki/apero-tools#lbl)

<br/><br/>
<br/><br/>

***

# 2. APERO raw data checks

These checks are designed to check the raw data for validity.
Results of running this code will be sent to: https://docs.google.com/spreadsheets/d/1zvU_XFA1ZOJE111qZKiav7v6ptYveWpDkMHjdhiN06M/edit#gid=2059069236

To run the raw data checks do the following:

    >> {apero conda env}
    >> apero_checks
    >> python apero_raw_data_checks.py {profile name}.yaml --obsdir={obs dir}

As a shortcut you can also use `--today` instead of `--obsdir={obs dir}` to run the tests on the most recent observation directory (assuming observation directories are named in the format `YYYY-MM-DD`).
e.g.

    >> {apero conda env}
    >> apero_checks
    >> python apero_raw_data_checks.py {profile name}.yaml --today

where:

- `{profile name}` is one of the yamls in the `apero_raw_data_check` directory
- `{obs dir}` is a single observation directory (e.g. 2023-01-01)
- `{apero conda env}` see [APERO setup for NIRPS](https://github.com/njcuk9999/apero-drs/wiki/nirps-general#apero-setup) or [APERO setup for SPIROU](https://github.com/njcuk9999/apero-drs/wiki/spirou-general#apero-setup)

Note one can leave off the `--obsdir={obs dir}` argument to get a prompt to run the tests for all nights (this will take longer).
Note only one observation directory can exist and the most recent will be used all others will be disregarded.
Note comments have to be added manually on a separate sheet (see the google sheet above).

## Running a single test (with logs on)

Similarly to above:

    >> {apero conda env}
    >> apero_checks
    >> python apero_raw_data_checks.py {profile name}.yaml --obsdir={obs dir} --test={test name}

where:

- `{test name}` is the name of the test (the column in the table)


## Adding new tests

All new tests must be added as follows:

1. make a branch from the developer branch of the apero_utils github: https://github.com/njcuk9999/apero-utils
2. Make you test in the same format as those in `apero-utils/nirps/apero_raw_data_check/apero_raw_tests/tests`
3. Copy your test into  `apero-utils/nirps/apero_raw_data_check/apero_raw_tests/tests` (make sure the test name is unique)
4. Update the `__init__.py` file in the `apero-utils/nirps/apero_raw_data_check/apero_raw_tests/tests` directory
5. git push to github and create a pull request, you will get feedback/comments or your test will be approved
6. Once your test is approved and you have verified it is in the master branch you must delete the branch you created for this test



<br/><br/>
<br/><br/>

***


# 3. APERO red checks


These checks are designed to check the APERO reduction data validity.
Results of running this code will be sent to: https://docs.google.com/spreadsheets/d/1zvU_XFA1ZOJE111qZKiav7v6ptYveWpDkMHjdhiN06M/edit#gid=2059069236

To run the reduced checks do the following:

    >> {apero conda env}
    >> apero_checks
    >> python apero_red_check.py {profile name}.yaml --obsdir={obs dir}

As a shortcut you can also use `--today` instead of `--obsdir={obs dir}` to run the tests on the most recent observation directory (assuming observation directories are named in the format `YYYY-MM-DD`).
e.g.

    >> {apero conda env}
    >> apero_checks
    >> python apero_red_check.py {profile name}.yaml --today

where:

- `{profile name}` is one of the yamls in the `apero_raw_data_check` directory
- `{obs dir}` is a single observation directory (e.g. 2023-01-01)
- `{apero conda env}` see [APERO setup for NIRPS](https://github.com/njcuk9999/apero-drs/wiki/nirps-general#apero-setup) or [APERO setup for SPIROU](https://github.com/njcuk9999/apero-drs/wiki/spirou-general#apero-setup)

Note one can leave off the `--obsdir={obs dir}` argument to get a prompt to run the tests for all nights (this will take longer).
Note only one observation directory can exist and the most recent will be used all others will be disregarded.
Note comments have to be added manually on a separate sheet (see the google sheet above).

## Running a single test (with logs on)

Similarly to above:

    >> {apero conda env}
    >> apero_checks
    >> python apero_red_check.py {profile name}.yaml --obsdir={obs dir} --test={test name}

where:

- `{test name}` is the name of the test (the column in the table)


## Adding new tests

All new tests must be added as follows:

1. make a branch from the developer branch of the apero_utils github: https://github.com/njcuk9999/apero-utils
2. Make you test in the same format as those in `apero-utils/nirps/apero_raw_data_check/apero_raw_tests/tests`
3. Copy your test into  `apero-utils/nirps/apero_raw_data_check/apero_raw_tests/tests` (make sure the test name is unique)
4. Update the `__init__.py` file in the `apero-utils/nirps/apero_raw_data_check/apero_raw_tests/tests` directory
5. git push to github and create a pull request, you will get feedback/comments or your test will be approved
6. Once your test is approved and you have verified it is in the master branch you must delete the branch you created for this test
