# APERO raw data checks

These checks are designed to check the raw data for validity.
Results of running this code will be sent to: https://docs.google.com/spreadsheets/d/1zvU_XFA1ZOJE111qZKiav7v6ptYveWpDkMHjdhiN06M/edit?usp=sharing


To run the raw data checks do the following:

    >> {apero conda env}
    >> apero_raw_data_checks
    >> python apero_raw_data_checks.py {profile name}.yaml --obsdir={obs dir}


where:

- `{profile name}` is one of the yamls in the `apero_raw_data_check` directory
- `{obs dir}` is a single observation directory (e.g. 2023-01-01)

Note one can leave off the `--obsdir={obs dir}` argument to get a prompt to run the tests for all nights (this will take longer).
Note only one observation directory can exist and the most recent will be used all others will be disregarded.

# Adding new tests

All new tests must be added as follows:

1. make a branch from the master branch of the apero_utils github: https://github.com/njcuk9999/apero-utils
2. Make you test in the same format as those in `apero-utils/nirps/apero_raw_data_check/apero_raw_tests/tests`
3. Copy your test into  `apero-utils/nirps/apero_raw_data_check/apero_raw_tests/tests` (make sure the test name is unique)
4. Update the `__init__.py` file in the `apero-utils/nirps/apero_raw_data_check/apero_raw_tests/tests` directory
5. git push to github and create a pull request, you will get feedback/comments or your test will be approved
6. Once your test is approved and you have verified it is in the master branch you must delete the branch you created for this test

