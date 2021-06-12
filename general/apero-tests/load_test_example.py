"""
Example script to load a test object
​
Compatible with the tests-general branch.
Run from apero-utils/general/apero-tests for imports to work.
​
Before running the code, run this in shell:
​
apero-env
setup_mini1_06132
​
First execution might be very long because we parse all fits files on disk.
The dataframes are cached so subsequent executions should be faster.
"""
# Load the list of test definitions
# (takes some time because we need to parse all files on disk)
from apero_tests.spirou.test_definitions import tests

# Get a specific test (thermal here)
mytest = tests[-3]
print(mytest.name)

# Get log/index dataframe
log_df = mytest.log_df
ind_df = mytest.ind_df

# See the subtest.py file for example subtests that count log and output contentes
# Example of subtest manipulation with the Log count

n = len(mytest.subtest_list)

for i in range(n):
    log_subtest = mytest.subtest_list[i]
    print(log_subtest.description)
    log_subtest.run()  # Generate subtest result
    print(log_subtest.result)  # view subtest result
    print()


