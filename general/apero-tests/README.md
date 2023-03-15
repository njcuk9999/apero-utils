# APERO Tests

*Note: This branch is a work in progress to make the tests compatible with
APERO version >= 0.7. Use tests-v06 if you want to add new tests or run the
tests with APERO 0.6 outputs*

Series of tests to make sure everything was properly reduced by APERO.

## Development and running the tests
The tests are meant to be integrated in APERO and are not a package on their own (i.e. are developed in apero-utils), but for development purposes, `apero_tests` is treated as a standalone package.

Make sure they are in your `PYTHONPATH` (with conda, this is just `conda develop /path/to/apero-tests`, or adding a `.pth` file in your virtual environment's site-package directory).
Or even simpler, run everyting from the parent directory.

## Structure
- The `main.py` script runs tests one by one and stores the output in an HTML report.
- HTML templates for jinja2 are in `templates`.
- The code for all tests is in the `tests` directory
  - A generic `DrsTest` class is in `apero_tests.drs_test`. It can generate its
    attribute based on a recipe, but this is not required.
  - The `DrsTest` class relies on `SubTest` objects to run the actual checks
    required here. The `SubTest` class is also flexible. You can either create
    a children class with a `.run()` method to run the test, or simply set the
    required attributes (`.result` for example) manually.
  - Tests for each SPIROU recipe are in `apero_tests.spirou.test_definitions`.
  - Utility functions are defined in `utils.py`.

### Creating a new tests
A new test can be created by creating an instance of `DrsTest`. Then, either
give it a recipe to generate its tests from or use it to store custom subtests.
There are many examples of subtests in `apero_tests.{global_subtests,subtest}`
and there are examples of `DrsTest` objects in
`apero_tests.spirou.test_definitions`.


## Outputs
A brief summary is generated in `summary.html` (in the output directory, `out` by default).
A more detailed overview of each test is accessible from there, including links to interactive
plots and tables.
