# APERO tests

Series of tests to make sure everything was properly reduced by APERO.

## Structure
- The `main.py` script runs tests one by one and stores the output.
- HTML templates created with jinja2 are in `templates`.
- The code for all tests is in the `tests` directory
  - Abstract classes are stored in `tests.py`. A general `Test` class defines
    properties and methods common to all tests and abstract subclasses for types of tests
    (e.g. calibration or science) are also defined.
  - There is a file for each recipe with the test objects corresponding to this recipe.
  - Utility functions are defined in `utils.py`.
  - A factory function is stored in `factory.py` to facilitate test creation in the main script.
  - Output and template paths are defined in the `__init__.py` file.


## Outputs
A brief summary is generated in `summary.html` (in the output directory, `out` by default).
A more detailed overview of each test is accessible from there, including links to interactive
plots and tables.

## List of Tests
1. Preprocessing Recipe Test #1
2. Dark Master Recipe Test #1
3. Bad Pixel Correction Recipe Test #1
4. Localisation Recipe Test #1
5. Shape Master Recipe Test #1
6. Shape (per night) Recipe Test #1
7. Flat/Blaze Correction Recipe Test #1
8. Thermal Correction Test #1
9. Master Leak Correction Recipe Test #1
10. Leak (per night) Correction Test #1
11. Master Wavelength Solution Recipe Test #1
12. Wavelength Solution (per night) Test #1
13. Extraction Recipe Test #1
14. Extraction Recipe Test #2
15. Extraction Recipe Test #3
16. Make Telluric Recipe Test #1
17. Fit Telluric Recipe Test #1
18. Make Template Recipe Test #1
19. CCF Recipe Test #1
