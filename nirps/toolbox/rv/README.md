# Radial velocity data working group 

author: Neil Cook
last update: 2023-03-13

Please edit this file when you add or change files in this directory.


## Index

- multi_object_validations.py
    
        Plots the RV rms vs median RV uncertainty and BERV coverage for all HA and HE lbl_target_template.rdb files. 
        The output is two figures (HA_multi-object_validation.pdf and HE_multi-object_validation.pdf)
        for which the top row shows all the points and the bottom row has the > 1 sigma outliers removed.

- [compare_apero_eso.py](https://github.com/njcuk9999/apero-utils/blob/olim_apero-eso-rv/nirps/toolbox/rv/compare_apero_eso.py)

Make plots to compare APERO and ESO radial velocities computed from LBL or CCF. 
The user should set up the variables in the "Define variables" section.

Special requirements: [dace_query](https://dace.unige.ch/pythonAPI/dace-query-doc/1.1.0/html/dace_introduction.html#installation)

To run the code from a terminal:
```
python compare_apero_eso.py
```
