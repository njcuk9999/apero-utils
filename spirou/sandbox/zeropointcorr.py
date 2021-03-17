import os

import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d
from tqdm import tqdm


def zpcorr(zeropoint_file, obj_sci, obj_template=None, force=True):
    """
    Function to correct the lbl rv curve with the latest Gaussian Process Zero
    point.
    """

    # zeropoint_file format
    # zeropoint_vx.csv
    version = zeropoint_file.split("_")[1]

    if obj_template is None:
        obj_template = obj_sci

    # input table
    inname = "lbl_{0}_{1}.rdb".format(obj_sci, obj_template)
    # output table
    outname = "lbl_{0}_{1}_zp{2}.rdb".format(obj_sci, obj_template, version)

    if (not os.path.isfile(outname)) or force:
        tbl = Table.read(inname, format="rdb")

        bjd, rv, erv = tbl["BJD"], tbl["vrad"], tbl["svrad"]
        bjd_zp, gp_zp, gp_ezp = np.genfromtxt(
            "{0}.csv".format(zeropoint_file), delimiter=",", unpack=True
        )

        # Zero-point correction (interpolation at tbl['BJD'] epochs)
        fint = interp1d(bjd_zp, gp_zp, fill_value="extrapolate")
        bjd_zp_resamp = bjd.copy()
        gp_zp_resamp = fint(bjd_zp_resamp)

        fint = interp1d(bjd_zp, gp_ezp, fill_value="extrapolate")
        bjd_zp_resamp = bjd.copy()
        gp_ezp_resamp = fint(bjd_zp_resamp)

        rv -= gp_zp_resamp
        erv = np.sqrt(erv ** 2 + gp_ezp_resamp ** 2)

        tbl_zp = tbl.copy()
        tbl_zp["vrad"] = rv
        tbl_zp["svrad"] = erv

        # Write to output
        tbl_zp.write(outname, overwrite=True)

        # output binned table
        outname2 = "lbl2_{0}_{1}_zp{2}.rdb".format(obj_sci, obj_template, version)

        udates = np.unique(tbl["DATE-OBS"])

        tbl_zp2 = Table(tbl[0 : len(udates)])  # create a table with a per-epoch value

        for i in tqdm(range(len(udates))):
            tbl_date = tbl_zp[udates[i] == tbl_zp["DATE-OBS"]]
            for key in tbl_date.keys():
                if "vrad" not in key:
                    try:
                        tbl_zp2[key][i] = np.mean(tbl_date[key])
                    except:
                        tbl_zp2[key][i] = tbl_date[key][0]

                if key[0:4] == "vrad":
                    rv = tbl_date[key]

                    err_rv = tbl_date["s" + key]
                    tbl_zp2[key][i] = np.nansum(rv / err_rv ** 2) / np.nansum(
                        1 / err_rv ** 2
                    )

                    tbl_zp2["s" + key][i] = np.sqrt(1 / np.nansum(1 / err_rv ** 2))
        tbl_zp2.write(outname2, overwrite=True)
