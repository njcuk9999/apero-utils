# %%
from astropy.io import fits
import pprint

FILE1 = "/nirps_raw/nirps/raw-data/nirps_he/2023-10-29/NIRPS_2023-10-30T00_15_11_266.fits"
FILE2 = "/nirps_raw/nirps/raw-data/nirps_ha/2023-10-29/NIRPS_2023-10-30T00_02_55_544.fits"

hdr1 = fits.getheader(FILE1)
hdr2 = fits.getheader(FILE2)

# %%
missing_from_2 = []
for key in hdr1:
    if key not in hdr2:
        missing_from_2.append(key)
missing_from_1 = []
for key in hdr2:
    if key not in hdr1:
        missing_from_1.append(key)

# %%
print(f"File 1: {FILE1}")
print(f"File 2: {FILE2}")
# %%
print("Keys in 2 not in 1:")
pprint.pprint(missing_from_1)

# %%

print("Keys in 1 not in 2:")
pprint.pprint(missing_from_2)
