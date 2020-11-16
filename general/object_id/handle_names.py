"""
Compare header names and CFHT
"""
import pandas as pd
import utils as ut

path_to_data = '/home/vandal/Documents/spirou/obj_db/object_info.csv'
sheet_id = '1jwlux8AJjBMMVrbg6LszJIpFJrk6alhbT5HA7BiAHD8'

# Load all data
df_sheet = ut.get_full_sheet(sheet_id)
df_local = pd.read_csv(path_to_data)

# For each name, see if in OBJECT