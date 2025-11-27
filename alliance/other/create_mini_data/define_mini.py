from astropy.table import Table, vstack
import numpy as np
import os
import shutil

# Define the CSV file containing NIRPS observation metadata
csv_file = 'nirps_mini_output.csv'

# Define input path where the raw NIRPS data is stored
input_path = '/project/6102120/apero/nirps_data/internal/nirps_raw/raw-data/nirps_he'

# Define output path where the mini dataset will be copied
output_path = '/project/6102120/apero/nirps_data/private/nirps_he_raw_minidata'

# Read the CSV file into an Astropy Table
tbl0 = Table.read(csv_file, format='csv')

# Filter out calibration files by removing rows where KW_OBJNAME is 'CALIB'
tbl = tbl0[tbl0['KW_OBJNAME'] != 'CALIB']

# Get unique observation directories (nights)
uobsdir = np.unique(tbl['OBS_DIR'])

# Initialize a new column to flag whether a star is a "hot star" (not PROXIMA or TOI4552)
tbl['HOTSTAR'] = True

# Mark PROXIMA and TOI4552 as NOT hot stars (False), all others remain True
tbl['HOTSTAR'] = np.where((tbl['KW_OBJNAME'] == 'PROXIMA') | (tbl['KW_OBJNAME'] == 'TOI4552'), False, True)

# Counter for good observation directories
ngood = 0

# List to store unique hot star names
name_hot = []

# Initialize the trimmed table
tbl_trimmed = None

# First pass: Find observation directories with multiple observations of PROXIMA, TOI4552, and hot stars
for u in uobsdir:
    # Create a mask for the current observation directory
    mask = tbl['OBS_DIR'] == u
    sub_tbl = tbl[mask]
   
    # Count observations of PROXIMA in this directory
    n_proxima = np.sum(sub_tbl['KW_OBJNAME'] == 'PROXIMA')
    
    # Count observations of TOI4552 in this directory
    n_toi4552 = np.sum(sub_tbl['KW_OBJNAME'] == 'TOI4552')

    # Count observations of hot stars in this directory
    n_hot = np.sum(sub_tbl['HOTSTAR'])

    # Only keep directories with >1 observation of each type
    if n_proxima > 1 and n_toi4552 > 1 and n_hot > 1:
        # Stack this night's data into the trimmed table
        if tbl_trimmed is None:
            tbl_trimmed = sub_tbl
        else:
            tbl_trimmed = vstack([tbl_trimmed, sub_tbl])

        # Print summary for this observation directory
        print(f'OBS_DIR: {u}, PROXIMA count: {n_proxima}, TOI4552 count: {n_toi4552}, HOTSTAR count: {n_hot}')
        print(sub_tbl)
        ngood += 1
        
        # Collect unique hot star names from this night
        name_hot_sub = np.array(sub_tbl['KW_OBJNAME'][sub_tbl['HOTSTAR']])
        for name in name_hot_sub:
            if name not in name_hot:
                name_hot.append(name)

# Get unique object names from the trimmed table
uobj = np.unique(tbl_trimmed['KW_OBJNAME'])

# Array to store the number of nights each object was observed
n_nights_obj = np.zeros(len(uobj), dtype=int)

# Count how many different nights each object was observed on
for i, obj in enumerate(uobj):
    # Mask for this specific object
    mask_obj = tbl_trimmed['KW_OBJNAME'] == obj
    
    # Find unique nights where this object was observed
    nights_obj = np.unique(tbl_trimmed['OBS_DIR'][mask_obj])
    
    # Store the count
    n_nights_obj[i] = len(nights_obj)

# Second pass: Remove objects that were only observed on a single night
final_tbl = None
for i, obj in enumerate(uobj):
    # Only keep objects observed on more than one night
    if n_nights_obj[i] > 1:
        mask_obj = tbl_trimmed['KW_OBJNAME'] == obj
        sub_tbl = tbl_trimmed[mask_obj]
        
        # Stack into final table
        if final_tbl is None:
            final_tbl = sub_tbl 
        else:
            final_tbl = vstack([final_tbl, sub_tbl])

# Third pass: Remove any nights that don't have at least one hot star
final_uobsdir = np.unique(final_tbl['OBS_DIR'])
cleaned_tbl = None
for u in final_uobsdir:
    # Mask for this observation directory
    mask = final_tbl['OBS_DIR'] == u
    sub_tbl = final_tbl[mask]
   
    # Count hot stars in this night
    n_hot = np.sum(sub_tbl['HOTSTAR'])

    # Only keep nights with at least one hot star
    if n_hot > 0:
        if cleaned_tbl is None:
            cleaned_tbl = sub_tbl
        else:
            cleaned_tbl = vstack([cleaned_tbl, sub_tbl])

# Fourth pass: Reorder the table by OBS_DIR using the original full table
# This ensures we get all files (including calibrations) for selected nights
final_uobsdir = np.unique(cleaned_tbl['OBS_DIR'])
ordered_tbl = None  
for u in final_uobsdir:
    # Use original table (tbl0) to get ALL files for this night, including calibrations
    mask = tbl0['OBS_DIR'] == u
    sub_tbl = tbl0[mask]
   
    # Stack into ordered table
    if ordered_tbl is None:
        ordered_tbl = sub_tbl
    else:
        ordered_tbl = vstack([ordered_tbl, sub_tbl])  

# Print summary of each selected night
for unight in np.unique(ordered_tbl['OBS_DIR']):
    g = ordered_tbl['OBS_DIR'] == unight  
    print(' -- Night: ', unight)
    print(ordered_tbl[g])

# Copy all files from selected nights to the output directory
for i in range(len(ordered_tbl)):
    # Get filename and observation directory for this entry
    fname = ordered_tbl['FILENAME'][i]
    OBS_DIR = ordered_tbl['OBS_DIR'][i]
    full_path = os.path.join(output_path, OBS_DIR)
    if not os.path.exists(full_path):
        print(f'Creating directory: {full_path} 🛠️')
        os.makedirs(full_path, exist_ok=True)

    # Construct full input and output file paths
    input_file = os.path.join(input_path, OBS_DIR, fname)
    output_file = os.path.join(output_path, OBS_DIR, fname)

    # check that input file exists and skip if not
    if not os.path.exists(input_file):
        print(f'Input file does not exist, skipping: {input_file} ❌')
        continue
    # check that the output file does not already exist
    if os.path.exists(output_file):
        print(f'Output file already exists, skipping: {output_file} ⚠️')
        continue
    # Print copy operation (with fun emojis for Neil and Lison!)
    print(f'Copying {input_file} to {output_file} 🚀✨')
    
    # Actually perform the copy (currently commented out)
    shutil.copyfile(input_file, output_file)
