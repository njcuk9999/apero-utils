Montreal APERO data sync (apero_mtl_sync.py)

===============================================
Options are:
===============================================

--name={NAME}       Change the directory name (i.e. mini_data_0_6_037)

--exclude           Night directories to exclude (separated by commas NO SPACES)
                    If spaces required use strings i.e.
                    'night 1, night 2, night 3'

--include           Night directories to include (separated by commas NO SPACES)
                    If spaces required use strings i.e.
                    'night 1, night 2, night 3'
                    Note wild cards may be used

--suffix            Define a set of characters that must be in any output
                    (can be used multiple times in one command but logic is as
                    follows: 
                        suffix1 OR suffix2 OR suffix3
                    Note: this complements individual recipes (i.e. it
                    does not filter outputs for a recipe)

--test              Perform a test (dry-run) without any actual copying
                    Recommended the first time this is run
                    Note wild cards may be used
                    
--debug             Performs a debug just printing the rsync command that will
                    be run when not in debug mode

--instrument        Change the instrument to download for (Default: SPIROU)
                    Must be one of the following SPIROU, NIRPS_HA

--help / -h         Display this help message

===============================================
Arguments are one or any of the following:
===============================================

Full data sets (warning large amount of data)

RAW                 Download all raw data
TMP                 Download all tmp data (preprocessed)
REDUCED             Download all reduced data 
CALIBDB             Download all calibration data (from calibDB)
TELLUDB             Download all telluric data (from telluDB)

Individual recipes (Warning these require APERO to be installed to use)

PP                  Download outputs for cal_preprocessing
BAD                 Download outputs for cal_badpix
DARK                Download outputs for cal_dark
DARKM               Download outputs for cal_dark_master
LOC                 Download outputs for cal_loc
SHAPEM              Download outputs for cal_shape_master
SHAPE               Download outputs for cal_shape
FF                  Download outputs for cal_ff
THERM               Download outputs for cal_thermal
EXT                 Download outputs for cal_extract
WAVE                Download outputs for cal_wave
WAVEM               Download outputs for cal_wave_master
CCF                 Download outputs for cal_ccf
MKTELL              Download outputs for obj_mk_tellu
FTELLU              Download outputs for obj_fit_tellu
MKTEMP              Download outputs for obj_mk_template
POLAR               Download outputs for pol

===============================================
Examples:
===============================================

The following example would download all calibration and telluric files
    for mini_data_0_6_037  (in debug mode) - prints only rsync command
    
>> apero_mtl_sync.py --debug --name=mini_data_0_6_037 CALIBDB TELLUDB

--------------------------------------------------------------------------------

The following example would download all extracted and ccf outputs
    for mini_data_0_6_037 in night directories 2019-04-20 and 2019-04-19
    in test mode (runs rsync but does not copy files)
    
>> apero_mtl_sync.py --test --name=mini_data_0_6_037 --include='2019-04-20,2019-04-19' EXT CCF 

--------------------------------------------------------------------------------

The following example would download badpix outputs
    for mini_data_0_6_037 except for night directories 2019-02-20 and 2019-02-19
    in test mode (runs rsync but does not copy files)
    
>> apero_mtl_sync.py --test --name=mini_data_0_6_037 --exclude='2019-04-20,2019-04-19' BAD 

--------------------------------------------------------------------------------

The following example would download all raw data for 2019-04-20

>> apero_mtl_sync.py --test --name=mini_data_0_6_037 --include=2019-04-20 RAW

--------------------------------------------------------------------------------

The following example would download all e2ds_AB files from the reduced directory
    for mini_data_0_6_037 in the night directory 2019-04-20
    in test mode (runs rsync but does not copy files)

>> apero_mtl_sync.py --test --name=mini_data_0_6_037 --suffix=e2ds_AB --include=2019-04-20 REDUCED

--------------------------------------------------------------------------------
