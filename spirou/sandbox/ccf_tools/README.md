# Getting started with the CCF codes

First, download all the CCFs from the DRS. These files will have names will look like _2434981o_pp_e2dsff_tcorr_AB_ccf_masque_sept18_andres_trans50_AB.fits_ and 
they should all be in the appropriate sub-folder (unless you want to tempt devil, they should be in the sub-folder
_all_ccfs_). 

First, dispatch all CCF files into per-object folders.
- Start python
- from ccf2rv import *
- Assuming that you want to create folders for **all** objects ever, just run 
_dispatch_object('all')_ (maybe with _in_dir = 'my_devil_folder'_ as a parameter)
- Maybe you just want one object, then run _dispatch_object('OBJECT-NAME')_

You are (or should) be all set to get some RV measurements!

Try getting an RV for a random object:

- _get_object_rv('TOI-1278', mask='MASK_NAME')_

You may force 

You can measure RVs with a number of 'methods' (3 are currently defined):

- A Gaussian fit to the CCF (parameter: method = 'gaussian')
- Mid-point of bisector between two heights. Let's say you want between the 30th and 70th 
percentile of the CCF, that would be method = 'bisector_30_70'

... you are invited to contribute other methods!

If you obtained CCFs with a number of masks, then you need to specify a mask name as 
an input parameter to the _ccf2rv_ function: mask = 'my_clever_mask', otherwise it uses the DRS default.

An alternative way to use _get_object_rv_ is to input a list of CCF file paths and the corresponding run ID (parameters: _ccf_files_, _run_id_). Run ID must have the format ("OBJECT__MASK__SANIT__DRSVERSION"). This is mainly useful for automated scripts that automatically run multiple batches of files.

The code determines the proper order weights, but this can be forced 
by providing a CSV file as an argument (parameter: _weight_table = "my_weights.csv"_). The 
program generates a CSV file when running in the default 'smart' mode where 
the weights are derived by the code. You can always use that file if you want to
have the same weights for another target.

Outputs : 

Lot's of plots and a CSV file that are named OBJECT_mask_MASKNAME_METHOD.csv.
This file contains tons of info about the data the velocity of your object. If
_outdir_ is specified, everything is saved there. Otherwise, the parent
directory of CCF files is used.



# All you ever wanted to know about get_object_rv but were too afraid to ask

### obj 
The name of your target, this tells the code in which folder to look for CCF files.
### mask 
The mask name used in your CCFs. This is used to search for files in the folder. For a file like _TOI-1278/2493616o_pp_e2dsff_tcorr_AB_ccf_gl846_neg_AB.fits_
you would enter _mask = "gl846_neg"_.
### sanitize = False
Look for CCFs of sanitized files with a name line _TOI-1452/2499239o_pp_e2dsff_sani_AB_ccf_gl846_neg_AB.fits_,
with a _sani_ name. This indicates that the sanitizing code has been passed on the data. Otherwise, the 
code expect a _tcorr_ in the name in place of the _sani_.
### ccf_parent
Prent directory of all CCF files. This directory should contain another directory called OBJECT-NAME where CCF files are stored.
### outdir
Directory where output should be saved.
### ccf_files
List of CCF files for  a given run. If this is provided, _run_id_ must be given to.
### run_id
ID of the run with format OBJECT__MASK__SANIT__DRSVERSION. If this is provided, _ccf_files_ must be given to. This overrides obj, mask, and sanitize arguments.
### method = 'template'
Method to compute the velocity. For now, it is better to stick to the method _template_ as it is the one
that has been tested more. We'll try to get _gaussian_ and _bisector_ up to speed at some point. If
you want to use the bisector method, then you need to specify a depth. The position of the lines 
is determined between a min and a max depth. To get the mean bisector value between the 30th and 70th 
percentile, you would use _method = 'bisector_30_70'_... of course, don't use bounds below ~5 or 
above ~95 or things will break!
### exclude_orders = [-1]
You can exclude orders. Not that the code may exclude other orders that have mis-behaving CCFs. You'll be
warned if other orders are rejected.
### weight_table = ''
You can provide a weight table. For now this is still experimental.
### force = True
If set to _False_ then if the CCF table exists, it is simply read and returned. This can be nice 
if you always  call _get_object_rv_ at the start of a program but just want the table returned.
### snr_min = 0.0
Set a threshold below which CCFs are rejected. We use the extracted SNR for order 35 as a threshold. 
This should be first set to 0 first, then you look at the plot of SNRs as a function of time for your
favorite target, then pick a reasonable threshold.
### weight_type = ''
Still exprerimental, leave this alone!
### bandpass = 'YJHK'
Which photometric bandpasses will be looked at. You can set 'YK' if you just want to get a mean RV for these
two bands (though that would be a weird pick). If you know that the target is red and has very weak Y and J features,
using bandpass = 'HK' would make sense.

Definition of photometric bandpasses (see [values here](http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=CFHT&gname2=Wircam&asttype=)):
* Y = 938.600	1113.400
* J = 1153.586	1354.422
* H = 1462.897	1808.544
* K = 1957.792	2400 (modified to get to the edge of the SPIRou domain)

### velocity_window = 10
Window (in km/s) used around the CCF minimum used to measure velocity. We go from -window to +window in km/s. 
When adjusting the median CCF to the CCF of one observation (in CCF _template_), this is the width over 
which the dot product between the derivative and the residuals is computed. When determining the Q value
of the CCF per band, this is also the width of the integral.
### dvmax_per_order = 1.0
Reject any order than has a CCF minimum that is beyond this dvmas (km/s). For very low SNR CCF in some ordres,
the median CCF of the order is not a (near) gaussian line, and the minimum can be way off. This is used 
to reject these orders. You get a warning when this happens.
### doplot = True
Set to _False_ if you don't want to see all the nice plots I prepared for you. Can be useful if you
want to put the code in a big batch and don't want to spend your weekends closing python windows.
### do_blacklist = True
Check if the files you input are blacklisted. Default is _False_ as this 
adds some overheads.
### verbose = False
Print some info and debugging messages.
