# Getting started with the CCF codes


Change the line at the top of _ccf2rf.py_ to point toward the base
folder where you will put everything. The folder will contain a sub-folder 
that may be called _all_ccfs_ (you can name it something else by giving a 
parameter to the function _dispatch_object_).

Download all the CCFs from the DRS. These files will have names will look like _2434981o_pp_e2dsff_tcorr_AB_ccf_masque_sept18_andres_trans50_AB.fits_ and 
they should all be in the appropriate sub-folder (unless you want to tempt devil, they should be in the sub-folder
_all_ccfs_). 

First, dispatch all CCF files into per-object folders.
- Start python
- from ccf2rv import *
- Assuming that you want to create folders for **all** objects ever, just run 
_dispatch_object('all')_ (maybe with _all_ccf_dir = 'my_devil_folder'_ as a parameter)
- Maybe you just want one object, then run _dispatch_object('TOI-1278')_

You are (or should) be all set to get some RV measurements!

Try getting an RV for a random object:

- _get_object_rv('TOI-1278')_

You may force 

You can measure RVs with a number of 'methods' (2 are currently defined):

- A gaussian fit to the CCF (parameter: method = 'gaussian')
- Mid-point of bisector between two heights. Let's say you want between the 30th and 70th 
percentile of the CCF, that would be method = 'bisector_30_70'

... you are invited to contribute other methods!

If you obtained CCFs with a number of masks, then you need to specify a mask name as 
an input parameter to the _ccf2rv_ function: mask = 'my_clever_mask', otherwise it uses the DRS default.

The code determines the proper order weights, but this can be forced 
by provinding a CSV file as an argument (_weight_table = "my_weights.csv"_). The 
program generates a CSV file when running in the default 'smart' mode where 
the weights are derived by the code. You can always use that file if you want to
have the same weights for another target.

Outputs : 

Lot's of plots and a CSV file that are named OBJECT_mask_MASKNAME_METHOD.csv. This file
contains tons of info about the data the velocity of your object.
