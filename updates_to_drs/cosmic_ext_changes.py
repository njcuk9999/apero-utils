# test on file 2503421o_pp_e2dsff_AB.fits

# insert into apero/science/extract/extraction.py
#   function = extraction(...)

# first guess of extracted flux
guess = np.nanmedian(sx/fx)
residual = sx - fx * guess

nsig_dev = residual / robust_sigma(residual)
bad = np.abs(nsig_dev) > cut
weights[bad] = 0



