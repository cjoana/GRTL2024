import h5py
import numpy as np


#### Instructions: 
#
#  Select your checkpoint-file to reinitialise, e.g. run0_999000.3d.hdf5
#  copy it with a new name as set in 'path'. e.g.  run0_rerun.3d.hdf5
#  run this script to rewrite the iteration attribute to 0. 
#
#  You are ready to rerun the file (run0_rerun.3d.hdf5), and GRChombo
#  will begin from iteration 0,   i.e. run0_000000.3d.hdf5


verbose = 2 
path = "./run0_rerun.3d.hdf5"

with h5py.File(path,  "a") as h5:
	
	# Basic print inspection.
	if verbose >1: 
		print("Print of attributes:")
		h5_attrs = h5.attrs.keys()
		print(f" > base atts are {h5_attrs}")
		print("")
		
		h_keys = h5.keys()
		for k in h_keys:
			sbh5 = h5[k]
			sbh5_att =sbh5.attrs.keys()
			print(" >> atts in ", k, " are ", sbh5_att)
		print("")


	print("Editing file...")
	# Modify the iteration attribute to zero.
	h5.attrs["iteration"] = 0
	
	print(f"    edited {path}") 
	print(f"    >> MODIFIED:  base-attr. 'iteration' to {h5.attrs['iteration']}")
        
