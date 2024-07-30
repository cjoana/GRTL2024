import h5py
import numpy as np


"""
This script changes the 'time' in a given Chombo checkpointfile.
by C. Joana 
"""
###############

verbose = 2 
path = "./myfile.hdf5"   # include full path

new_time = 3.1415


###################


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


    num_levels = h5.attrs["num_levels"]

    print("Editing file...")
    # Modify the iteration attribute to zero.
    h5.attrs["time"] = new_time

    for lev in range(num_levels):
        h5[f'level_{lev}'].attrs["time"] = new_time


    print(f"  Done.") 

    
        
