"""
This script generates a Chombo checkpoint file that can be used 
to start GRChombo simulations. The specificities of the initial data 
(grid params, initial values, etc) are set up in a separated
params python-file, e.g. write_params_[myscience].py  imported below.

Please note that current script is restricted to seetings with: 
	* Periodic Boundary Conditions 
	* With no ghosts-cells. 

author: Cristian Joana
"""

################## Settings ########################
verbose = 1
overwrite = True
filename = "./test.hdf5"

# import params file (problem specific).
from write_params_1SF_minimal import *

#####################################################


############  Make Checkpoint file ##################

import numpy as np
import h5py as h5
import os


if overwrite:
    if os.path.exists(filename):
        os.remove(filename)
else:
    if os.path.exists(filename):
        raise Exception(">> The file already exists, and set to not overwrite.")

h5file = h5.File(filename, 'w')  # New hdf5 file I want to create

# base attributes
if verbose > 0: print("Setting base attrs...")
for key in base_attrb.keys():
    h5file.attrs[key] = base_attrb[key]

# group: Chombo_global
if verbose > 0: print("Setting 'chombo global' attrs...")
chg = h5file.create_group('Chombo_global')
for key in chombogloba_attrb.keys():
    chg.attrs[key] = chombogloba_attrb[key]

# group: levels
if verbose > 0: print("Setting levels...")
metadic = dict()
for il in range(base_attrb['num_levels']):
    level_id = 'level_{}'.format(int(il))
    if verbose > 0: print("     creating ", level_id)
    lev = h5file.create_group(level_id)
    level_attrb = levels_attrb[level_id]
    for key in level_attrb.keys():
        if verbose>2: print( "      key-att is", key,  "      value-att is  ", level_attrb[key])
        lev.attrs[key] = level_attrb[key]
    sl = lev.create_group('data_attributes')
    sl.attrs['ghost'] = data_attributes['ghost']
    sl.attrs['outputGhost'] = data_attributes['outputGhost']
    sl.attrs['comps'] = base_attrb['num_components']
    sl.attrs['objectType'] = data_attributes['objectType']

    lev.create_dataset("Processors", data=np.array([0]))
    if verbose > 2: print("     boxes is", boxes)
    boxes_lev = np.array(boxes[level_id])
    lev.create_dataset("boxes", data=boxes_lev)


    Nlev = params['N'] * 2 ** (il)
    dd_lev = params['L'] / Nlev
    boxes_lev = np.array(boxes[level_id].tolist())
    offsets = [0]
    fdset = []  # list containing all box-data (flatten)
    #lev.create_dataset("data:datatype=0", data=np.array([]), chunks=True , maxshape=(None, 1)  )
    lev.create_dataset("data:datatype=0", shape=(0,), chunks=True , maxshape=(None, )  )
    for ib, lev_box in enumerate(boxes_lev):
        if verbose > 2:
            print("  {} box of level {}".format(ib, il))
        X = np.arange(lev_box[0], lev_box[3]+1)
        Y = np.arange(lev_box[1], lev_box[4]+1)
        Z = np.arange(lev_box[2], lev_box[5]+1)
        cord_grid_check = False
        for ic, comp in enumerate(components):
            comp_grid = np.zeros((len(X), len(Y), len(Z)))
            try:
                cid = np.where(components_vals[:, 0] == comp)[0][0]
                eval = components_vals[cid, 1]
            except Exception as e:
                print(" !! component {} not found, values set to zero".format(comp))
                print("   Execption: ", e)
                eval = 0
            if callable(eval):
                # Create coordinate grids
                if not cord_grid_check:
                    cord_grid_check = True
                    x_cord_grid = comp_grid.copy()
                    y_cord_grid = comp_grid.copy()
                    z_cord_grid = comp_grid.copy()
                    # loop over all coords
                    for ix, px in enumerate(X):
                        for iy, py in enumerate(Y):
                            for iz, pz in enumerate(Z):
                                dcnt = 0.5  # cell centering 
                                x_cord_grid[ix, iy, iz] = (px + dcnt) * dd_lev
                                y_cord_grid[ix, iy, iz] = (py + dcnt) * dd_lev
                                z_cord_grid[ix, iy, iz] = (pz + dcnt) * dd_lev
                comp_grid = eval(x_cord_grid, y_cord_grid, z_cord_grid)
                if comp_grid.size != x_cord_grid.size:
                    print("data size does not agree with mesh size")
                    print(" grid size {}, data size {}".format(x_cord_grid.size, comp_grid.size))
                    raise ValueError
            else:
                try:
                    eval = float(eval)
                except ValueError:
                    print("data eval is not a function or digit  --> ", eval)
                    raise
                comp_grid = comp_grid.copy() + eval

            fc = comp_grid.T.flatten()
            
            levdat = lev["data:datatype=0"]
            levdat.resize((levdat.shape[0] + fc.shape[0]), axis=0)
            levdat[-fc.shape[0]:] = fc
                  

            if il ==0:
                metadic[comp] = [np.mean(fc), np.min(fc), np.max(fc)]
            else:
                metadic[comp] = [metadic[comp][0], \
                                np.min([np.min(fc),metadic[comp][1]]), \
                                np.max([np.max(fc),metadic[comp][2]])]				

            # Cleanining 
            # del  comp_grid
        offsets.extend([len(fdset)])

    offsets = np.array(offsets)
    lev.create_dataset("data:offsets=0", data=offsets)

h5file.close()


################# Print a Summary ######################

if verbose > 0:
	print("")
	print("Summary checkfile:  [average, min, max]")
	for comp in metadic:
		print(comp, "\t\t-->\t\t", metadic[comp])


###################### END #########################
