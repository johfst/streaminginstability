import numpy as np
import pandas as pd
import sys
import os
import plotutils as pu
import argparse
from time import time

### space-time plots, as in Carrera fig. 4 ###
start_time = time()

parser = argparse.ArgumentParser(
        description=("Save a csv containing as rows the surface density "
        "at each specified time step.")
        )
parser.add_argument("simfolder", help="Path to folder containing id folders")
parser.add_argument("savepath", help="Name of output csv file")
parser.add_argument("x1_blocks", type=int, help="Number of blocks in the x1 direction")
parser.add_argument("x2_blocks", type=int, help="Number of blocks in the x2 direction")
parser.add_argument("--tmin", type=int, help="Starting tabfile index, default 0", default=0)
parser.add_argument("--tmax", type=int, 
        help="Final tabfile index, default is highest one present", default=None)
parser.add_argument("--dt", type=int, 
        help="Number of time steps between each tabfile output, default 1", default=1)
parser.add_argument("--posfile", 
        help="Output a csv containing the position of density values along the x1 axis",
        default=None)

args = parser.parse_args()


simfolder = args.simfolder
savepath = args.savepath
x1_blocks = args.x1_blocks
x2_blocks = args.x2_blocks

tmin = args.tmin
tmax = args.tmax
if tmax is None:
    tmax = pu.get_last_tab_time(os.path.join(simfolder, "id0"))
if tmin > tmax:
    print("Error: tmin must be less than or equal to tmax!")
    print("Exiting...")
    sys.exit(1)
index_dt = args.dt
posfile = args.posfile

pos_vals = None
# these are the true timesteps, for each step that has tab data
t_vals = np.array(list(range(tmin, index_dt*(tmax - tmin)+tmin + 1, index_dt)))
print(t_vals)

for i, t_index in enumerate(range(tmin, tmax+1)):
    print(f"Importing t_index = {t_index}")
    tab_df = pu.get_tab_df_multicore(simfolder, x1_blocks, x2_blocks, t_index)

    if t_index == tmin:
        pos_vals = np.array(tab_df["x1"])

        if posfile is not None:
            np.savetxt(posfile, pos_vals, delimiter=",")
            print(f"Saved position values in {posfile}.")

    print(f"Calculating surface density for time index {t_index}")
    surface_density_t = pu.compute_surface_dens(tab_df)
    print(surface_density_t)
    print(np.array(surface_density_t["dpar"]))
    print(np.array(surface_density_t["dpar"]).shape)
    # prepend row with timestep
    row = np.insert(np.array(surface_density_t["dpar"]), 0, t_vals[i])

    with open(savepath, "ab") as surfdensfile:
        np.savetxt(surfdensfile, row, delimiter=",")
    print(f"Saved surface density for time index {t_index} successfully.")

end_time = time()
print(f"Done ({end_time - start_time:.1f} seconds)")
