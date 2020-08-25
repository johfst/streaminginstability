import argparse
import numpy as np
import matplotlib.pyplot as plt
import plotutils as pu
from sys import exit

parser = argparse.ArgumentParser(
        description="Tool to plot the particle density at a given timestep"
        )

parser.add_argument("simfolder", help="Path to folder containing id folders")
parser.add_argument("savepath", help="Name of output plot file")
parser.add_argument("timestep", type=int, help="Timestep")
parser.add_argument("--dt", type=int, help="Timesteps between tab outputs, default 1", default=1)
parser.add_argument("x1_blocks", type=int, help="Number of blocks in the x1 direction")
parser.add_argument("x2_blocks", type=int, help="Number of blocks in the x2 direction")

args = parser.parse_args()

# for now, error at any timestep that doesn't actually exist
# in the future, could instead find closest real timestep and output that
if args.timestep % args.dt != 0:
    print(f"Error: timestep {args.timestep} doesn't exist in tab format.")
    exit(1)

t_index = int(args.timestep / args.dt)

tab_df = pu.get_tab_df_multicore(args.simfolder, args.x1_blocks, args.x2_blocks, t_index)
print(f"Got tab dataframe for tabfile index {t_index}")
