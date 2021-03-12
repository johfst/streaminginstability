import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.ticker as tck
import pandas as pd
import argparse
import plotutils as pu

parser = argparse.ArgumentParser(
        description="Plot particle density from tab file."
        )

parser.add_argument("datapath", help="Path to directory containing id directories")
parser.add_argument("tabindex", help="Tabfile index", type=int)
parser.add_argument("savepath", help="Name of output image file")
parser.add_argument("x1_blocks", help="Number of blocks in the x1 direction", type=int)
parser.add_argument("x2_blocks", help="Number of blocks in the x2 direction", type=int)
parser.add_argument("title", help="Title of plot")
parser.add_argument("vmin", help="Colorbar minumum", type=float)
parser.add_argument("vmax", help="Colorbar maximum", type=float)

args = parser.parse_args()

# load tabfile dataframe
# sort by rows, then columns
tab_df = pu.get_tab_df_multicore(
        args.datapath, 
        args.x1_blocks, 
        args.x2_blocks, 
        timestep=args.tabindex,
        ).sort_values(by=["j-zone", "i-zone"])

x_vals = np.array(tab_df["x1"])
z_vals = np.array(tab_df["x2"])
dpar_array = np.array(tab_df["dpar"])

i_vals = tab_df["i-zone"]
j_vals = tab_df["j-zone"]
i_length = max(i_vals) - min(i_vals) + 1
j_length = max(j_vals) - min(j_vals) + 1

dpar_matrix = np.flip(np.reshape(dpar_array, (j_length, i_length)), axis=0)

### plotting ###
plt.rcParams.update({"text.usetex" : True})

base = 2
print("Creating plot...")
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_xlabel("$x$")
ax.set_ylabel("$z$")
ax.set_title(args.title)

print("Setting locator & formatter...")
locator = tck.LogLocator(base=base)
formatter = tck.LogFormatter(base=base)

print("Plotting...")
plot = ax.imshow(
            dpar_matrix,
            extent=(
                x_vals.min(),
                x_vals.max(),
                z_vals.min(),
                z_vals.max(),
                ),
            aspect="auto",
            interpolation="none",
            norm=col.LogNorm(
                vmin=args.vmin,
                vmax=args.vmax,
                ),
            )

cbar = fig.colorbar(plot, ax=ax, ticks=locator, format=formatter)
cbar.set_label(
        r"$\frac{\Sigma_p}{\left< \Sigma_p \right>}$",
        rotation="horizontal",
        fontsize=16,
        labelpad=15,
        )

vmin_order_of_magnitude = np.ceil(np.log10(args.vmin))
vmax_order_of_magnitude = np.floor(np.log10(args.vmax))

ticks = np.append(
        args.vmin,
        [np.logspace(
            vmin_order_of_magnitude,
            vmax_order_of_magnitude,
            6,
            base=10,
            ),
            args.vmax]
        )
plt.tight_layout()

print("Saving figure...")
fig.savefig(args.savepath)
print("Done.")
