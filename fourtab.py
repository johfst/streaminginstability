import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colorbar import Colorbar
import matplotlib.colors as col
import matplotlib.ticker as tck
import argparse
import plotutils as pu

parser = argparse.ArgumentParser(
        description="Plot particle density from tab file."
        )

parser.add_argument("datapath", help="Path to directory containing id directories")
parser.add_argument("tabindices", help="Four tabfile indices",
        type=int, nargs=4)
parser.add_argument("savepath", help="Name of output image file")
parser.add_argument("x1_blocks", help="Number of blocks in the x1 direction", type=int)
parser.add_argument("x2_blocks", help="Number of blocks in the x2 direction", type=int)
parser.add_argument("--x3_blocks", help="Number of blocks in the x3 direction (for 3d data)",
        type=int)
parser.add_argument("titles", help="Titles of plots", nargs=4)
parser.add_argument("--vmin", help="Colorbar minumum", type=float, default=0.01)
parser.add_argument("--vmax", help="Colorbar maximum", type=float, default=50)
"""
parser.add_argument("--peakfile",
        help="Output file from PLAN; will draw circles around peak locations")
        """
args = parser.parse_args()

tab_dfs = []
if args.x3_blocks is None:
    # 2d only

    for tabindex in args.tabindices:
        # load tabfile dataframe
        # sort by rows, then columns
        tab_df = pu.get_tab_df_multicore(
                args.datapath, 
                args.x1_blocks, 
                args.x2_blocks, 
                timestep=tabindex,
                ).sort_values(by=["j-zone", "i-zone"])
        tab_dfs.append(tab_df)


else:
    # 3d

    for tabindex in args.tabindices:
        tab_df_3d = pu.get_tab_df_multicore_3d(
                args.datapath,
                args.x1_blocks,
                args.x2_blocks,
                args.x3_blocks,
                timestep=tabindex
                ).sort_values(by=["k-zone", "j-zone", "i-zone"])

        tab_df = pu.integrate_df_3d(tab_df_3d, "x3").sort_values(by=["j-zone", "i-zone"])
        tab_dfs.append(tab_df)

x_vals = np.array(tab_dfs[0]["x1"])
y_vals = np.array(tab_dfs[0]["x2"])
dpar_arrays = [np.array(tab_df["dpar"]) for tab_df in tab_dfs]

indx_length = lambda arr : max(arr) - min(arr) + 1
i_length = indx_length(tab_dfs[0]["i-zone"])
j_length = indx_length(tab_dfs[0]["j-zone"]) 

"""
i_vals = tab_df["i-zone"]
j_vals = tab_df["j-zone"]
i_length = max(i_vals) - min(i_vals) + 1
j_length = max(j_vals) - min(j_vals) + 1
"""
dpar_matrices = [np.flip(np.reshape(dpar_array, (j_length, i_length)), axis=0) 
        for dpar_array in dpar_arrays]



### plotting ###
plt.rcParams.update({"text.usetex" : True})

fig = plt.figure()
grid = ImageGrid(
        fig,
        111,
        nrows_ncols=(2,2),
        axes_pad=0.15,
        share_all=True,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="7%",
        cbar_pad=0.15,
        )

base = 2
locator = tck.LogLocator(base=base)
formatter = tck.LogFormatter(base=base)

for ax_index, ax in enumerate(grid):

    dpar_matrix = dpar_matrices[ax_index]
    plot = ax.imshow(
                dpar_matrix,
                extent=(
                    x_vals.min(),
                    x_vals.max(),
                    y_vals.min(),
                    y_vals.max(),
                    ),
                aspect="auto",
                interpolation="none",
                norm=col.LogNorm(
                    vmin=args.vmin,
                    vmax=args.vmax,
                    ),
                )
    ax.set_xlabel("$x \;(H)$")
    ax.set_ylabel("$y \;(H)$")
    ax.set_title(args.titles[ax_index])

    ax.tick_params(axis="both", which="major", length=4, labelsize=9)

ax.cax.cla()
cbar = Colorbar(ax.cax, plot, ticks=locator, format=formatter)
ax.cax.toggle_label(True)

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

#plt.tight_layout(h_pad=10)

fig.savefig(args.savepath, dpi=200, bbox_inches="tight")
