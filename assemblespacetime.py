import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid, ImageGrid
import spacetime as st
import matplotlib.ticker as tck
import matplotlib.colors as col

VMIN = 0.01
VMAX = 50
prefixpath = "/home/theia/jatkins5/2021-02-04/"

plt.rcParams.update({"text.usetex" : True})
#fig, axs = plt.subplots(ncols=3, nrows=1)
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 8))
#fig = plt.figure(figsize=(6, 4))
"""
axgrid = AxesGrid(
        fig,
        111,
        nrows_ncols=(1, 3),
        cbar_mode="single",
        cbar_location="right",
        cbar_pad=0.1,
        )
grid = ImageGrid(
        fig,
        111,
        nrows_ncols=(1, 3),
        )
        """

titles = [r"$Z = 0.01$", r"$Z = 0.02$", r"$Z = 0.04$"]
plot = None
for i in range(1, 4):
    path = prefixpath + f"run_{i}/"
    surfdenscsv = path + "surfdens.csv"
    posvalsfile = path + "positions.csv"

    surf_dens_arr, pos_vals, orbit_vals = st.make_surf_dens_arr(surfdenscsv, posvalsfile)
    pos_min = min(pos_vals)
    pos_max = max(pos_vals)
    orbits_min = min(orbit_vals)
    #orbits_max = max(orbit_vals)
    orbits_max = 180

    plot = ax[i-1].imshow(
            surf_dens_arr,
            vmin=VMIN,
            vmax=VMAX,
            extent=(
                pos_min,
                pos_max,
                orbits_min,
                orbits_max
                ),
            aspect='auto',
            interpolation="none",
            norm=col.LogNorm(
                 vmin=VMIN,
                 vmax=VMAX,
                ),
            )
    ax[i-1].set_title(titles[i-1])
    ax[i-1].set_xlabel(r"$x/H$")
    """
    plot = st.make_spacetimeax(
        grid[i-1],
        surf_dens_arr,
        pos_vals.min(),
        pos_vals.max(),
        orbit_vals.min(),
        orbit_vals.max(),
        title,
        VMIN,
        VMAX,
        )
        """

#st.make_colorbar(fig, axs, plot)
#cbar = axgrid[i-1].cax.colorbar(plot)
#cbar = axgrid.cbar_axes[0].colorbar(plot)

#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.9, 0.15, 0.05, 0.7])
#fig.colorbar(plot, ax=ax.ravel().tolist())

"""
grid[-1].cax.colorbar(plot)
grid[-1].cax.toggle_label(True)
"""

base=2
locator = tck.LogLocator(base=base)
formatter = tck.LogFormatter(base=base)

cbar = fig.colorbar(plot, ax=ax[-1], ticks=locator, format=formatter)
cbar.set_label(
    r"$\frac{\Sigma_p}{\left< \Sigma_p \right>}$",
    rotation='horizontal',
    fontsize=16,
    labelpad=15,
    )
ax[0].set_ylabel("Number of Orbits")
plt.tight_layout()
fig.savefig("3stplot.png")
