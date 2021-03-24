import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.ticker as tck
from matplotlib import rc
import plotutils as pu
import argparse

### space-time plots, as in Carrera fig. 4 ###

def make_spacetimefig(surfdenscsv, posvalsfile, title, vmin, vmax):
    print("Making spacetime plot...")

    surf_dens_arr = np.loadtxt(surfdenscsv, delimiter=",")
    pos_vals = np.loadtxt(posvalsfile, delimiter=",")
    t_vals = surf_dens_arr[:,0] # the first column
    orbits_vals = t_vals / (2 * np.pi)
    surf_dens_arr = surf_dens_arr[:,1:] # all but the first column

    print("Calculating mean & scaling...")
    mean_arr = np.array([np.average(row) for row in surf_dens_arr])
    surf_dens_arr_scaled = np.flip(np.array([surf_dens_arr[i] / mean_arr[i] for i in range(len(mean_arr))]), 0)

    ### space-time plots, continued ###
    plt.rcParams.update({"text.usetex" : True})

    base = 2
    print("Creating plot...")
    fig, ax = plt.subplots(figsize=(4,8))
    ax.set_xlabel("$x$")
    ax.set_ylabel("revolutions")
    ax.set_title(title)

    print("Setting locator & formatter...")
    locator = tck.LogLocator(base=base)
    formatter = tck.LogFormatter(base=base)

    print("Plotting...")
    plot = ax.imshow(
            surf_dens_arr_scaled,
            extent=(
                pos_vals.min(),
                pos_vals.max(),
                orbits_vals.min(),
                orbits_vals.max()
                ),
            aspect='auto',
            interpolation="none",
            norm=col.LogNorm(
            #    vmin=surf_dens_arr_scaled.min(),
            #    vmax=surf_dens_arr_scaled.max()
                 vmin=vmin,
                 vmax=vmax,
                ),
            #vmin=2**-0.5,
            #vmax=2**5
            )

    cbar = fig.colorbar(plot, ax=ax, ticks=locator, format=formatter)
    cbar.set_label(
        r"$\frac{\Sigma_p}{\left< \Sigma_p \right>}$",
        rotation='horizontal',
        fontsize=16,
        labelpad=15,
        )

    vmin_order_of_magnitude = np.ceil(np.log10(vmin))
    vmax_order_of_magnitude = np.floor(np.log10(vmax))
    ticks = np.append(
            vmin,
            [np.logspace(
                vmin_order_of_magnitude,
                vmax_order_of_magnitude,
                6,
                base=10,
                ),
                vmax
                ])
    #cbar.set_ticks(ticks)
    #cbar.set_ticklabels(list(map("{:.2f}".format, ticks)))
    plt.tight_layout()

    print("Spacetime figure done.")
    return fig

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description="Save a plot showing the radius on the x axis, time on the y axis, and \
                    surface density on the color axis."
    )

    parser.add_argument("surfdenscsv", help="Path to surface density csv file")
    parser.add_argument("savepath", help="Name of output image file")
    parser.add_argument("posvalsfile", help="Path to file containing x1 position values")
    parser.add_argument("title", help="Title of plot")
    parser.add_argument("vmin", help="Colorbar minimum", type=float)
    parser.add_argument("vmax", help="Colorbar maximum", type=float)

    args = parser.parse_args()

    surfdenscsv = args.surfdenscsv
    savepath = args.savepath
    posvalsfile = args.posvalsfile
    title = args.title

    fig = make_spacetimefig(surfdenscsv, posvalsfile, title, args.vmin, args.vmax)

    print("Saving figure...")
    fig.savefig(savepath)
    print("Done.")
