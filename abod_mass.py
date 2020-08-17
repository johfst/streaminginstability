import plotutils as pu
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from astropy import units as u
import argparse
from sys import exit

LENGTH_SCALE = 3.84e6 * u.meter
GAS_DENSITY = 12.0 * u.kilogram / u.meter**3
ANGULAR_FREQUENCY = 2.47e-4 / u.second
GRAVITATIONAL_CONSTANT = 6.67408e-11 * u.meter**3 / u.kilogram / u.second**2

parser = argparse.ArgumentParser(
        description="Calculating mass of largest lunar seed, as in Abod"
        )
parser.add_argument("surfdenscsv", help="Location of surface density csv file")
parser.add_argument("--plotdir", 
        help="Directory to save plots, default is to not do so. Requires --posfile")
parser.add_argument("--posfile", help="Location of x positions csv file")
parser.add_argument("--timestep", 
        help="Timestep to compute mass for, default is the last one",
        type=int,
        default=-1,
        )

args = parser.parse_args()

if args.plotdir is not None and args.posfile is None:
    print("Error: plot output requires a position values file")
    exit(1)

t_vals, surf_dens_df = pu.load_surfdenscsv(args.surfdenscsv)

timestep_index = None
if args.timestep > -1:
    if args.timestep not in t_vals:
        print("Error: provided timestep not present in data.")
        exit(1)

    # get the index of surface_density corresponding to the given timestep
    timestep_index = np.where(t_vals == args.timestep)[0][0]
else:
    timestep_index = -1

surface_density = surf_dens_df[timestep_index] * LENGTH_SCALE * GAS_DENSITY

"""
surf_dens_df = np.loadtxt(args.surfdenscsv, delimiter=",")
x1 = np.loadtxt(args.posfile, delimiter=",") * LENGTH_SCALE
surface_density = surf_dens_df[:, 1:] * LENGTH_SCALE * GAS_DENSITY
"""

#surf_dens_df = pu.get_tab_df_multicore(foldername, x1_blocks=20, x2_blocks=2)
#print("Got surface density data frame")
#x1 = surf_dens_df["x1"] * LENGTH_SCALE
#surface_density = surf_dens_df["dpar"] * LENGTH_SCALE * GAS_DENSITY


"""
This is equation 23 from Abod et. al. 2019
"""
mass_distribution = 4 * np.pi**5 * (
        GRAVITATIONAL_CONSTANT**2 * surface_density**3 /
        ANGULAR_FREQUENCY**4
        )

if args.plotdir is not None:
    x1 = np.loadtxt(args.posfile, delimiter=",") * LENGTH_SCALE
    plt.rcParams.update({"text.usetex" : True})

    surface_density_fig, surface_density_ax = plt.subplots()
    surface_density_ax.set_ylabel(r"$\Sigma_p$ (kg m$^{-2}$)")
    surface_density_ax.set_xlabel("$x$ (m)")
    surface_density_ax.set_title("Disk Surface Density vs. Radius")
    surface_density_ax.tick_params(direction="in")
    surface_density_ax.plot(x1, surface_density)
    surface_density_fig.tight_layout()
    pu.prompt_savefig(surface_density_fig, f"{args.plotdir}/surfacedensity.png")


    mass_distribution_fig, mass_distribution_ax = plt.subplots()
    mass_distribution_ax.set_xlabel("$x$ (m)")
    mass_distribution_ax.set_ylabel("$M_G$ (kg)")
    mass_distribution_ax.set_title("Mass Distribution Over Disk Radius")
    mass_distribution_ax.tick_params(direction="in")
    mass_distribution_ax.plot(x1, mass_distribution)
    mass_distribution_fig.tight_layout()
    pu.prompt_savefig(mass_distribution_fig, f"{args.plotdir}/massdistribution.png")

max_mass = max(mass_distribution)
print(f"max mass: {max_mass}")
