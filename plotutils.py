import numpy as np
import subprocess
import pandas as pd
import os
import sys

def prompt_savefig(figure, filename):
    """
    Take a filename for a figure, and prompt the user to overwrite an already existing
    file with that name if needed
    """
    if os.path.exists(filename):
        prompt = input(f"{filename} already exists, overwrite? [y/N]> ")
        if not (prompt.lower() == "y" or prompt.lower() == "yes"):
            prompt_savefig(figure, input("Enter new filename: "))
    
    figure.savefig(filename)


def get_sorted_id_dirs(foldername):
    """
    Take a path to a multicore run folder and return a sorted list of id directories
    """

    return sorted(
            [directory for directory in os.listdir(foldername) if os.path.isdir(
                os.path.join(foldername, directory)
                )],
            key=lambda directory: int(directory[2:]),
            )

def get_last_tab_file(id_foldername, id_directory_num=0):
    """
    Take a path to a(n) (id) directory, return the name of the tab file at the last
    time index
    """
    
    # final index on the tab files
    last_tab_time = get_last_tab_time(id_foldername)

    # build the tabfile name; the ones in the 0th id directory
    # aren't labeled by id, but the rest are
    id_filename_modifier = ""
    if id_directory_num > 0:
        id_filename_modifier = f"-id{id_directory_num}"
    return f"Par_Strat2d{id_filename_modifier}.{last_tab_time:04}.tab"

def get_tab_file(idfolder, timestep_id):
    """
    Take a path to a(n) (id) directory, return the name of the tab file
    with the given timestep id
    """

    filelist = [f for f in os.listdir(idfolder) if os.path.isfile(os.path.join(idfolder, f))]
    tabfilelist = [f for f in filelist if "tab" in f]
    return [f for f in tabfilelist if f"{timestep_id:04}" in f][0]



def get_tab_df_multicore(foldername, x1_blocks, x2_blocks, timestep=None):
    """
    Take a path to a multicore run folder and return a dataframe containing
    all of the tabfile data
    x1_blocks: Number of blocks (for MPI) in the x1 direction
    x2_blocks: Number of blocks (for MPI) in the x2 direction
    timestep: Tab file index
    """

    # This is a sorted list of the id directories in the parent folder
    id_directories = get_sorted_id_dirs(foldername)

    i_offset = 0
    j_offset = 0
    previous_x2_block_index = 0
    master_df = None
    for id_directory_num, id_directory in enumerate(id_directories):
        tab_filename = None
        if timestep is None:
            tab_filename = get_last_tab_file(
                    os.path.join(foldername, id_directory),
                    id_directory_num,
                    )
        else:
            tab_filename = get_tab_file(os.path.join(foldername, id_directory), timestep)

        df = get_df(os.path.join(foldername, id_directory, tab_filename))

        """
        Need to add an offset to the cell indices to put them in the right place, since
        cell indices are computed by athena relative to the CPU space (blocks), but we need
        them relative to the entire domain.
        """
        x1_block_index = id_directory_num % x1_blocks
        x2_block_index = np.floor(float(id_directory_num) / x1_blocks)

        # compute j_offset for this block
        if x2_block_index > previous_x2_block_index:
            j_offset += max(df["j-zone"]) - min(df["j-zone"]) + 1
        # add accumulated offset
        df["i-zone"] += i_offset
        df["j-zone"] += j_offset

        """
        Compute the row (x2_block_index) and column (x1_block_index) of the current file's
        block. Athena creates id folders first along x1, then along x2, which is why these
        formulae work.
        """

        # compute i offset for next block
        if x1_block_index < x1_blocks - 1:
            i_offset += max(df["i-zone"]) - min(df["i-zone"]) + 1
        else:
            i_offset = 0

        previous_x2_block_index = x2_block_index

        if master_df is None:
            master_df = df
        else:
            master_df = master_df.append(df)

    return master_df

def compute_surface_dens(df):
    """
    Takes a dataframe containing the density at each point, returns
    a dataframe containing the surface density for each x
    """
    trim = df[["i-zone", "j-zone", "x1", "x2", "d", "dpar"]]
    #print(len(trim))
    imin = min(trim["i-zone"])
    imax = max(trim["i-zone"])

    dens = []
    for i in range(imin, imax+1):
        i_rows = trim.loc[(trim["i-zone"] == i)]
        row = [i, i_rows.iloc[0]["x1"]]
        

        dzarr = np.array([i_rows.iloc[j]["x2"] - i_rows.iloc[j-1]["x2"] for j in range(1,len(i_rows))])
        par_surf_dens = sum(np.array(i_rows["dpar"])[1:] * dzarr)
        #print(surf_dens)
        row.append(par_surf_dens)
        
        #print(dzarr)
        #print(i_rows["d"])
        #print(np.array(i_rows["d"])[1:])
        gas_surf_dens = sum(np.array(i_rows["d"])[1:] * dzarr)
        row.append(gas_surf_dens)
        
        dens.append(row)

    densdf = pd.DataFrame(dens, columns=["i-zone", "x1", "dpar", "d"])
    
    return densdf

def get_last_tab_time(folder):
    lsproc = subprocess.Popen(["ls", folder], stdout=subprocess.PIPE)
    grepproc = subprocess.Popen(["grep", "tab"], stdin=lsproc.stdout, stdout=subprocess.PIPE)
    tailproc = subprocess.Popen(["tail", "-1"], stdin=grepproc.stdout, stdout=subprocess.PIPE)

    stdout, stderr = tailproc.communicate()
    #print(stdout.decode("utf-8"))
    return int( stdout.decode("utf-8").split(".")[1] )

def load_surfdenscsv(filename):
    filevals = np.loadtxt(filename, delimiter=",")
    t_vals = filevals[:,0] # the first column
    surf_dens_arr = filevals[:,1:] # all but the first column
    return t_vals, surf_dens_arr
    
def get_df(filename):
    # i == x1 == r, j == x2 == z
    head = ["i-zone", "j-zone", "x1", "x2", "d", "M1", "M2", "M3", "dpar", "M1par", "M2par", "M3par"]
    df = pd.read_csv(filename, sep="\s+", names=head, header=5, usecols=range(0,12))
    df = df.apply(pd.to_numeric, errors="coerce")
    return df

def notify(msg):
    notifystr = "notify-send -i /home/jeremy/Pictures/icons/python-icon-32.png Python '{0}'".format(msg)
    prc = subprocess.Popen(shlex.split(notifystr), stdout=subprocess.PIPE)

def Q(z, jmax):
    ans = 0
    for j in range(1, jmax+1):
        ans += (-1)**(j-1) * np.exp(-2 * z**2 * j**2)
    return 2*ans

def get_max_dist(cumdist):
    uniformline = np.array(range(0,len(cumdist)))/(len(cumdist)+0.0)
    return max(cumdist - uniformline)

def get_p(D, n):
    return Q(D*np.sqrt(n), 10000)
