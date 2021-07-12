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
                )
                and directory[:2] == "id"
                and directory[2:].isdigit()],
            key=lambda directory: int(directory[2:]),
            )

def get_last_tab_file(id_foldername, id_directory_num=0, dim3=False):
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

    dimstr = None
    if dim3:
        dimstr = "3d"
    else:
        dimstr = "2d"

    return f"Par_Strat{dimstr}{id_filename_modifier}.{last_tab_time:04}.tab"

def get_tab_file(idfolder, timestep_id):
    """
    Take a path to a(n) (id) directory, return the name of the tab file
    with the given timestep id
    """

    filelist = [f for f in os.listdir(idfolder) if os.path.isfile(os.path.join(idfolder, f))]
    tabfilelist = [f for f in filelist if "tab" in f]
    return [f for f in tabfilelist if f"{timestep_id:04}" in f][0]

def get_iddir_2d_chunk(id_dirs, id_offset, chunk2dsize):
    return [id_dir for id_dir in id_dirs 
            if id_offset <= get_idnum_from_foldername(id_dir) < id_offset+chunk2dsize]

def get_idnum_from_foldername(id_foldername):
    return int(id_foldername[2:])

def get_tab_df_multicore(foldername, x1_blocks, x2_blocks, timestep=None, x3_block_offset=None):
    """
    Take a path to a multicore run folder and return a dataframe containing
    all of the tabfile data
    x1_blocks: Number of blocks (for MPI) in the x1 direction
    x2_blocks: Number of blocks (for MPI) in the x2 direction
    timestep: Tab file index
    """

    # This is a sorted list of the id directories in the parent folder
    id_directories = get_sorted_id_dirs(foldername)
    if x3_block_offset is not None:
        chunk2dsize = x1_blocks * x2_blocks
        id_offset = x3_block_offset * chunk2dsize
        id_directories = get_iddir_2d_chunk(id_directories, id_offset, chunk2dsize)

    i_offset = 0
    j_offset = 0
    previous_x2_block_index = 0
    master_df = None
    for id_directory in id_directories:
        id_directory_num = get_idnum_from_foldername(id_directory)
        tab_filename = None
        if timestep is None:
            tab_filename = get_last_tab_file(
                    os.path.join(foldername, id_directory),
                    id_directory_num,
                    dim3=(x3_block_offset is not None),
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
        x2_block_index = np.floor(float(id_directory_num % chunk2dsize) / x1_blocks)

        # compute j_offset for this block
        if x2_block_index > previous_x2_block_index:
            j_offset += max(df["j-zone"]) - min(df["j-zone"]) + 1
        # add accumulated offset
        df["i-zone"] += i_offset
        df["j-zone"] += j_offset

        # compute k_offset for this block (if needed)
        if x3_block_offset is not None:
            df["k-zone"] += ( max(df["k-zone"]) - min(df["k-zone"]) + 1 ) * x3_block_offset

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

def get_tab_df_multicore_3d(foldername, x1_blocks, x2_blocks, x3_blocks, timestep=None):
    master_df = None
    for x3_index in range(x3_blocks):
        # get 2d slice (think of a plane in the x1 and x2 directions)
        slice_2d_df = get_tab_df_multicore(
                foldername, x1_blocks, x2_blocks, timestep=timestep, x3_block_offset=x3_index
                )
        
        if master_df is None:
            master_df = slice_2d_df
        else:
            master_df = master_df.append(slice_2d_df)

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
        

        dzarr = np.array([
            i_rows.iloc[j]["x2"] - i_rows.iloc[j-1]["x2"] for j in range(1,len(i_rows))
            ])
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

def integrate_df_3d(df, direction):
    """
    Integrate a 3d dataframe over an arbitrary direction
    honestly, this would probably be better if it turned the
    input df into a numpy cube, because i bet you could vectorize
    it and it would be faster, but whatever

    df : dataframe containing positions, dpar, etc.
    direction : one of "x1", "x2", "x3"
    """
    def direction_to_index(direction):
        return int(direction[1:]) - 1
    def direction_to_index_label(direction):
        return chr( direction_to_index(direction) + 105 )
    def get_direction_index_compl(dir_indx):
        return sorted([ (dir_index + n) % 3 for n in (1,2) ])
    def index_to_direction(dir_indx):
        return "x" + str(dir_indx + 1)

    dir_index = direction_to_index(direction)

    # directions orthogonal to integrating direction
    dir1_index, dir2_index = get_direction_index_compl(dir_index)
    dir1 = index_to_direction(dir1_index)
    dir2 = index_to_direction(dir2_index)
    dir1_index_label = direction_to_index_label(dir1) + "-zone"
    dir2_index_label = direction_to_index_label(dir2) + "-zone"

    dir1_min = min(df[dir1_index_label])
    dir1_max = max(df[dir1_index_label])
    dir2_min = min(df[dir2_index_label])
    dir2_max = max(df[dir2_index_label])

    new_rows = []
    for n in range(dir1_min, dir1_max+1):
        for m in range(dir2_min, dir2_max+1):
            locator = (df[dir1_index_label] == n) & (df[dir2_index_label] == m)
            nm_location_df = df.loc[locator]
            row = [n, m, nm_location_df.iloc[0][dir2], nm_location_df.iloc[0][dir1]]

            d_dir = np.abs(nm_location_df.iloc[1][direction] - nm_location_df.iloc[0][direction])
            compute_avgs = lambda arr : np.array([
                    0.5*(arr[l] + arr[l+1]) for l in range(0, len(arr)-1)
                    ])
            def compute_integral(col):
                vals = np.array(nm_location_df[col])
                avgs = compute_avgs(vals)
                return sum(avgs * d_dir)

            dpar_integral = compute_integral("dpar")
            d_integral = compute_integral("d")
            row.append(dpar_integral)
            row.append(d_integral)

            new_rows.append(row)

    columns = [dir1_index_label, dir2_index_label, dir1, dir2, "dpar", "d"]
    integrated_df = pd.DataFrame(new_rows, columns=columns)
    return integrated_df

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
    header_line = get_tab_header_line(filename)
    sep = "(?<!#)\s+"
    df = pd.read_csv(filename, sep=sep, header=header_line, engine="python")
    
    # column names in the file suck, so renaming them
    newcols = {}
    for col in list(df):
        newname = col[col.find("=")+1:]
        newcols[col] = newname
    df = df.rename(columns=newcols)
    df = df.apply(pd.to_numeric, errors="coerce")
    return df

def get_tab_header_line(filename):
    header_line = 0
    with open(filename, "r") as f:
        line = f.readline()
        while line != "":
            if line.find("[1]") != -1:
                break
            else:
                header_line += 1
                line = f.readline()
    return header_line

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
