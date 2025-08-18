#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib
import os

def read_command_line():
    parser = argparse.ArgumentParser(
        description="usage: ./plasmonX_analysis [options]",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-i", "--input_file", type=str, required=True,
                        help="input file (tar.gz from plasmonX execution)")
    parser.add_argument("-w", "--what", type=str, required=True,
                        help="what to extract. Possible values:\n"
                             "  1. xyz  (xyz file)\n"
                             "  2. pqr  (PQR file for charges)\n"
                             "  3. field (PLT/GNU file for field)\n"
                             "  4. density (PLT/GNU file for charge density)")
    parser.add_argument("-n", "--num_ex_freq", type=int, default=0,
                        help="[integer]\n"
                             "number of frequencies for PQR or FIELD or DENSITY")
    parser.add_argument("-freq", "--frequencies", type=str,
                        help="[float]\n"
                             "frequencies (eV) for PQR or FIELD or DENSITY")
    parser.add_argument("-min_freq", "--min_freq", type=float, default=-1.0,
                        help="[float]\n"
                             "min. freq (eV) \n"
                             "for lore or gaus (start freq for conv. [default = 0.0d0])\n"
                             "for PQR or FIELD or DENSITY (alternative to -freq)")
    parser.add_argument("-max_freq", "--max_freq", type=float, default=-1.0,
                        help="[float]\n"
                             "max. freq (eV) \n"
                             "for lore or gaus (end freq for conv. [default = 5.0d0])\n"
                             "for PQR or FIELD or DENSITY (alternative to -freq)")
    parser.add_argument("-step", "--step", type=float, default=-1.0,
                        help="[float]\n"
                             "step (eV) \n"
                             "for lore or gaus (graining for conv. [default = 1.0d-3])\n"
                             "for PQR or FIELD or DENSITY (alternative to -freq)")
    parser.add_argument("-scale_e0", "--scale_e0", type=float, default=1.0,
                        help="[float]\n"
                             "scale electric field intensity in input\n"
                             "this affects some properties, such as force [default = 1.0d0]")
    parser.add_argument("-plane", "--plane", type=str, default='null',
                        help="which plane to plot. Possible values:\n"
                             "  1. xy (xy plane)\n"
                             "  2. xz (xz plane)\n"
                             "  3. yz (yz plane)")
    parser.add_argument("-n_plane", "--n_plane", type=int, default=10,
                        help="number of planes. [default = 10]")
    parser.add_argument("-step_plane", "--step_plane", type=float, default=1.0,
                        help="distance btw planes. (Ang) [default = 1.0Ang]")
    parser.add_argument("-start", "--start", type=float, default=0.0,
                        help="starting coord for planes. (Ang) [default = 0.0Ang]")
    parser.add_argument("-field_dir", "--field_dir", type=str, default='all',
                        help="direction of the external field. Possible values:\n"
                             "  1. x (field // x)\n"
                             "  2. y (field // y)\n"
                             "  3. z (field // z)")
    parser.add_argument("-nx_points", "--nx_points", type=int, default=200,
                        help="Grid Points on X [default = 200]")
    parser.add_argument("-ny_points", "--ny_points", type=int, default=200,
                        help="Grid Points on Y [default = 200]")
    parser.add_argument("-nz_points", "--nz_points", type=int, default=200,
                        help="Grid Points on Z [default = 200]")
    parser.add_argument("-volume", "--volume", action="store_true",
                        help="Request calculation of effective volume")
    parser.add_argument("-min_grid", "--min_grid", type=str, default='0.0,0.0,0.0',
                        help="x,y,z min. values of the grid.\n" 
                             "For plane calculations, specify the values in the 2 required axis (-plane).\n"
                             'Syntax: -min_grid="0.0,0.0,0.0".\n'
                             "[mandatory if -volume]")
    parser.add_argument("-max_grid", "--max_grid", type=str, default='0.0,0.0,0.0',
                        help="x,y,z max. values of the grid. \n"
                             "For plane calculations, specify the values in the 2 required axis (-plane).\n"
                             'Syntax: -max_grid="0.0,0.0,0.0".\n'
                             "[mandatory if -volume]")
    parser.add_argument("-offset_grid", "--offset_grid", type=str, default='10.0',
                        help="Offset for constructing the grid [default = +- 10 Ang. from max/min atom coords]")
    parser.add_argument("-format_grid", "--format_grid", type=str, default='cube',
                        help="output format for calculations involving grid [default = cube]")
    parser.add_argument("-separate_q_mu", "--separate_q_mu", action="store_true",
                        help="Separate q and mu contributions in wfqfmu density [only for what=density]")
    parser.add_argument("-omp", "--omp_threads", type=int, default=0,
                        help="Number of OMP threads. Default: maximum available")
    parser.add_argument("-mem", "--memory_available", type=float, required=False,
                        help="Optional: available memory in GB [Default: the maximum available]")

    args = parser.parse_args()

    # Convert comma-separated strings into lists of floats
    min_grid = list(map(float, args.min_grid.split(',')))
    max_grid = list(map(float, args.max_grid.split(',')))

    # --- Input file check ---
    if not os.path.isfile(args.input_file):
        print(f"\nError: Input file '{args.input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    if not os.access(args.input_file, os.R_OK):
        print(f"\nError: Input file '{args.input_file}' is not readable.", file=sys.stderr)
        sys.exit(1)
    
    if not args.input_file.endswith('.tar.gz'):
        print(f"\nError: Input file '{args.input_file}' is not a .tar.gz archive.", file=sys.stderr)
        sys.exit(1)

    # Handling defaults and conditions
    if args.what not in ['xyz', 'pqr', 'field', 'density']:
        print(f"\nWhat option '{args.what}' not recognized"
              f"Valid options are: xyz, pqr, field, density")
        exit(1)

    # field_dir
    valid_directions = {'x', 'y', 'z', 'all'}
    if args.field_dir not in valid_directions:
        print(f"\nInvalid value for --field_dir: '{args.field_dir}'. "
              f"Valid options are: x, y, z, all")
        exit(1)

    if args.what in ['field', 'density']:
        if not args.num_ex_freq:
            print("\nYou did not specify the number of frequency [-n option]")
            exit(1)

        if args.volume and not args.min_grid:
            print("\n-volume but no -min_grid in input")
            exit(1)

        if args.volume and not args.max_grid:
            print("\n-volume but no -max_grid in input")
            exit(1)

    if args.separate_q_mu: 
        if args.what != "density": 
            print("\nSeparate q and mu only available for density")
            exit(1)
        if args.plane != "null": 
            print("\nSeparate q and mu currently not available for plane calculations")
            exit(1)

    frequencies = []
    if args.frequencies:
        frequencies = list(map(float, args.frequencies.split(',')))


    # Generate frequencies if step, min_freq, and max_freq are provided
    if args.step > 0 and args.min_freq < args.max_freq and args.num_ex_freq > 0:
        frequencies = [args.min_freq + i * args.step for i in range(int((args.max_freq - args.min_freq) / args.step) + 1)]

    return args, min_grid, max_grid, frequencies

def write_parameters_to_file(args, min_grid, max_grid, frequencies, n_omp_threads):
    with open('parameters.txt', 'w') as file:
        file.write(f"input_file = {args.input_file}\n")
        file.write(f"what = {args.what}\n")
        file.write(f"num_ex_freq = {args.num_ex_freq}\n")
        if frequencies and args.num_ex_freq <= 30:
            file.write(f"frequencies = {frequencies}\n")
        file.write(f"step = {args.step}\n")
        file.write(f"scale_e0 = {args.scale_e0}\n")
        file.write(f"plane = {args.plane}\n")
        file.write(f"n_plane = {args.n_plane}\n")
        file.write(f"step_plane = {args.step_plane}\n")
        file.write(f"start = {args.start}\n")
        file.write(f"field_dir = {args.field_dir}\n")
        file.write(f"nx_points = {args.nx_points}\n")
        file.write(f"ny_points = {args.ny_points}\n")
        file.write(f"nz_points = {args.nz_points}\n")
        file.write(f"volume = {args.volume}\n")
        if(min_grid != max_grid):
           file.write(f"min_grid = {min_grid}\n")
           file.write(f"max_grid = {max_grid}\n")
        else:
           file.write(f"min_grid = default\n")
           file.write(f"max_grid = default\n")
        file.write(f"offset_grid = {args.offset_grid}\n")
        file.write(f"format_grid = {args.format_grid}\n")
        file.write(f"separate_q_mu = {args.separate_q_mu}\n")
        file.write(f"n_omp_threads = {n_omp_threads}\n")


def print_banner(output_file, configurations):
    """
    Print the banner of plasmonX analysis

    Args:
        output_file (str): Path del file di output.
        configurations (dict): Dizionario delle configurazioni
    """
    sticks = "-" * 80  # Definisce 80 trattini
    # Calcola la lunghezza massima delle chiavi per l'allineamento
    max_key_length = max(len(key) for key in configurations.keys())

    with open(output_file, "a") as f:  # Apre il file in modalità scrittura
        f.write(f" {sticks}\n")
        f.write("                         _                                __  __ \n")
        f.write("                   _ __ | | __ _ ___ _ __ ___   ___  _ __ \ \/ / \n")
        f.write("                  | '_ \| |/ _` / __| '_ ` _ \ / _ \| '_ \ \  /  \n")
        f.write("                  | |_) | | (_| \__ \ | | | | | (_) | | | |/  \  \n")
        f.write("                  | .__/|_|\__,_|___/_| |_| |_|\___/|_| |_/_/\_\ \n")
        f.write("                  |_|                     _           _          \n")
        f.write("                         __ _ _ __   __ _| |_   _ ___(_)___      \n")
        f.write("                        / _` | '_ \ / _` | | | | / __| / __|     \n")
        f.write("                       | (_| | | | | (_| | | |_| \__ \ \__ \     \n")
        f.write("                        \__,_|_| |_|\__,_|_|\__, |___/_|___/     \n")
        f.write("                                            |___/                ")
        f.write("\n\n")
        f.write(f" {sticks}\n")
        f.write(" Tommaso Giovannini    University of Rome Tor Vergata, Rome, Italy\n")
        f.write(f" {sticks}\n")
        for key, value in configurations.items():
            f.write(f" {key.ljust(max_key_length)} : {value}\n")
        f.write(f" {sticks}\n")

def write_plot_script(file_python, csv_file, plane, what):
    svg_output = csv_file.replace(".csv", ".svg")
    png_output = csv_file.replace(".csv", ".png")
    dim_1 = plane[0]
    dim_2 = plane[1]

    if (what == "field") :
        sel_palette = "jet"
    elif (what == "density") :
        sel_palette = "seismic"

    with open(file_python, "w") as f:
        f.write("import numpy as np\n")
        f.write("import matplotlib.pyplot as plt\n")
        f.write("from scipy.interpolate import griddata\n\n")

        f.write(f"csv_file = '{csv_file}'\n")
        f.write(f"output = '{svg_output}'\n")
        f.write(f"output_png = '{png_output}'\n\n")

        f.write("data = []\n")
        f.write("with open(csv_file, 'r') as file:\n")
        f.write("    for line in file:\n")
        f.write("        if line.strip():\n")
        f.write("            data.append([float(line[:25]), float(line[28:53]), float(line[56:81])])\n")
        f.write("data = np.array(data)\n\n")

        f.write("x = data[:, 0]\n")
        f.write("y = data[:, 1]\n")
        f.write("z = data[:, 2]\n\n")

        f.write("point = 500\n")
        f.write("xi = np.linspace(min(x), max(x), point)\n")
        f.write("yi = np.linspace(min(y), max(y), point)\n")
        f.write("X, Y = np.meshgrid(xi, yi)\n")
        f.write("Z = griddata((x, y), z, (X, Y), method='linear')\n\n")

        f.write("fig = plt.figure(figsize=(10, 10))\n")
        f.write(f"plt.pcolormesh(X, Y, Z, cmap='{sel_palette}', shading='auto')\n")
        f.write("cbar = plt.colorbar(extend='both')\n")
        f.write("cbar.formatter.set_useMathText(True)\n")
        f.write("cbar.update_ticks()\n")
        f.write(f"plt.xlabel('{dim_1} [Å]')\n")
        f.write(f"plt.ylabel('{dim_2} [Å]')\n")
        f.write("plt.savefig(output_png, dpi=300)\n")
        f.write("plt.savefig(output, format='svg')\n")
        f.write("plt.close()\n")

def plot_2d_map(csv_file, plane, what): 

    output = csv_file[:-4] + ".png"
    dir_1 = plane[0]
    dir_2 = plane[1]

    if (what == "field") :
        sel_palette = "jet"
    elif (what == "density") :
        sel_palette = "seismic"

    # Load the data, skipping blank lines and considering the specific format
    data = []
    with open(csv_file, 'r') as file:
        for line in file:
            if line.strip():  # Ignore blank lines
                data.append([float(line[:25]), float(line[28:53]), float(line[56:81])])

    data = np.array(data)

    # Separate the data into x, y, z components
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    # Prepare the grid
    point = 300
    xi = np.linspace(min(x), max(x), point)
    yi = np.linspace(min(y), max(y), point)
    X, Y = np.meshgrid(xi, yi)

    # Interpolate the data
    Z = griddata((x, y), z, (X, Y), method='linear')


    # Determine the aspect ratio of the data
    x_range = max(x) - min(x)
    y_range = max(y) - min(y)
    aspect_ratio = x_range / y_range

    # Set figure size based on aspect ratio
    base_size  = 5
    fig_width  = base_size * aspect_ratio
    fig_height = base_size
    fig = plt.figure(figsize=(fig_width, fig_height))

    plt.pcolormesh(X, Y, Z, cmap=sel_palette, shading='auto')

    # Color bar
    cbar = plt.colorbar(extend='both')
    tick_font_size = 14
    cbar.ax.tick_params(labelsize=tick_font_size)
    cbar.formatter.set_powerlimits((0, 0))
    cbar.formatter.set_useMathText(True)
    cbar.ax.yaxis.set_offset_position('left')
    cbar.update_ticks()

    # Set axis limits and labels
    x_min, x_max = min(x), max(x)
    y_min, y_max = min(y), max(y)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.xlabel(dir_1+r' [$\AA$]', fontsize=20)
    plt.ylabel(dir_2+r' [$\AA$]', fontsize=20)
    plt.tick_params(axis='x', length=3, width=1, direction="out", bottom=True, labelbottom=True, labelsize=18)
    plt.tick_params(axis='y', length=3, width=1, direction="out", labelsize=18)

    # Adjust plot and save to file
    fig.subplots_adjust(left=0.07, bottom=0.11, right=1.05, top=0.95, wspace=0.0, hspace=0.0)
    plt.savefig(output, dpi=300)
    plt.close(fig)

def assign_polarization_name(field_dir):
    
    n_pol = 0
    pol_name = ["E", "E", "E"]
    if(field_dir == "all") :
       n_pol = 3
       pol_name[0] = "Ex"
       pol_name[1] = "Ey"
       pol_name[2] = "Ez"
    elif(field_dir == "x"):
       n_pol = 1
       pol_name[0] = "Ex"
    elif(field_dir == "y"):
       n_pol = 1
       pol_name[0] = "Ey"
    elif(field_dir == "z"):
       n_pol = 1
       pol_name[0] = "Ez"

    return n_pol, pol_name

def plot_planes_maps(args, frequencies): 

    root_filename = args.input_file[:-7]
    root_folder = "post_process_"+root_filename+"/planes/"
    folder_planes = root_folder + args.plane +"/"

    n_pol, pol_name = assign_polarization_name(args.field_dir)

    for i in range(args.num_ex_freq): 
        form_freq = f"{frequencies[i]:12.10f}".strip()
        for j in range(n_pol): 
            for k in range(args.n_plane):
                coord_plane = args.start + args.step_plane*k
                form_coord_plane = f"{coord_plane:20.2f}".strip()
                if(args.what == "field"): 
                   #plot 2D
                   csv_file = folder_planes+root_filename+"-"+pol_name[j]+"-"+form_freq+"-p-"+form_coord_plane+".csv"
                   plot_2d_map(csv_file,args.plane,args.what)
                   #create Python
                   file_python = folder_planes+root_filename+"-"+pol_name[j]+"-"+form_freq+"-p-"+form_coord_plane+".py"
                   csv_file_python = root_filename+"-"+pol_name[j]+"-"+form_freq+"-p-"+form_coord_plane+".csv"
                   write_plot_script(file_python,csv_file_python,args.plane, args.what)
                elif(args.what == "density"):
                   what1 = "-densityRe"
                   what2 = "-densityIm"
                   #plot 2D Re
                   csv_file = folder_planes+root_filename+"-"+pol_name[j]+"-"+form_freq+"-p-"+form_coord_plane+what1+".csv"
                   plot_2d_map(csv_file,args.plane,args.what)
                   #create python Re
                   file_python = folder_planes+root_filename+"-"+pol_name[j]+"-"+form_freq+"-p-"+form_coord_plane+what1+".py"
                   csv_file_python = root_filename+"-"+pol_name[j]+"-"+form_freq+"-p-"+form_coord_plane+what1+".csv"
                   write_plot_script(file_python, csv_file_python,args.plane, args.what)
                   #plot 2D Im
                   csv_file = folder_planes+root_filename+"-"+pol_name[j]+"-"+form_freq+"-p-"+form_coord_plane+what2+".csv"
                   plot_2d_map(csv_file,args.plane,args.what)
                   #create python Im
                   file_python = folder_planes+root_filename+"-"+pol_name[j]+"-"+form_freq+"-p-"+form_coord_plane+what2+".py"
                   csv_file_python = root_filename+"-"+pol_name[j]+"-"+form_freq+"-p-"+form_coord_plane+what2+".csv"
                   write_plot_script(file_python, csv_file_python, args.plane, args.what)
