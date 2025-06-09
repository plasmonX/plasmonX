#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import os
import sys


if __name__ == "__main__":

    build_path = "@BUILD_PATH@" #this is modified by CMake
    project_root = "@PROJECT_ROOT@"#this is modified by CMake

    # Add the parent directory of the project to PYTHONPATH
    parent_dir = os.path.abspath(project_root+"/python/")
    if parent_dir not in sys.path:
        sys.path.insert(0, parent_dir)

    #import my project
    from analysis import read_command_line, write_parameters_to_file, plot_planes_maps, print_banner
    from input_ import start_timers, print_execution_summary
    from input_ import get_memory_in_gb, run_fortran_code

    # Start timers
    start_cpu_time , start_wall_time = start_timers()

    args, min_grid, max_grid, frequencies = read_command_line()

    #configurations
    configurations = {
        "Configuration date": "@CONFIGURATION_DATE@",
        "Git branch": "@GIT_BRANCH@",
        "Fortran compiler": "@FORTRAN_COMPILER@",
        "Lapack type": "@LAPACK_TYPE@",
        "Blas type": "@BLAS_TYPE@",
        "Integers 64bit": "@INT64_STATUS@",
        "OpenMP": "@OMP_STATUS@",
    }

    #define number of OMP threads
    omp_enabled = "@OMP_STATUS@" == "ON"
    if omp_enabled:
        if args.omp_threads:
            n_omp_threads = args.omp_threads
        else:
            n_omp_threads = os.cpu_count()  # Default to all available cores
    else:
        if args.omp_threads:
            if args.omp_threads > 1:
                print("Warning: OMP Threads > 1 provided in input, but OpenMP compilation is disabled.")
        n_omp_threads = 1  # Single thread as fallback when OMP is disabled

    #define the number of GB available
    if args.memory_available:
        available_GB = float(args.memory_available)
    else:
        available_GB = get_memory_in_gb()

    write_parameters_to_file(args, min_grid, max_grid, frequencies, n_omp_threads)

    root_filename = args.input_file[:-7]
    output_ = "post_process_"+root_filename+"/"+root_filename+"_analysis.log"

    os.makedirs(os.path.dirname(output_), exist_ok=True)

    print_banner(output_, configurations)

    stdout, stderr, success, fortran_cpu_time = run_fortran_code(build_path + "/plasmonX_analysis", "parameters.txt")

    if(args.plane != "null") : 
       plot_planes_maps(args, frequencies)

    # summary
    print_execution_summary('plasmonX analysis', start_cpu_time, start_wall_time, fortran_cpu_time, success, stderr, output_)

    # Redirect stdout and stderr
    sys.stdout.write(stdout)
    sys.stderr.write(stderr)
