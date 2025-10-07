#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import sys
import time
from datetime import datetime

def main():
    parser = argparse.ArgumentParser(description="Run plasmonX program")
    parser.add_argument("-i", "--input_file", type=str, required=True,
                        help="input file [.yaml file]")
    parser.add_argument("-o", "--output_file", type=str, required=False,
                        help="Optional: output file [Default: input_file with log extension")
    parser.add_argument("-omp", "--omp_threads", type=int, required=False,
                        help="Optional: number of OMP threads to exploit [Default: the maximum available]")
    parser.add_argument("-mem", "--memory_available", type=float, required=False,
                        help="Optional: available memory in GB [Default: the maximum available]")
    args = parser.parse_args()

    build_path = "@BUILD_PATH@"    #this is modified by CMake
    project_root = "@PROJECT_ROOT@"#this is modified by CMake

    # Add the parent directory of the project to PYTHONPATH
    parent_dir = os.path.abspath(project_root+"/python/")
    if parent_dir not in sys.path:
        sys.path.insert(0, parent_dir)

    #import my project
    from input_ import yaml_to_fortran_input, get_memory_in_gb, check_input_file, run_fortran_code
    from input_ import start_timers, print_execution_summary

    # Start timers
    start_cpu_time , start_wall_time = start_timers()

    #check if the file is readable
    check_input_file(args.input_file)

    #define output file
    if args.output_file:
        output_ = args.output_file
        # Check if the output file has a .log or .out extension
        _, ext = os.path.splitext(output_)
        if ext.lower() not in ('.log', '.out'):        
            print(f"Error: The output file '{output_}' must have a '.log' or '.out' extension.", file=sys.stderr)
            sys.exit(1)
    else: 
        output_ = args.input_file[:-5]+".log" #default


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

    fortran_file = args.input_file[:-5]+".tmp"

    citations = []

    #create the fortran input file
    yaml_to_fortran_input(start_cpu_time, start_wall_time, args.input_file, fortran_file, output_, n_omp_threads, available_GB, project_root, configurations, citations)

    #run the calculation 
    stdout, stderr, success, fortran_cpu_time = run_fortran_code(build_path + "/plasmonX", fortran_file)
    
    #final output
    print_execution_summary('plasmonX', start_cpu_time, start_wall_time, fortran_cpu_time, success, stderr, output_, citations=citations)

    #rm file fortran input file
    if os.path.exists(fortran_file):
        os.remove(fortran_file)

    # Redirect stdout and stderr
    sys.stdout.write(stdout)
    sys.stderr.write(stderr)

if __name__ == "__main__":
    main()

