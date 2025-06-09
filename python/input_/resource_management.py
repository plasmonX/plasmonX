#!/usr/bin/env python

import psutil

def get_memory_in_gb():
    """
    Retrieve the total available system memory in gigabytes.

    Returns:
        float: Total system memory in GB.
    """
    total_memory = psutil.virtual_memory().total
    return total_memory / (1024 ** 3)  # Converti da byte a GB

def check_best_algorithm(data, natoms, n_omp, memory):
    """
    Determines the optimal algorithm method based on the number of atoms (natoms),
    available memory, and the static force field type.

    Args:
        data (dict): Input data containing algorithm and force field configurations.
        natoms (int): Number of atoms in the system.
        n_omp (int): Number of OpenMP threads available.
        memory (float): Available memory in GB.
    """
    algorithm = data.get("algorithm", {})
    if algorithm["adaptive tuning"] == "no":
        return 

    # grep the field intensity
    field = data.get("field", {})
    if ("field intensity" in field):
        field_intensity = field["field intensity"]

    # Case 1: BEM calculation (natoms == 0)
    if natoms == 0:
        if algorithm.get("method") not in ["inversion", "iterative"]:
            algorithm["method"] = "inversion"
            algorithm["parallel execution"] = "frequencies" if data.get("field", {}).get("nfreq", 0) > n_omp else "matrix"
            algorithm["number of iterations"] = 0
            algorithm["gmres dimension"] = 0
            algorithm["tolerance"] = 0.0
            algorithm["rmse convergence"] = False
        else:
            if (algorithm.get("method") == "inversion"):
                algorithm["parallel execution"] = "frequencies" if data.get("field", {}).get("nfreq", 0) > n_omp else "matrix"
                algorithm["number of iterations"] = 0
                algorithm["gmres dimension"] = 0
                algorithm["tolerance"] = 0.0
                algorithm["rmse convergence"] = False
        return


    # Case 2: Atomistic calculation (natoms > 0)
    forcefield = data.get("forcefield", {})
    static_forcefield = forcefield.get("static", "none").lower()
    dynamic_forcefield = forcefield.get("dynamic", "none").lower()

    if dynamic_forcefield == "none":
        # Static or energy calculation, no action required
        return

    # Determine nvar based on the static force field type
    if static_forcefield in ["fq", "fq_pqeq"]:
        nvar = 1
    elif static_forcefield in ["fqfmu", "fqfmu_pqeq"]:
        nvar = 4

    # Compute approximate required memory for inversion 
    # 1 matrix real and 1 complex
    matrix_memory_gb = 3* ((nvar * natoms) ** 2) * 8.0 / (1024.0 ** 3)  # GB for matrices
    # 1 array complex
    vector_memory_gb = 2 * (nvar * natoms) * 3 * 8.0 / (1024.0 ** 3)  # GB for vectors
    total_memory_gb = matrix_memory_gb + vector_memory_gb

    # Decide the algorithm method based on memory availability
    if total_memory_gb <= memory/2.0:
        # Use inversion method
        if algorithm["method"] == "inversion":
            algorithm["parallel execution"] = "frequencies" if data.get("field", {}).get("nfreq", 0) > n_omp else "matrix"
            return
        else:
            algorithm["method"] = "inversion"
            algorithm["parallel execution"] = "frequencies" if data.get("field", {}).get("nfreq", 0) > n_omp else "matrix"
            algorithm["number of iterations"] = 0
            algorithm["gmres dimension"] = 0
            algorithm["tolerance"] = 0.0
            algorithm["rmse convergence"] = False
    else:
        if algorithm["method"] == "iterative on the fly":
            algorithm["parallel execution"] = "matrix"
            if not algorithm.get("number of iterations"):  # Se è 0, None o False
                algorithm["number of iterations"] = 1000
            
            if not algorithm.get("gmres dimension"):  # Se è 0, None o False
                algorithm["gmres dimension"] = 1000
            
            if "rmse convergence" not in algorithm:
                algorithm["rmse convergence"] = True
            if algorithm["tolerance"] > field_intensity/10000.0:
                algorithm["tolerance"] = field_intensity/10000.0
            if not algorithm.get("tolerance"):  # Se è 0.0, None o False
                algorithm["tolerance"] = field_intensity / 10000.0
        else:
            # Use iterative method
            algorithm["method"] = "iterative on the fly"
            algorithm["parallel execution"] = "matrix"
            if not algorithm.get("number of iterations"):  # Se è 0, None o False
                algorithm["number of iterations"] = 1000
            
            if not algorithm.get("gmres dimension"):  # Se è 0, None o False
                algorithm["gmres dimension"] = 1000
            
            if not algorithm.get("tolerance"):  # Se è 0.0, None o False
                algorithm["tolerance"] = field_intensity / 10000.0
            
            if "rmse convergence" not in algorithm:
                algorithm["rmse convergence"] = True
