#!/usr/bin/env python

import psutil
import os

def get_memory_in_gb():
    """
    Retrieve the total available system memory in gigabytes.

    Returns:
        float: Total system memory in GB.
    """
    total_memory = psutil.virtual_memory().total
    return total_memory / (1024 ** 3)  # from byte a GB

def get_tesserae_msh_file(msh_file):
    if not msh_file or msh_file == "none" or not os.path.isfile(msh_file):
        return 0
    count, in_elements = 0, False
    with open(msh_file, "r") as f:
        for line in f:
            s = line.strip()
            if s.startswith("$Elements"):
                in_elements = True
                continue
            if in_elements:
                if s.startswith("$EndElements"):
                    break
                if len(s.split()) == 8:
                    count += 1
    return count    

def nvar_from_static_ff(ff, natoms, nts):
    if ff == "fq":
        return natoms
    elif ff == "fqfmu":
        return 4*natoms
    else:
        if nts > 0:
            return nts  # 
        else:
           return 0

def memory_per_frequency_gb(nvar):
    # estimate memory per frequency
    # we estimate only for atomistic calculations
    if nvar == 0:
        return 0.0

    # Compute approximate required memory for inversion 
    # 1 matrix real and 1 complex
    matrix_gb = 6.0 * (nvar** 2) * 8.0 / (1024.0 ** 3)
    # 1 array complex
    vector_gb = 4.0 * nvar* 3.0 * 8.0 / (1024.0 ** 3)

    return matrix_gb + vector_gb

def check_best_algorithm(data, atomtypes, natoms, n_omp, memory):
    """
    Determines the optimal algorithm method based on the number of atoms (natoms),
    available memory, and the static force field type. Performs an early safety
    check on frequency-parallel memory usage even when adaptive tuning is off.

    Args:
        data (dict): Input data containing algorithm and force field configurations.
        atomtypes (list): List of atom types extracted from geometry.
        natoms (int): Number of atoms in the system.
        n_omp (int): Number of OpenMP threads available.
        memory (float): Available memory in GB.
    """

    algorithm = data.get("algorithm", {}) or {}
    forcefield = data.get("forcefield", {}) or {}
    field = data.get("field", {}) or {}
    bem = data.get("bem", {}) or {}

    errors = []

    # --- Heterogeneity check ---
    heterogeneous = False
    if natoms > 0:
        if len(atomtypes) > 1:
            heterogeneous = True
            if (algorithm.get("parallel execution") == "frequencies") and algorithm.get("adaptive tuning") == "no":
                errors.append("Heterogeneous calculation and parallel execution: frequencies not allowed")
                return errors  

    static_forcefield = (forcefield.get("static", "none") or "none").lower()
    dynamic_forcefield = (forcefield.get("dynamic", "none") or "none").lower()
    requested_parallel = algorithm.get("parallel execution")

    msh_file = (bem.get("mesh file", "none") or "none").lower()
    nts = get_tesserae_msh_file(msh_file)

    nvar = nvar_from_static_ff(static_forcefield, natoms, nts)

    # let us check the usage memory for wMM methods or BEM at maximum
    if nvar > 0:
        nfreq = int(field.get("nfreq", 0) or 0)

        concurrent_freq = 1
        concurrent_freq = max(1, min(nfreq, int(n_omp or 1)))

        per_freq_gb = memory_per_frequency_gb(nvar)
        required_gb_freq_parallel = per_freq_gb * concurrent_freq

        safe_budget_gb = float(memory or 0.0) 

        # do not overcome the available memory
        if requested_parallel == "frequencies" and per_freq_gb > 0.0:
            if required_gb_freq_parallel > safe_budget_gb:
                if algorithm.get("adaptive tuning") == "no":
                    errors.append(
                        (f"Insufficient memory for frequency-parallel execution\n"
                         f"   Estimated: {required_gb_freq_parallel:.2f} GB \n"
                         f"   Estimated: {per_freq_gb:.2f} GB / frequency\n"
                         f"   Available: {safe_budget_gb:.2f} GB \n"
                         f"Switch to 'matrix' parallelization or use the iterative method.")
                    )
                    return errors
        elif requested_parallel == "matrix" and per_freq_gb > 0.0:
            if per_freq_gb > safe_budget_gb:
                if algorithm.get("adaptive tuning") == "no":
                    errors.append(
                        (f"Insufficient memory for inversion algorithm\n"
                         f"   Estimated: {per_freq_gb:.2f} GB \n"
                         f"   Available: {safe_budget_gb:.2f} GB \n"
                         f"Switch to the iterative method.")
                    )
                    return errors

    # If adaptive tuning and we arrived here, we can come back
    if algorithm.get("adaptive tuning") == "no":
        return errors

    #grep the field intensity
    field_intensity = field.get("field intensity", None)

    # Case 1: BEM calculation 
    if nts > 0:
        if algorithm.get("method") not in ["inversion", "iterative"]:
            algorithm["method"] = "inversion"
            if required_gb_freq_parallel <= safe_budget_gb:
                algorithm["parallel execution"] = "frequencies" if data.get("field", {}).get("nfreq", 0) > n_omp else "matrix"
            else:
                if per_freq_gb <= safe_budget_gb :
                    algorithm["parallel execution"] = "matrix"
                else: 
                    errors.append(
                        (f"Insufficient memory for BEM execution\n"
                         f"   Estimated: {per_freq_gb:.2f} GB \n"
                         f"   Available: {safe_budget_gb:.2f} GB \n"
                         f"To run BEM, reduce the number of the tesserae")
                    )
            algorithm["number of iterations"] = 0
            algorithm["gmres dimension"] = 0
            algorithm["tolerance"] = 0.0
            algorithm["rmse convergence"] = False
        else:
            if (algorithm.get("method") == "inversion"):
                if required_gb_freq_parallel <= safe_budget_gb:
                    algorithm["parallel execution"] = "frequencies" if data.get("field", {}).get("nfreq", 0) > n_omp else "matrix"
                else:
                    if per_freq_gb <= safe_budget_gb :
                        algorithm["parallel execution"] = "matrix"
                    else : 
                        errors.append(
                            (f"Insufficient memory for BEM execution\n"
                             f"   Estimated: {per_freq_gb:.2f} GB \n"
                             f"   Available: {safe_budget_gb:.2f} GB \n"
                             f"To run BEM, reduce the number of the tesserae")
                        )
                algorithm["number of iterations"] = 0
                algorithm["gmres dimension"] = 0
                algorithm["tolerance"] = 0.0
                algorithm["rmse convergence"] = False
        return errors

    # Case 2: Atomistic calculation (natoms > 0)
    if dynamic_forcefield == "none":
        # Static or energy calculation, no action required
        return errors  # statico: nessuna azione

    # Now decide the algorithm method based on per_freq_gb
    if required_gb_freq_parallel <= safe_budget_gb:
        # Use inversion method
        if algorithm["method"] == "inversion":
            algorithm["parallel execution"] = (
                "matrix" if heterogeneous
                else "frequencies" if data.get("field", {}).get("nfreq", 0) > n_omp
                else "matrix"
            )
            return errors
        else:
            algorithm["method"] = "inversion"
            algorithm["parallel execution"] = (
                "matrix" if heterogeneous
                else "frequencies" if data.get("field", {}).get("nfreq", 0) > n_omp
                else "matrix"
            )
            algorithm["number of iterations"] = 0
            algorithm["gmres dimension"] = 0
            algorithm["tolerance"] = 0.0
            algorithm["rmse convergence"] = False
    else:
        #we can handle the inversion for single frequencies
        if per_freq_gb <= safe_budget_gb :
            algorithm["method"] = "inversion"
            algorithm["parallel execution"] = "matrix"
            algorithm["number of iterations"] = 0
            algorithm["gmres dimension"] = 0
            algorithm["tolerance"] = 0.0
            algorithm["rmse convergence"] = False
        #we cannot handle the inversion algorithm --> iterative on the fly
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
    return errors
