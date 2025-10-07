#!/usr/bin/env python

import os
import math
import re
from .check_sections import create_keywords, check_section, check_unknown_section

def convert_freq_to_au(value, unit):
    """
    Function to convert input frequencies to atomic units (a.u.).

    Currently supporting: eV [default], nm, micron, um, cm⁻¹, THz
    """
    if unit is None or unit.lower() == "ev":
        # eV → a.u. : CODATA https://physics.nist.gov/cgi-bin/cuu/Value?hrev
        return value / 27.211386245981 
    elif unit.lower() == "nm":
        # nm → a.u. : nm -> cm-1 -> a.u.
        return (10000000.0/219474.6313632) / value 
    elif unit.lower() == "micron" or unit.lower() == "um":
        # micron → a.u. : micron -> cm-1 -> a.u.
        return (10000.0/219474.6313632) / value  
    elif unit.lower() in ["cm-1", "cm^-1"]:
        # cm⁻¹ → a.u. : wikipedia
        return value / 219474.63136320  
    elif unit.lower() == "thz":
        # THz → a.u. : wikipedia
        return value / 6579.683920502 
    elif unit.lower() == "au" or unit.lower() == "a.u.":
        return value 
    else:
        raise ValueError(f"Unsupported unit '{unit}'. Allowed: au (a.u.), eV, nm, micron (um), cm-1 (cm^-1), THz.")

def convert_atomtype_value(value, unit, param_type):
    """
    Convert a parameter value to atomic units (a.u.) based on its type.

    - 'chi', 'eta': eV → a.u.
    - 'alpha': cm³, nm³, Å³, SI → a.u.
    - 'rq', 'rmu': nm, Å → a.u.

    Parameters:
        value (float): Numerical value to convert.
        unit (str): Unit of the input value.
        param_type (str): Type of the parameter (e.g., 'chi', 'alpha').
    """    
    # chi e eta
    if param_type in ["chi", "eta"]:
        if unit is None or unit.lower() in ["au", "a.u."]: 
            return value 
        elif unit.lower() == "ev":
            # eV → a.u. : CODATA https://physics.nist.gov/cgi-bin/cuu/Value?hrev
            return value / 27.211386245981 
        else:
            raise ValueError(f"Unsupported unit '{unit}' for '{param_type}'. Allowed: eV, a.u..")

    # alpha
    elif param_type == "alpha":
        if unit is None or unit.lower() in ["au", "a.u."]:
            return value  # Già in a.u.
        elif unit.lower() in ["nm^3", "nm3"]:
            return value / (0.0529177210544**3)  # nm³ → a.u.
        elif unit.lower() in ["angstrom^3", "ang^3", "angstrom3", "ang3"]:
            return value / (0.529177210544**3)  # Å³ → a.u.
        else:
            raise ValueError(f"Unsupported unit '{unit}' for 'alpha'. Allowed: nm3, ang3, a.u.")

    # rq e rmu 
    elif param_type in ["rq", "rmu"]:
        if unit is None or unit.lower() in ["au", "a.u."]:
            return value
        elif unit.lower() in ["angstrom", "ang"]:
            return value / 0.529177210544  # Å → a.u.
        elif unit.lower() == "nm":
            return value /0.0529177210544  # nm → a.u.
        else:
            raise ValueError(f"Unsupported unit '{unit}' for '{param_type}'. Allowed: ang, nm, a.u.")

    return value

def convert_parameter_to_au(value, unit, param_type):
    """
    Convert parameters to atomic units (a.u.).

    - 'tau': s, fs → a.u.
    - 'sigma0': S/m → a.u.
    - 'RI': nm², Å² → a.u.
    - 'gamma': s⁻¹, fs⁻¹ → a.u.
    - 'w_p': any frequency unit → a.u.

    Parameters:
        value (float): Numerical value to convert.
        unit (str): Unit of the input value.
        param_type (str): Type of the parameter to convert.
    """    
    # tau (fs, s → a.u.)
    if param_type == "tau":
        if unit is None or unit.lower() in ["au", "a.u."]:
            return value
        elif unit.lower() == "fs":
            return value * 41.3413733351702  # fs → a.u.
        elif unit.lower() == "s":
            return value * 41.3413733351702 * 1e15  # s → a.u.
        else:
            raise ValueError(f"Unsupported unit '{unit}' for 'tau'. Allowed: fs, s.")

    # sigma0 (S/m → a.u.)
    elif param_type == "sigma0":
        if unit is None or unit.lower() == "s/m":
            return value * 0.00000021739  # S/m → a.u.
        elif unit in ["au", "a.u."]:
            return value
        else:
            raise ValueError(f"Unsupported unit '{unit}' for 'sigma0'. Allowed: S/m, a.u.")

    # Aij (nm², Å² → a.u.)
    elif param_type == "a_ij":
        if unit is None or unit.lower() in ["au", "a.u."]:
            return value
        elif unit.lower() in ["nm^2", "nm2"]:
            return value / (0.0529177210544 ** 2)  # nm² → a.u.
        elif unit.lower() in ["angstrom^2", "ang2", "ang^2", "angstrom2"]:
            return value / (0.529177210544 ** 2)  # Å² → a.u.
        else:
            raise ValueError(f"Unsupported unit '{unit}' for 'RI'. Allowed: nm2, ang2.")

    # gamma (s⁻¹, fs⁻¹ → a.u.)
    elif param_type == "gamma":
        if unit is None or unit.lower() in ["au", "a.u."]:
            return value
        elif unit.lower() == "fs^-1" or unit.lower() == "fs-1":
            return value / 41.3413733351702  # fs⁻¹ → a.u.
        elif unit.lower() == "s^-1" or unit.lower() == "s-1":
            return value / (41.3413733351702 * 1e15)  # s⁻¹ → a.u.
        else:
            raise ValueError(f"Unsupported unit '{unit}' for 'gamma'. Allowed: s-1, fs-1.")

    # w_p (frequenza → a.u.)
    elif param_type == "w_p":
        return convert_freq_to_au(value, unit)  

    return value

def convert_length_to_au(value, unit):
    """
    Convert a length value to atomic units (a.u.).

    Accepted units:
    - "angstrom" → Conversion: value / 0.529177210544
    - "nm" → Conversion: value / 0.0529177210544
    - "au", "a.u." → No conversion needed (already in a.u.)

    Parameters:
        value (float): Numeric length value.
        unit (str): Unit of measurement ("angstrom", "nm", "au").
    """    

    if unit is None or unit.lower() in ["angstrom", "ang"]:
        return value / 0.529177210544
    elif unit.lower() in ["au", "a.u."]:
        return value  
    elif unit.lower() in ["nm"]:
        return value / 0.0529177210544
    else:
        raise ValueError(f"Unsupported unit '{unit}' for length. Allowed: ang, nm, au.")


def convert_length_to_angstrom(value, unit):
    """
    Convert a length value to Angstroms.

    Accepted units:
    - "angstrom" → No conversion
    - "nm" → Conversion: value * 10
    - "au" → Conversion: value * 0.529177210544

    Parameters:
        value (float): Numeric length value.
        unit (str): Unit of measurement ("angstrom", "nm", "au").
    """    

    if unit is None or unit.lower() in ["angstrom", "ang"]:
        return value 
    elif unit.lower() in ["au", "a.u."]:
        return value * 0.529177210544  
    elif unit.lower() in ["nm"]:
        return value * 10.0
    else:
        raise ValueError(f"Unsupported unit '{unit}' for length. Allowed: ang, angstrom, nm, au, a.u.")

def parse_value_with_unit(value):
    """
    Extract the numeric value and unit from a string input.

    Parameters:
        value (str | int | float): Input value, either as a number or a string with a unit.
    """    
    if isinstance(value, (int, float)):  # If it is already a number, do nothing
        return value, None, ""

    if isinstance(value, str):  # If it is a string let's' separate
        match = re.match(r"([0-9.eE+-]+)\s*(\S*)", value)
        if match:
            number = float(match.group(1))  # Convert number in float
            unit = match.group(2) if match.group(2) else None  # Unit
            return number, unit, ""
        else: 
            errors = f"Invalid string: '{value}'. Do not contain a number and a unit"

    return value, None, errors  # If we cannot understand anything, go back

def validate_what(data, keywords):
    """
    Validates and sets defaults for the 'what' section based on provided data.

    Args:
        data (dict): Input data.

    Returns:
        tuple: A list of errors identified during validation and the updated 'what'.
    """
    # Initialize errors
    errors = []

    # Get or set default for 'what'
    what = data.setdefault("what", [])

    # check the section types
    errors.extend(check_section("what", what, keywords))
    if errors != []:
        return errors, what

    # 1) Default to 'energy' if 'what' is empty or missing
    if not what:
        what.append("energy")

    # 2) Validate 'what'
    elif len(what) > 2:
        errors.append("'what' can contain a maximum of two keywords.")
    elif len(what) == 2:
        # 2.1) Check if one of the keywords is 'restart'
        if "restart" not in what:
            errors.append("If two keywords are provided in 'what' section, one must be 'restart'.")
        # 2.2) Validate combination of 'restart' with 'dynamic response'
        elif "dynamic response" not in what:
            errors.append("Keyword 'restart' is only valid with 'dynamic response'. Not yet implemented for other cases.")
        # 2.3) Reorder to place the useful keyword first, 'restart' second
        if what[0] == "restart":
            what[0], what[1] = what[1], what[0]

    data["what"] = what.copy()
    return errors, what

def validate_algorithm(data, what, keywords):
    """
    Validates and sets defaults for the 'algorithm' section based on 'what'.

    Args:
        data (dict): Input data.
        what (str): The type of analysis ('energy', 'static response', 'dynamic response', etc.).

    Returns:
        list: Errors identified during validation.
    """
    errors = []

    algorithm = data.setdefault("algorithm", {})

    # check the section types
    errors.extend(check_section("algorithm", algorithm, keywords))
    if errors != []:
        return errors

    if algorithm:
        if "method" not in algorithm:
            errors.append("Section 'algorithm' must contain 'method'.")
        if algorithm["method"] in ["iterative","iterative on the fly"] and what in ["energy", "static response"]:
            errors.append(f"Method: '{algorithm['method']}' not yet implemented for what: {what}. Must be inversion.")
        if "parallel execution" not in algorithm:
            if what in ["energy", "static response"]:
                algorithm["parallel execution"] = "matrix"
            else:
                if "method" in algorithm:
                    if "iterative on the fly" in algorithm["method"]:
                        algorithm["parallel execution"] = "matrix"
                    else:
                        algorithm["parallel execution"] = "frequencies"
    else: #default inversion
        algorithm["method"] = "inversion"
        algorithm["parallel execution"] = "matrix"

    #set other variables defaults
    if algorithm ["method"] == "inversion":
        algorithm["number of iterations"] = 0
        algorithm["gmres dimension"] = 0
        algorithm["tolerance"] = 0.0
        algorithm["rmse convergence"] = False
    else:
        if "number of iterations" not in algorithm:
            algorithm["number of iterations"] = 1000
        if "gmres dimension" not in algorithm:
            algorithm["gmres dimension"] = 1000
        if "tolerance" not in algorithm:
            algorithm["tolerance"] = 1.0e-9
        else:
            # convert format
            tolerance_value = algorithm["tolerance"]
            if isinstance(tolerance_value, str):
                if "d" in tolerance_value:
                    tolerance_sub = tolerance_value.replace("d", "e")
                    try:
                        tolerance_sub = float(tolerance_sub)
                        algorithm["tolerance"] = tolerance_sub
                    except ValueError:
                        errors.append(f"Invalid tolerance value: {tolerance_value}.")
            elif not isinstance(tolerance_value, float):
                errors.append(f"Invalid tolerance value: {tolerance_value}. Must be a string or a float.")                
        if "rmse convergence" not in algorithm:
            algorithm["rmse convergence"] = True

    # logical --> string
    if "adaptive tuning" not in algorithm:
        algorithm["adaptive tuning"] = "yes" #default

    data["algorithm"]["method"]               = algorithm["method"]              
    data["algorithm"]["parallel execution"]   = algorithm["parallel execution"]  
    data["algorithm"]["number of iterations"] = algorithm["number of iterations"]
    data["algorithm"]["gmres dimension"]      = algorithm["gmres dimension"]     
    data["algorithm"]["tolerance"]            = algorithm["tolerance"]           
    data["algorithm"]["rmse convergence"]     = algorithm["rmse convergence"]    
    data["algorithm"]["adaptive tuning"]      = algorithm["adaptive tuning"]    

    return errors

def assign_epsilon_solvent(solvent):
    """
    Assigns the dielectric constant (epsilon) based on the solvent name.

    Args:
        solvent (str): The name of the solvent.

    Returns:
        float: The dielectric constant for the given solvent.
    """

    #taken from https://macro.lsu.edu/howto/solvents/refractive%20index.htm
    refractive_indices = {
        "vacuum": 1.0,
        "trifluoroacetic acid": 1.2850,
        "methanol": 1.3284,
        "water": 1.3330,
        "acetonitrile": 1.3441,
        "ethyl ether": 1.3524,
        "1,1,2-trichlorotrifluoroethane": 1.3557,
        "pentane": 1.3575,
        "acetone": 1.3586,
        "ethyl alcohol": 1.3614,
        "petroleum ether": 1.3650,
        "methyl t-butyl ether": 1.3689,
        "ethyl acetate": 1.3724,
        "hexane": 1.3749,
        "isopropyl alcohol": 1.3772,
        "methyl ethl ketone": 1.3788,
        "glyme": 1.3796,
        "n-propyl alcohol": 1.3856,
        "heptane": 1.3876,
        "methyl n-propyl ketone": 1.3901,
        "iso-octane": 1.3914,
        "n-butyl acetate": 1.3942,
        "methyl isobutyl ketone": 1.3957,
        "isobutyl alcohol": 1.3959,
        "n-butyl alcohol": 1.3993,
        "triethylamine": 1.4010,
        "n-butyl chloride": 1.4021,
        "2-methoxyethanol": 1.4021,
        "cyclopentane": 1.4064,
        "methyl isoamyl ketone": 1.4070,
        "tetrahydrofuran": 1.4072,
        "propylene carbonate": 1.4210,
        "1,4-dioxane": 1.4224,
        "dichloromethane": 1.4241,
        "cyclohexane": 1.4262,
        "n,n-dimethylformamide": 1.4305,
        "isopropyl myristrate": 1.4332,
        "dimethyl acetamide": 1.4384,
        "ethylene dichloride": 1.4448,
        "chloroform": 1.4458,
        "n-methylpyrrolidone": 1.4700,
        "dimethyl sulfoxide": 1.4793,
        "toluene": 1.4969,
        "o-xylene": 1.5054,
        "pyridine": 1.5102,
        "chlorobenzene": 1.5248,
        "o-dichlorobenzene": 1.5514,
        "1,2,4-trichlorobenzene": 1.5717
    }

    # Normalize solvent name (case insensitive)
    solvent_normalized = solvent.lower()
    if solvent_normalized in refractive_indices:
        epsilon = refractive_indices[solvent_normalized]**2
    else:
        epsilon = None  
        raise ValueError(f"Unknown solvent '{solvent}'. Please provide a valid solvent.")

    return epsilon


def validate_bem(yaml_file, data, what, atomtypes, keywords, project_root):
    """
    Validates and sets defaults for the 'bem' section based on 'what'.

    Args:
        yaml_file: YAML file
        data (dict): Input data.
        what (str): The type of analysis ('energy', 'static response', 'dynamic response', etc.).
        atomtypes (set): Unique atom types in the system.
        project_root (str) : path to src.

    Returns:
        list: Errors identified during validation.
    """
    errors = []

    # Get or set default for 'bem'
    bem = data.setdefault("bem", {})

    if not bem:
        return errors

    # check the section types
    errors.extend(check_section("bem", bem, keywords))
    if errors != []:
        return errors

    # 1) BEM is only valid for dynamic calculations
    if what in ["energy", "static response"]:
        errors.append("BEM can only be run for what: dynamic response.")
        return errors

    # 2.1) Check atomtypes
    if atomtypes:
        errors.append("Atomistic-Continuum multiscale calculations not yet implemented. 'atomtypes' must be empty for BEM calculations.")

    # 2.2) Validate keywords
    # mesh file
    msh_file = bem.get("mesh file", yaml_file.replace(".yaml", ".msh"))
    if not os.path.isfile(msh_file):
        errors.append(f"BEM: Mesh file '{msh_file}' does not exist")
    else:
        with open(msh_file, 'r') as f:
            for line in f:
                if line.strip() == "$MeshFormat":
                    version_line = next(f).strip()
                    parts = version_line.split()
                    if len(parts) >= 2:
                        version_major = parts[0].split(".")[0]
                        file_type = parts[1]
                        if version_major == "2" and file_type == "0":
                            break  # Valid format
                    errors.append(f"BEM: Mesh file '{msh_file}' has unsupported format '{version_line}'. Only ASCII-2 format is supported.")
                    break  # Stop checking after $MeshFormat section
    bem["mesh file"] = msh_file

    # normal scalar factor
    if "normal scalar factor" in bem:
        nsf = bem["normal scalar factor"]
        if isinstance(nsf, str) and "d" in nsf:
            nsf = float(nsf.replace("d", "e"))
            bem["normal scalar factor"] = nsf
        if not isinstance(nsf, float) or nsf not in [1.0, -1.0]:
            errors.append("BEM 'normal scalar factor' must be a float with value 1.0 or -1.0.")
    else:
        bem["normal scalar factor"] = 1.0

    # Ensure either 'permittivity' or 'permittivity file' is present, but not both
    if "permittivity" in bem and "permittivity file" in bem:
        errors.append("Only one of 'permittivity' or 'permittivity file' can be provided. Remove one.")
    elif "permittivity" not in bem and "permittivity file" not in bem:
        errors.append("Either 'permittivity' or 'permittivity file' must be provided for BEM calculation.")

    if "permittivity" in bem:
        path_to_csv = project_root+"/src/parameters/"
        if bem["permittivity"] in ["silver etchegoin","silver jc","silver johnson-christy","silver johnson christy"]: 
            bem["permittivity"] = "silver etchegoin"
        elif bem["permittivity"] in ["gold etchegoin","gold jc","gold johnson-christy","gold johnson christy"]: 
            bem["permittivity"] = "gold etchegoin"
        elif bem["permittivity"] in ["silver bb","silver brendel-bormann","silver brendel bormann"]:
            bem["permittivity"] = "none"
            bem["permittivity file"] = path_to_csv + "bb_permittivity_ag.csv"
        elif bem["permittivity"] in ["silver palik"]:
            bem["permittivity"] = "none"
            bem["permittivity file"] = path_to_csv + "palik_permittivity_ag.csv"
        elif bem["permittivity"] in ["gold bb","gold brendel-bormann","gold brendel bormann"]:
            bem["permittivity"] = "none"
            bem["permittivity file"] = path_to_csv + "bb_permittivity_au.csv"
        elif bem["permittivity"] in ["gold palik"]:
            bem["permittivity"] = "none"
            bem["permittivity file"] = path_to_csv + "palik_permittivity_au.csv"

    # permittivity file
    if "permittivity file" in bem:
        if not os.path.isfile(bem["permittivity file"]):
            errors.append(f"BEM permittivity file '{bem['permittivity file']}' does not exist.")
        errors.extend(check_freq_in_file(bem["permittivity file"], data))

    # green function
    if "green function" not in bem:
        bem["green function"] = "approximate"

    # variant
    if "variant" not in bem:
        bem["variant"] = "dpcm"

    # sphere radius
    if "sphere radius" not in bem and bem["green function"] == "approximate": #default: read msh_file and calculate it
        if os.path.isfile(msh_file):
            try:
                with open(msh_file, "r") as f:
                    lines = f.readlines()
                    
                # Trova la sezione $Nodes e il numero di nodi
                node_section = False
                min_coords = [float("inf")] * 3
                max_coords = [float("-inf")] * 3
                
                for line in lines:
                    if "$Nodes" in line:
                        node_section = True
                        continue
                    if "$EndNodes" in line:
                        break
                    if node_section:
                        parts = line.split()
                        if len(parts) >= 4:  # Deve contenere almeno 4 valori: ID e (x, y, z)
                            x, y, z = map(float, parts[1:4])
                            min_coords = [min(min_coords[i], coord) for i, coord in enumerate([x, y, z])]
                            max_coords = [max(max_coords[i], coord) for i, coord in enumerate([x, y, z])]
    
                # Calcola il massimo intervallo tra min e max nelle tre direzioni
                if float("inf") not in min_coords and float("-inf") not in max_coords:
                    max_diff = max(max_coords[i] - min_coords[i] for i in range(3)) #diameter
                    # we round to the largest integer, we divide by 2.0 -> radius, and we convert to au (assuming input in Ang)
                    bem["sphere radius"] = convert_length_to_au(float(round(max_diff, 1))/2.0, None) 
                else:
                    errors.append("BEM: Failed to parse node coordinates from mesh file.")
            except Exception as e:
                errors.append(f"BEM: Error reading mesh file '{msh_file}': {str(e)}")
    elif "sphere radius" in bem:
        if bem["green function"] == "approximate":
            if isinstance(bem["sphere radius"],str):
                value, unit, errors_parse = parse_value_with_unit(bem["sphere radius"])  # Default angstrom
                if errors_parse == "":
                    try:
                        bem["sphere radius"] = convert_length_to_au(value, unit)
                    except ValueError as e:
                        errors.append(str(e))            
                else:
                    errors.append(errors_parse + " (BEM section)")
            else:
                bem["sphere radius"] = convert_length_to_au(bem["sphere radius"], None)
        else: 
            errors.append(f"Sphere radius and accurate green function in BEM section is not allowed")

    # solvent
    if "solvent" not in bem:
        bem["solvent"] = "vacuum"
    try:
        # Assign the epsilon value based on the solvent
        bem["epsilon solvent"] = assign_epsilon_solvent(bem["solvent"])
    except ValueError as e:
        errors.append(str(e))        

    return errors

def validate_output(data, what, keywords):
    """
    Validates and sets defaults for the 'output' section based on 'what'.

    Args:
        data (dict): Input data.
        what (str): The type of analysis ('energy', 'static response', 'dynamic response', etc.).

    Returns:
        list: Errors identified during validation.
    """
    errors = []

    # output
    output = data.setdefault("output", {})

    # check the section types
    errors.extend(check_section("output", output, keywords))
    if errors != []:
        return errors

    #Handle verbose
    if "verbose" in output:
        if output["verbose"] < 0 : 
            errors.append(f"Invalid 'verbose': {output['verbose']}. Must be >= 0")

    # Handle maxima analysis validation and defaulting
    if "maxima analysis" in output:
        if (what in ["energy", "static response"]) and output["maxima analysis"] != "none":
            errors.append(f"Invalid 'maxima analysis' keyword: {output['maxima analysis']} for what: {what}. Must be 'none'.")
    else:
        # Set defaults based on 'what'
        if what in ["dynamic response"]:
            output["maxima analysis"] = "absorption"
        else:
            output["maxima analysis"] = "none"    
           
    data["output"]["maxima analysis"] = output["maxima analysis"]

    return errors

def validate_control(data, keywords):
    """
    Validates and sets defaults for the 'output' section based on 'what'.

    Args:
        data (dict): Input data.

    Returns:
        list: Errors identified during validation.
    """
    errors = []

    # control
    control = data.setdefault("control", {})

    # check the section types
    errors.extend(check_section("control", control, keywords))
    if errors != []:
        return errors

    return errors

def validate_field(data, what, keywords):
    """
    Validates the 'field' section based on the value of 'what'.

    Args:
        data (dict): Data for the 'field' section.
        what (str): Type of response ('static response' or 'dynamic response').

    Returns:
        list: Errors identified during validation.
    """
    errors = []

    #check if present
    field = data.setdefault("field", {})

    # check the section types
    errors.extend(check_section("field", field, keywords))
    if errors != []:
        return errors


    if field and what == "energy": #inconsistent keywords
        errors.append(f"Required {what} calculation, but field section present")
    else:
        if what != "energy" and not field: #inconsistent keywords
            if what != "-": #it means there is an error
                errors.append(f"Required {what} calculation, but not field section")

    unit_map = {}  # units init
    for key in ["min freq", "max freq", "step freq"]:
        if key in field: 
            value, unit, errors_parse = parse_value_with_unit(field[key])  
            if errors_parse == "":
                unit_map[key] = unit  # save unit 
                # Conversion in atomic units
                try:
                    field[key] = convert_freq_to_au(value, unit)
                except ValueError as e:
                    errors.append(str(e))  # error
            else:
                errors.append(errors_parse + "(field section)")
            if not isinstance(field[key], float):
                errors.append(f"'{key}' must be a float.")
    #swap min freq and max freq if nm in input
    if "min freq" in field and "max freq" in field:
        if unit_map.get("min freq") and unit_map.get("max freq"):  
            if unit_map["min freq"].lower() == "nm" and unit_map["max freq"].lower() == "nm":
                if field["min freq"] > field["max freq"]:
                    field["min freq"], field["max freq"] = field["max freq"], field["min freq"]  # Swap
    if "external freq" in field:
        # external freq --> able to recognize very general string, e.g: "1.0, 400 nm, 800 THz, 2 eV"
        external_freq = field["external freq"]
    
        if isinstance(external_freq, str):
            try:
                # Split string by commas, then parse each value
                freq_list = [freq.strip() for freq in external_freq.split(",")]
                parsed_freqs = []
                for freq in freq_list:
                    value, unit, errors_parse = parse_value_with_unit(freq)  
                    if (errors_parse == ""):
                        # Conversion unit
                        try:
                            converted_value = convert_freq_to_au(value, unit)  # a.u.
                            parsed_freqs.append(converted_value)
                        except ValueError as e:
                            errors.append(f"Invalid unit in 'external freq': {e}")
                    else:
                        errors.append(errors_parse + " (field section)")
                #Update the external freq key
                field["external freq"] = parsed_freqs  
            except ValueError:
                errors.append("'external freq' contains invalid numbers.")
    
        elif isinstance(external_freq, float):
            # If it is a single float, we put on a string and convert from eV (default) to a.u.
            field["external freq"] = [convert_freq_to_au(external_freq, None)]
    
        elif isinstance(external_freq, list):
            converted_list = []
            # This is standard case (everything in eV)
            for freq in external_freq:
                if isinstance(freq, (int, float)):
                    converted_list.append(convert_freq_to_au(freq, None))  
                elif isinstance(freq, str):
                    value, unit, errors_parse = parse_value_with_unit(freq)
                    if errors_parse == "":
                        try:
                            converted_list.append(convert_freq_to_au(value, unit))
                        except ValueError as e:
                            errors.append(f"Invalid unit in 'external freq': {e}")
                    else:
                        errors.append(errors_parse + "(field section)")
                else:
                    errors.append(f"Invalid type in 'external freq': {freq}")
    
            field["external freq"] = converted_list
    
        else:
            errors.append("'external freq' must be a list of floats, a single float, or a comma-separated string of values.")

    if "polarization" not in field:
        if what in ["energy", "static response"]:
            field["polarization"] = "none"
        elif what == "dynamic response":
            field["polarization"] = "all"

    if what in ["energy", "static response"]:
        if field.get("polarization") != "none":
            errors.append("'polarization' must be 'none' for 'energy' or 'static response'.")

    #validate the rest
    if what == "static response":
        # 1.1) Ensure 'type' is 'static'
        if "type" in field and field["type"] != "static":
            errors.append("Field type must be 'static' for static response.")
        # 1.2) If 'type' is missing, set it to 'static'
        elif "type" not in field:
            field["type"] = "static"

        # 1.3) Error if forbidden parameters are present
        forbidden_keys = ["field intensity", "nfreq", "min freq", "max freq", "step freq", "external freq"]
        for key in forbidden_keys:
            if key in field:
                errors.append(f"'{key}' is not allowed for static response.")

        # 1.4) Set all unspecified parameters to 0
        default_keys = ["field intensity", "min freq", "max freq", "step freq", "external freq"]
        for key in default_keys:
            if key not in field:
                field[key] = 0.0

        # 1.4) Set all unspecified parameters to 0
        default_keys = ["nfreq"]
        for key in default_keys:
            if key not in field:
                field[key] = 0

    elif what == "dynamic response":
        # 2.1) Ensure 'type' is 'dynamic'
        if "type" in field and field["type"] != "dynamic":
            errors.append("Field type must be 'dynamic' for dynamic response.")
        # 2.2) If 'type' is missing, set it to 'dynamic'
        elif "type" not in field:
            field["type"] = "dynamic"

        if "field intensity" not in field:
            field["field intensity"] = 1.0e-4 #new default

        # 2.3) Handle 'external freq'
        if "external freq" in field:
            external_freq = field["external freq"]
            if isinstance(external_freq, str):
                try:
                    # Convert string to list of floats
                    external_freq = [float(freq.strip()) for freq in external_freq.split(",")]
                    field["external freq"] = external_freq
                except ValueError:
                    errors.append("'external freq' contains invalid numbers.")
            elif isinstance(external_freq, float):
                # Handle the case where a single float is provided
                field["external freq"] = [external_freq]
            elif not isinstance(external_freq, list):
                errors.append("'external freq' must be a list or a comma-separated string of frequencies.")

            #check if too small
            if any(freq < 1.0e-10 for freq in field["external freq"]):
                errors.append("'external freq' contains values that are too small (< 1.0e-10).")

            # Check for incompatibility with other frequency parameters
            conflicting_keys = ["min freq", "max freq", "step freq"]
            for key in conflicting_keys:
                if key in field:
                    errors.append(f"'external freq' is incompatible with '{key}'. Please use only one method to specify frequencies.")

            # Calculate 'nfreq' based on the list length if not provided
            if "nfreq" in field and len(field["external freq"]) != field["nfreq"]:
                errors.append(
                    f"'nfreq' ({field['nfreq']}) does not match the number of external frequencies provided ({len(field['external freq'])})."
                )
            else:
                field["nfreq"] = len(field["external freq"])

        # 2.4) Handle frequency ranges
        else:
            required_keys = ["min freq", "max freq", "step freq", "nfreq"]
            present_keys = [key for key in required_keys if key in field]
            if len(present_keys) < 3:
                errors.append("At least three of 'min freq', 'max freq', 'step freq', 'nfreq' must be provided.")
            else:
                # first check that min freq and max freq > 0.0
                if "min freq" in field and field["min freq"] < 1.0e-10:
                    errors.append("'min freq' is too small (< 1.0e-10).")
                if "max freq" in field and field["max freq"] < 1.0e-10:
                    errors.append("'max freq' is too small (< 1.0e-10).")
                
                #now check that everything is consistent
                if "min freq" in field and "max freq" in field and "step freq" in field:
                    if "nfreq" in field:
                        n_freq_input = field["nfreq"]
                    else:
                        n_freq_input = 0
                    field["nfreq"] = math.ceil((field["max freq"] - field["min freq"]) / field["step freq"]) + 1
                    if n_freq_input != 0 and n_freq_input != field["nfreq"]:
                        errors.append("'nfreq' given as input is not coherent with the other parameters")
                elif "min freq" in field and "step freq" in field and "nfreq" in field:
                    field["max freq"] = field["min freq"] + (field["nfreq"] - 1) * field["step freq"]
                elif "max freq" in field and "step freq" in field and "nfreq" in field:
                    field["min freq"] = field["max freq"] - (field["nfreq"] - 1) * field["step freq"]
                elif "nfreq" in field and "min freq" in field and "max freq" in field:
                    field["step freq"] = (field["max freq"] - field["min freq"]) / (field["nfreq"] - 1)

    data["field"] = field.copy()

    return errors

def validate_forcefield(data, what, keywords):
    """
    Validates the 'forcefield' section based on 'what' and 'field'.

    Args:
        data: Data for the 'forcefield' section.
        what (str): The type of response ('energy', 'static response', 'dynamic response', etc.).

    Returns:
        list: Errors identified during validation.
    """
    errors = []

    #check if present
    forcefield = data.setdefault("forcefield", {})
    if not forcefield:
        errors.append(f"Forcefield section is mandatory for atomistic calculations")
        return errors

    # check the section types
    errors.extend(check_section("forcefield", forcefield, keywords))
    if errors != []:
        return errors

    # Define defaults
    kernel_defaults = {
        "fq": "ohno",
        "fqfmu": "gaussian",
    }

    # Case 1: Static responses
    if what in ["energy", "static response"]:
        # 1.1) Ensure forcefield is defined as static
        if "dynamic" in forcefield:
            errors.append(f"'dynamic' forcefield is not allowed for '{what}' responses.")
        elif "static" not in forcefield:
            errors.append(f"'static' forcefield must be defined for '{what}' responses.")
        else:
            static_forcefield = forcefield["static"].lower()

        # 1.2) Handle kernel defaults
        kernel = forcefield.get("kernel", kernel_defaults.get(forcefield.get("static", ""), "ohno"))
        if forcefield.get("static") == "fq" and kernel not in ["ohno", "gaussian", "coulomb"]:
            errors.append(f"Invalid kernel for 'fq': {kernel}. Must be 'ohno', 'gaussian', or 'coulomb'.")
        elif forcefield.get("static") == "fqfmu" and kernel != "gaussian":
            errors.append(f"Invalid kernel for '{forcefield['static']}': {kernel}. Must be 'gaussian'.")
        else:
            forcefield["kernel"] = kernel  # Set the kernel default if missing

    # Case 2: Dynamic responses
    elif what in ["dynamic response"]:
        # 2.1) Ensure forcefield is defined as dynamic
        if "dynamic" not in forcefield:
            errors.append(f"'dynamic' forcefield must be defined for '{what}' responses.")
            dynamic_forcefield = ""
        else:
            dynamic_forcefield = forcefield["dynamic"].lower()

        # 2.2) Handle optional static forcefield
        if "static" in forcefield:
            static_forcefield = forcefield["static"].lower()
            if dynamic_forcefield == "wfq" and static_forcefield != "fq":
                errors.append(f"Invalid 'static' forcefield for 'wfq': {static_forcefield}. Must be 'fq'.")
            elif dynamic_forcefield == "wfqfmu" and static_forcefield != "fqfmu":
                errors.append(f"Invalid 'static' forcefield for 'wfqfmu': {static_forcefield}. Must be 'fqfmu'.")
        else:
            # 2.3) Set defaults for static forcefield if missing
            if dynamic_forcefield == "wfq":
                forcefield["static"] = "fq"
            elif dynamic_forcefield == "wfqfmu":
                forcefield["static"] = "fqfmu"

        # 2.4) Handle kernel defaults
        kernel = forcefield.get("kernel")
        if kernel is None:
            forcefield["kernel"] = "gaussian"
        elif kernel.lower() != "gaussian":
            errors.append(f"Invalid kernel for 'wfq' or 'wfqfmu': {kernel}. Must be 'gaussian'.")        

    data["forcefield"] = forcefield.copy()

    return errors


def validate_atomtypes(data, unique_atomtypes, keywords):
    """
    Validate atom types from the geometry and assign default values if needed.

    Parameters:
        data (dict): Input data.
        unique_atomtypes (set): Set of unique atom types.
        keywords (dict): Reference keyword definitions.
    """    
    errors = []
    used_defaults = False

    predefined_atomtypes = {
        "Au": {
            "chi": 0.0,
            "eta": 0.52685,
            "alpha": 39.5297,
            "rq": 3.4996572071713126,
            "rmu": 3.0710227258450864,
        },
        "Ag": {
            "chi": 0.0,
            "eta": 0.4829838,
            "alpha": 49.9843,
            "rq": 3.3775793621738925,
            "rmu": 2.368922250085755,
        },
        "C": {
            "chi": 0.0,
            "eta": 0.372124,
            "rq": 0.0,
        },
        "Na": {
            "chi": 0.0,
            "eta": 0.292000,
            "rq": 0.0,
        },
    }

    # Required keys depending on forcefield if something is given in input
    required_keys_by_forcefield = {
        "fq": ["chi", "eta"],
        "fqfmu": ["chi", "eta", "alpha"],
        "wfq": ["eta"],
        "wfqfmu": ["eta", "alpha"],
    }

    # Allowed keys depending on static forcefield
    allowed_keys_by_forcefield = {
        "fq": ["chi", "eta", "rq"],
        "fqfmu": ["chi", "eta", "alpha", "rq", "rmu"]
    }

    atom_types = data.setdefault("atom_types", {})

    # check the section types
    errors.extend(check_section("atom_types", atom_types, keywords))
    if errors != []:
        return errors, used_defaults

    # Normalize keys in atom_types based on length
    atom_types = {
        (key.lower().capitalize() if len(key) <= 2 else key.upper()): value
        for key, value in atom_types.items()
    }
    lowercase_keys = {key.lower(): key for key in atom_types.keys()}
    lowercase_keys = {
        (key.lower().capitalize() if len(key) <= 2 else key.upper()): key
        for key in atom_types.keys()
    }

    for atomtype in unique_atomtypes:
        if atomtype not in lowercase_keys:
            if atomtype in predefined_atomtypes:
                atom_types[atomtype] = predefined_atomtypes[atomtype].copy()
                used_defaults = True
            else:
                errors.append(f"Atomtype parameters for atomtype '{atomtype}' are missing and no predefined values exist.")
        # Here we convert units if needed
        else: 
            for param, value in atom_types[atomtype].items():
                parsed_value, unit, errors_parse = parse_value_with_unit(value)
                if (errors_parse == ""):
                    try:
                        atom_types[atomtype][param] = convert_atomtype_value(parsed_value, unit, param)
                    except ValueError as e:
                        errors.append(str(e))  
                else:
                    errors.append(errors_parse + "(forcefield section)")

    # forcefield check input
    forcefield = data.get("forcefield", {})
    if "dynamic" in forcefield: 
        static_forcefield = forcefield["static"].lower()
        dynamic_forcefield = forcefield["dynamic"].lower()
        required_keys = required_keys_by_forcefield.get(dynamic_forcefield, [])
        ff_name = dynamic_forcefield
    else: 
        static_forcefield = forcefield["static"].lower()
        required_keys = required_keys_by_forcefield.get(static_forcefield, [])
        ff_name = static_forcefield

    allowed_keys = allowed_keys_by_forcefield.get(static_forcefield, [])

    # Determine required keys based on the static forcefield

    for atomtype_key, properties in atom_types.items():
        atomtype = (
            atomtype_key.capitalize() if len(atomtype_key) <= 2 else atomtype_key.upper()
        )
        for required_key in required_keys:
            if required_key not in properties:
                errors.append(f"Missing required key '{required_key}' for atomtype '{atomtype}' with forcefield '{ff_name}'.")
                if ff_name == "wfq" or ff_name == "wfqfmu":
                    errors.append("If you want to use the default values for wFQ or wFQFmu models, delete the entire 'atom_types' section.")        
            elif required_key != "chi" and properties[required_key] <= 0:
                errors.append(f"Value of '{required_key}' for atomtype '{atomtype}' must be greater than 0.")
        for key in properties: 
            if key not in allowed_keys: 
                errors.append(f"Key '{key}' is not allowed for atomtype '{atomtype}' with forcefield '{ff_name}'.")
        # Set unnecessary keys to 0.0 for consistency
        for key in predefined_atomtypes.get(atomtype, {}):
            if key not in required_keys and key not in properties:
                atom_types[atomtype][key] = 0.0

    # Additional calculations if kernel is gaussian
    kernel = data.get("forcefield", {}).get("kernel", "gaussian").lower()
    if kernel == "gaussian":
        for atomtype_key, properties in atom_types.items():
            eta = properties.get("eta", 0.0)
            alpha = properties.get("alpha", 0.0)
            r_q = properties.get("rq", 0.0)
            r_mu = properties.get("rmu", 0.0)

            if eta <= 0:
                errors.append(f"'eta' must be defined and greater than 0 for atomtype '{atomtype_key}' with kernel 'gaussian'.")
                continue
            if static_forcefield == "fqfmu" and alpha <= 0:
                errors.append(f"'alpha' must be defined and greater than 0 for atomtype '{atomtype_key}' with kernel 'gaussian'.")
                continue

            # Calculate r_q and r_mu based on the forcefield
            if static_forcefield == "fq":
                if r_q == 0.0: #define rq based on limit 
                    properties["rq"] = (2 / math.pi)**0.5 / eta

            if static_forcefield == "fqfmu":
                if r_q == 0.0: #define rq based on limit
                    properties["rq"] = (2 / math.pi)**0.5 / eta
                if r_mu == 0.0: #define rmu based on limit
                    properties["rmu"] = (alpha / (3.0 * math.sqrt(math.pi /2.0)) )**(1.0 / 3)

    #check R_q and R_mu gt 0
    for atomtype_key, properties in atom_types.items():
        if "rq" in properties:
            if properties["rq"] < 0.0 : 
                errors.append(f"'R_q less than 0 for atomtype '{atomtype_key}'.")
        if "rmu" in properties:
            if properties["rmu"] < 0.0 : 
                errors.append(f"'R_mu less than 0 for atomtype '{atomtype_key}'.")
       
    for atomtype in atom_types.keys():
        if atomtype not in unique_atomtypes:
            errors.append(f"Parameters defined for unused atomtype '{atomtype}'.")

    atom_types["number"] = len(unique_atomtypes)  # Usa la chiave corretta per 'number'

    # Update the data dictionary to ensure atom_types is saved back properly
    data["atom_types"] = atom_types

    return errors, used_defaults


def validate_parameters(data, unique_atomtypes, keywords):
    """
    Validate that each atom type in the geometry has properly defined parameters.

    Parameters:
        data (dict): Main input dictionary.
        unique_atomtypes (set): Set of atom types present in the geometry.
        keywords (dict): Reference keyword definitions.
    """    
    errors = []
    used_defaults = False

    predefined_values_by_atomtype = {
        "Au": { #ACS Photonics, 2022, 9, 3025−3034
            "tau": 318.018,
            "sigma0": 11834849.81,
            "a_ij": 9.61,
            "fermi_function": {"d": 12.0, "s": 1.1},
            "permittivity": "gold etchegoin",
        },
        "Ag": { #ACS Photonics, 2022, 9, 3025−3034
            "tau": 1633.608,
            "sigma0": 6.5313e7,
            "a_ij": 9.61,
            "fermi_function": {"d": 12.0, "s": 1.1},
            "permittivity": "silver etchegoin",
        },
        "C": { #J. Phys. Chem. Lett. 2020, 11, 7595−7602
            "tau": 170.0,
            "a_ij": 1.7424,
            "fermi_function": {"d": 100.0, "s": 1.2},
            "fermi energy": 0.4,
        },
        "Na": { #Nanoscale, 2019, 11, 6004–6015
            "tau": 1323.0,
            "sigma0": 2.39494681e7,
            "scaling sigma0-tau": 10.0,
            "a_ij": 12.07910025,
            "fermi_function": {"d": 12.0, "s": 1.1},
        },
    }

    # Required keys depending on forcefield if something is given in input
    required_keys_by_forcefield = {
        "wfq": {
            "two_of": ["tau", "sigma0", "w_p", "gamma"],
            "requires_group": {
                "fermi_function": ["d", "s"]
            },
            "required": {"a_ij"},
            "conditional": {
                "C": ["fermi energy"]
            }
        },
        "wfqfmu": {
            "two_of": ["tau", "sigma0", "w_p", "gamma"],
            "required": {"a_ij"},
            "requires_group": {
                "fermi_function": ["d", "s"]
            },
            "one_of": ["permittivity", "wfqfmu file"]
        }
    }
    
    parameters = data.setdefault("parameters", {})

    # check the section types
    errors.extend(check_section("parameters", parameters, keywords))
    if errors != []:
        return errors, used_defaults

    # Normalize keys in parameters based on length
    parameters = {
        (key if key.lower().startswith("interaction")
         else (key.lower().capitalize() if len(key) <= 2 else key.upper())): value
        for key, value in parameters.items()
    }
    

    # forcefield check input
    forcefield = data.get("forcefield", {})
    if "dynamic" in forcefield: 
        dynamic_forcefield = forcefield["dynamic"].lower()
    else: 
        errors.append("Parameters section given by no dynamic forcefield")
        return errors, used_defaults

    for atomtype in unique_atomtypes:
        predefined_values = predefined_values_by_atomtype.get(atomtype, {})

        # Check if the atomtype exists in parameters, considering case-insensitivity
        parameter_key = atomtype if atomtype in parameters else None
        if not parameter_key:
            if not predefined_values:
                errors.append(f"Parameters for atomtype '{atomtype}' are missing and no predefined values exist.")
            else:
                parameters[f"{atomtype}"] = predefined_values.copy()
                used_defaults = True
        else:
            param_data = parameters[parameter_key]

            rules = required_keys_by_forcefield.get(dynamic_forcefield, {})

            # 1. check the autoexclusive rules
            keys = rules.get("two_of", [])
            present_keys = [k for k in keys if k in param_data]
            if atomtype == "C":
                if len(present_keys) != 1:
                    errors.append(
                        f"Atomtype '{atomtype}' in forcefield '{dynamic_forcefield}' requires only one between 'tau' and 'gamma', found {present_keys}"
                    )
            else:
                if len(present_keys) != 2:
                    errors.append(
                        f"Atomtype '{atomtype}' in forcefield '{dynamic_forcefield}' requires at least 2 of {keys}, found {present_keys}"
                    )
            # 2. mandatory keys
            for key in rules.get("required", []):
                if key not in param_data:
                    errors.append(f"Missing required key '{key}' for atomtype '{atomtype}' and forcefield '{dynamic_forcefield}'")
            # 3. One of 
            one_of = rules.get("one_of", [])
            found = [k for k in one_of if k in param_data]
            if one_of and len(found) != 1:
                errors.append(f"One of {one_of} must be specified for atomtype '{atomtype}', found {found}")
            # 4. mandatory groups (fermi_function)
            for group_key, subkeys in rules.get("requires_group", {}).items():
                group = param_data.get(group_key, {})
                if not isinstance(group, dict):
                    errors.append(f"Expected group '{group_key}' to be a dictionary in atomtype '{atomtype}'")
                    continue
                for subkey in subkeys:
                    if subkey not in group:
                        errors.append(f"Missing '{subkey}' in group '{group_key}' for atomtype '{atomtype}'")
            # 5. If Graphene, needed fermi energy
            conditional = rules.get("conditional", {})
            if atomtype in conditional:
                for key in conditional[atomtype]:
                    if key not in param_data:
                        errors.append(f"Missing conditionally required key '{key}' for atomtype '{atomtype}'")
        
    # Validate parameter types and defaults
    for atomtype in unique_atomtypes:
        if atomtype in parameters:
            param = parameters[atomtype]

            if dynamic_forcefield == "wfq":
                for key in param:
                    if key in ["permittivity", "wfqfmu file"]:
                        errors.append(f"Invalid key '{key}' for atomtype '{atomtype}' and forcefield 'wfq'")
                    if atomtype == "C":
                        if key in ["sigma0", "w_p"]:
                            errors.append(f"Invalid key '{key}' for atomtype 'C' and forcefield 'wfq'")            

            # Convert tau, sigma0, gamma, w_p, RI in a.u.
            for key in ["tau", "sigma0", "gamma", "w_p", "a_ij"]:
                if key in param:
                    value, unit, errors_parse = parse_value_with_unit(param[key])
                    if errors_parse == "": 
                        try:
                            param[key] = convert_parameter_to_au(value, unit, key)
                        except ValueError as e:
                            errors.append(str(e))
                    else:
                        errors.append(errors_parse + "(parameters section)")
    
            # sigma0 and tau
            has_tau = "tau" in param
            has_sigma0 = "sigma0" in param
            has_gamma = "gamma" in param
            has_w_p = "w_p" in param
            # if atomtype ("C"), we check only tau
            if atomtype == "C":
                if not (has_tau or has_gamma):
                    errors.append(f"Missing required 'tau' or 'gamma' for '{atomtype}'.")
                elif has_gamma and not has_tau:
                    param["tau"] = 1 / param["gamma"]  # calculate tau
            else:
                if has_sigma0 and has_tau:
                    pass  
                elif has_w_p and has_tau:
                    param["sigma0"] = (param["w_p"] ** 2 * param["tau"]) / (4 * math.pi)  
                elif has_w_p and has_gamma:
                    param["tau"] = 1 / param["gamma"]
                    param["sigma0"] = (param["w_p"] ** 2 * param["tau"]) / (4 * math.pi) 
                elif has_sigma0 and has_gamma:
                    param["tau"] = 1 / param["gamma"]  
                else:
                    errors.append(
                        f"Invalid parameter combination in '{atomtype}'. "
                        "You must provide either (sigma0, tau), (w_p, tau), (w_p, gamma), or (sigma0, gamma)."
                    )            

            param.pop("gamma", None)
            param.pop("w_p", None)

            # Ensure wfqfmu file and permittivity are mutually exclusive
            if "wfqfmu file" in param and param["wfqfmu file"] != "none":
                # If 'wfqfmu file' is defined, set 'permittivity' to 'none'
                param["permittivity"] = "none"
            elif "permittivity" not in param or param["permittivity"] == "none":
                # Default case: set 'permittivity' to its default value
                param["permittivity"] = predefined_values.get("permittivity", "none")
            # Check if 'wfqfmu file' exists
            wfqfmu_file = param.get("wfqfmu file", "none")
            if wfqfmu_file != "none":
                if not os.path.isfile(wfqfmu_file):
                    errors.append(f"'wfqfmu file' {wfqfmu_file} in atomtype '{atomtype}' section does not exist.")
                else:
                    errors.extend(check_freq_in_file(wfqfmu_file, data))

            if param["tau"] != 0.0:
                scaling = param.get("scaling sigma0-tau", 0.0)
                if scaling == 0.0:
                    param["scaling sigma0-tau"] = 1.0

    # Rename atomtypes, excluding those starting with "interaction"
    normalized_parameters = {
        (f"atomtype {atomtype}" if not atomtype.lower().startswith("interaction") else atomtype): values
        for atomtype, values in parameters.items()
    }
    
    
    # Update the data dictionary to ensure parameters is saved back properly
    data["parameters"] = normalized_parameters

    return errors, used_defaults

def check_freq_in_file(file_read, data):
    """
    Validates that the frequencies in the data match those defined in the wfqfmu file.

    Args:
        file_read (str): Path to the file containing frequency data.
        data (dict): Input data containing frequency definitions.

    Returns:
        list: Errors encountered during validation.
    """
    errors = []

    # Extract frequencies from data
    field = data.get("field", {})
    freq = []

    if "nfreq" in field and field["nfreq"] == 0:
       return

    if "external freq" in field and field["external freq"] != 0.0:
        external_freq = field["external freq"]
        if not isinstance(external_freq, list):
            errors.append("'external freq' must be a list of values.")
        else:
            freq = external_freq
    else:
        min_freq = field.get("min freq", 0.0)
        step_freq = field.get("step freq", 0.0)
        nfreq = field.get("nfreq", 0)
        freq = [min_freq + i * step_freq for i in range(nfreq)]

    # Open and parse the wfqfmu file
    try:
        with open(file_read, "r") as file:
            freq_read = []
            for line in file:
                line = line.strip()
                if line.startswith("#") or not line:  # Skip comments and empty lines
                    continue
                parts = line.split()
                if len(parts) != 3:
                    errors.append(f"Invalid format in {file_read}: {line}. Expected 3 columns of float values.")
                    continue
                try:
                    freq_value = float(parts[0])
                    # Validate the remaining columns
                    float(parts[1])  # Just to ensure they're valid floats
                    float(parts[2])
                    freq_read.append(freq_value)
                except ValueError:
                    errors.append(f"Invalid numeric values in {file_read}: {line}.")

    except IOError as e:
        errors.append(f"Error reading the file '{file_read}': {str(e)}.")
        return errors

    tolerance = 1e-12  # Tolerance for comparing permittivity file and input freq
    
    for f in freq:
        freq_eV = f * 27.211386245981
        if not any(abs(freq_eV - fr) < tolerance for fr in freq_read):
            errors.append(f"Frequency {freq_eV:.12f} not found in {file_read} [tolerance 10^-12].")    

    return errors


def validate_interactions(data, unique_atomtypes):
    """
    Validate interactions between atom types and assign default values if needed.

    Parameters:
        data (dict): Input data.
        unique_atomtypes (set): Set of unique atom types.
    """
    errors = []
    used_defaults = False

    # predefined interactions for Ag->Au and Au->Ag
    predefined_interactions = {
        "Ag->Au": {"fermi_function": {"d": 12.0, "s": 0.95}},
        "Au->Ag": {"fermi_function": {"d": 12.0, "s": 0.95}},
    }

    #these are the required interactions for input geometry
    required_interactions = {
        f"{a}->{b}" for a in unique_atomtypes for b in unique_atomtypes if a != b
    }

    parameters = data.setdefault("parameters", {})

    lowercase_keys = {key.lower(): key for key in parameters.keys()}  # Map lowercase keys to original keys

    normalized_parameters = {}
    for key, value in parameters.items():
        if key.startswith("interaction"):
            interaction = key.split(" ", 1)[1].lower()
            # Normalize based on the length of atomtypes in the interaction
            parts = interaction.split("->")
            normalized_interaction = "->".join(
                part.capitalize() if len(part) <= 2 else part.upper()
                for part in parts
            )
            normalized_parameters[f"interaction {normalized_interaction}"] = value    

    parameters.update(normalized_parameters)

    defined_interactions = {
        key.split()[1] for key in parameters.keys() if key.startswith("interaction")
    }

    # Extract fermi_function values from atomic parameters
    atomic_fermi = {}
    for atom_key, values in parameters.items():
        if atom_key.startswith("atomtype "):
            atom = atom_key.split(" ", 1)[1]  # Get "Au" from "atomtype Au"
            fermi = values.get("fermi_function", {})
            d = fermi.get("d")
            s = fermi.get("s")
            if d is not None and s is not None:
                atomic_fermi[atom] = {"d": d, "s": s}

    for interaction in required_interactions:
        reverse_interaction = "->".join(reversed(interaction.split("->")))

        # 1: Both interactions are missing
        if f"interaction {interaction}" not in parameters and f"interaction {reverse_interaction}" not in parameters:
            if interaction in predefined_interactions and reverse_interaction in predefined_interactions: 
                # we take the default values
                parameters[f"interaction {interaction}"] = {
                    "fermi_function": {"d": predefined_interactions[interaction]["fermi_function"]["d"], 
                                       "s": predefined_interactions[interaction]["fermi_function"]["s"]}
                }
                parameters[f"interaction {reverse_interaction}"] = {
                    "fermi_function": {"d": predefined_interactions[reverse_interaction]["fermi_function"]["d"], 
                                       "s": predefined_interactions[reverse_interaction]["fermi_function"]["s"]}
                }
                used_defaults = True
            else: 
                # case 1, we have both d and s for both atomtypes defined in parameters -> arithmetic mean
                if interaction.split("->")[0] in atomic_fermi and interaction.split("->")[1] in atomic_fermi:
                    atom1, atom2 = interaction.split("->")
                    d1, s1 = atomic_fermi[atom1]["d"], atomic_fermi[atom1]["s"]
                    d2, s2 = atomic_fermi[atom2]["d"], atomic_fermi[atom2]["s"]
                    avg_d = (d1 + d2) / 2
                    avg_s = (s1 + s2) / 2
                
                    parameters[f"interaction {interaction}"] = {
                        "fermi_function": {"d": avg_d, "s": avg_s}
                    }
                    parameters[f"interaction {reverse_interaction}"] = {
                        "fermi_function": {"d": avg_d, "s": avg_s}
                    }                
                else:
                    # case 2, error -> we do not known how to calculate the d and s parameters
                    errors.append(f"Missing interaction: {interaction} and {reverse_interaction}")

        # 2: Only one interaction is missing
        elif f"interaction {interaction}" in parameters and f"interaction {reverse_interaction}" not in parameters:
            parameters[f"interaction {reverse_interaction}"] = {
                key: (value.copy() if isinstance(value, dict) else value)
                for key, value in parameters[f"interaction {interaction}"].items()
            }
            used_defaults = True

        elif f"interaction {reverse_interaction}" in parameters and f"interaction {interaction}" not in parameters:
            parameters[f"interaction {interaction}"] = {
                key: (value.copy() if isinstance(value, dict) else value)
                for key, value in parameters[f"interaction {reverse_interaction}"].items()
            }
            used_defaults = True    

    # Update the data dictionary to ensure parameters is saved back properly
    data["parameters"] = parameters

    return errors, used_defaults

def validate_atomtypes_and_parameters(data, unique_atomtypes, what, keywords):
    """
    Main function to validate atom types and parameters based on the calculation type.

    Parameters:
        data (dict): Input data.
        unique_atomtypes (set): Set of atom types in the geometry.
        what (str): Type of calculation ("dynamic response", "energy", etc.).
        keywords (dict): Reference keyword definitions.
    """

    # atomtypes
    errors_atomtypes, used_defaults_atomtypes = validate_atomtypes(data, unique_atomtypes, keywords)

    # parameters and interactions
    parameters = data.get("parameters", {})
    if what == "dynamic response":
        # parameters
        errors_parameters, used_defaults_parameters = validate_parameters(data, unique_atomtypes, keywords)
        # interactions
        errors_interactions, used_defaults_interactions = validate_interactions(data, unique_atomtypes)
    elif parameters and (what == "energy" or what == "static response"): #static calculation
        errors.append(f"Required {what} calculation, but paramters section present")
        errors_parameters = []
        errors_interactions = []
        used_defaults_parameters = False
        used_defaults_interactions = False
    else: 
        errors_parameters = []
        errors_interactions = []
        used_defaults_parameters = False
        used_defaults_interactions = False

    used_defaults = used_defaults_atomtypes or used_defaults_parameters or used_defaults_interactions

    errors = []
    errors.extend(errors_atomtypes)
    errors.extend(errors_parameters)
    errors.extend(errors_interactions)

    return errors, used_defaults

def validate_input(yaml_file, data, atomtypes, project_root):
    """
    Validate and prepare data for generating the Fortran input file.

    Parameters:
        yaml_file (str): Path to the input YAML file.
        data (dict): Parsed input data.
        atomtypes (list): List of atom types extracted from geometry.
        project_root (str): Root directory of the project.
    """    
    errors = []

    #create 
    updated_keywords = create_keywords(data, atomtypes)

    # Check the names of the sections
    errors.extend(check_unknown_section(data, updated_keywords))
    if errors != []:
        return errors, False

    # what section
    errors_what, what = validate_what(data, updated_keywords)
    errors.extend(errors_what)

    #algorithm section
    errors.extend(validate_algorithm(data, what[0], updated_keywords))
    # field
    errors.extend(validate_field(data, what[0], updated_keywords))

    if atomtypes:
        # forcefield
        errors_forcefield = validate_forcefield(data, what[0], updated_keywords)
        if errors_forcefield != []:
            errors.extend(errors_forcefield)
            used_defaults = False
        else:
            # atomtypes, parameters
            errors_atomtypes_parameters, used_defaults = validate_atomtypes_and_parameters(data, atomtypes, what[0], updated_keywords)
            errors.extend(errors_atomtypes_parameters)
    else:
        used_defaults = False

    # bem
    errors.extend(validate_bem(yaml_file, data, what[0],atomtypes, updated_keywords, project_root))

    # control
    errors.extend(validate_control(data, updated_keywords))

    # output
    errors.extend(validate_output(data, what[0], updated_keywords))

    return errors, used_defaults
