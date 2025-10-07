#!/usr/bin/env python

import yaml
import os
from pathlib import Path
import sys
import re
import platform
import subprocess
import psutil
import time
from ase.io import read
from .validation import validate_input
from .resource_management import check_best_algorithm
from .timing_utils import print_execution_summary
from .geom_interface import geom_generation
from .check_sections import create_starting_keywords, check_unknown_section, check_section
from .references import define_references_to_be_printed
from collections import OrderedDict

def load_yaml_with_duplicate_check(path):

    class DuplicateKeyLoader(yaml.SafeLoader):
        duplicate_keys = []

    def construct_mapping(loader, node, deep=False):
        mapping = OrderedDict()
        for key_node, value_node in node.value:
            key = loader.construct_object(key_node, deep=deep)
            if key in mapping:
                DuplicateKeyLoader.duplicate_keys.append(f"Duplicate section found: '{key}'")
            value = loader.construct_object(value_node, deep=deep)
            mapping[key] = value
        return mapping

    DuplicateKeyLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping
    )

    with open(path, "r") as f:
        DuplicateKeyLoader.duplicate_keys.clear()
        data = yaml.load(f, Loader=DuplicateKeyLoader)

    return data, DuplicateKeyLoader.duplicate_keys

def convert_to_lowercase(data, exclude_keys=None):
    """
    Recursively converts all string keys and values in a nested dictionary or list to lowercase,
    except for specified keys.

    Args:
        data (dict | list): The data structure to convert.
        exclude_keys (set | None): Keys to exclude from lowercase conversion.

    Returns:
        dict | list: The data structure with all strings converted to lowercase, except for excluded keys.
    """
    if exclude_keys is None:
        exclude_keys = set()

    if isinstance(data, dict):
        return {
            key.lower(): (
                convert_to_lowercase(value, exclude_keys)
                if isinstance(value, (dict, list)) and key.lower() not in exclude_keys
                else value
            ) if key.lower() in exclude_keys else (
                convert_to_lowercase(value, exclude_keys)
                if isinstance(value, (dict, list)) else value.lower()
                if isinstance(value, str) else value
            )
            for key, value in data.items()
        }
    elif isinstance(data, list):
        return [
            convert_to_lowercase(item, exclude_keys)
            if isinstance(item, (dict, list)) else item.lower()
            if isinstance(item, str) else item
            for item in data
        ]
    return data

def read_geometry(data, yaml_file, keywords):
    """
    Read and process geometry from the 'input_geometry' section.

    Parameters:
        data (dict): Input data containing geometry information.
        yaml_file (str): Path to the input YAML file.
        keywords (dict): Dictionary of expected section keywords for validation.
    """
    input_geometry = data.get("input_geometry", {})

    processed_geometry = []
    unique_atomtypes = set()
    natoms, nmol = 0, 0
    errors = []
    string_input = ""

    # check the section types
    errors.extend(check_section("input_geometry", input_geometry, keywords))
    if errors != []:
        return processed_geometry, unique_atomtypes, natoms, nmol, string_input, errors

    ## READING XYZ
    if isinstance(input_geometry, dict) and "external xyz file" in input_geometry:
        xyz_file = input_geometry["external xyz file"]
        if not os.path.isfile(xyz_file):
            errors.append(f"Input geometry file '{xyz_file}' does not exist.")
        else:
            atoms = read(xyz_file)
            natoms = len(atoms)
            nmol = 1  # Default for single-molecule systems
            for atom in atoms:
                processed_geometry.append((atom.symbol, 1, atom.position[0], atom.position[1], atom.position[2]))
                unique_atomtypes.add(atom.symbol)
            string_input = "read from "+xyz_file

    ## READING INPUT FILE
    elif isinstance(input_geometry, str):
        geometry_lines = input_geometry.strip().split("\n")
        imol_pattern = re.compile(r"^(.*?)\s*\[\s*IMol\s*=\s*(\d+)\s*\]\s*(.*)$", re.IGNORECASE)
        natoms = len(geometry_lines)
        nmol_set = set()
        max_imol = 0
        for i, line in enumerate(geometry_lines, start=1):
            line = line.strip()
            match = imol_pattern.match(line)

            if match:
                atomtype = match.group(1).strip()
                imol_raw = match.group(2)
                coords = match.group(3).strip()
                if not imol_raw.isdigit() or int(imol_raw) <= 0:
                    errors.append(f"Line {i}: Invalid IMol value '{imol_raw}'. Must be a positive integer.")
                    continue
        
                imol = int(imol_raw)
                coord_parts = coords.split()
                if len(coord_parts) != 3:
                    errors.append(f"Line {i}: IMol structure found but expected 3 coordinates, got {len(coord_parts)}.")
                    continue
            else:
                parts = line.split()
                if len(parts) == 4:
                    atomtype, x, y, z = parts
                    imol = 1
                elif len(parts) > 4:
                    # Let see if there are IMol-like structure
                    has_bracket_open = "[" in line
                    has_bracket_close = "]" in line
                    has_equal = "=" in line
                
                    if not (has_bracket_open or has_bracket_close or has_equal):
                        errors.append(f"Line {i}: Too many elements and invalid format for atomtype '{parts[0]}'")
                        continue
                
                    # Let see which [ = ] are we missing
                    missing = []
                    if not has_bracket_open: missing.append("'['")
                    if not has_equal: missing.append("'='")
                    if not has_bracket_close: missing.append("']'")
                    if missing:
                        errors.append(f"Line {i}: Invalid IMol format. Missing {', '.join(missing)}")
                        continue
                
                    try:
                        between_brackets = line[line.index("[")+1:line.index("]")]
                        key_value = between_brackets.split("=")
                        if len(key_value) != 2:
                            errors.append(f"Line {i}: Failed IMol assignment. Expected format '[IMol=X]'")
                            continue
                
                        key, value = key_value[0].strip().lower(), key_value[1].strip()
                        if key != "imol":
                            errors.append(f"Line {i}: Expected 'IMol' as key in IMol info, got '{key}'")
                            continue
                
                        if not value.isdigit() or int(value) <= 0:
                            errors.append(f"Line {i}: IMol value must be a positive integer, got '{value}'")
                            continue
                
                    except Exception as e:
                        errors.append(f"Line {i}: Could not parse IMol info. Error: {str(e)}")
                        continue
                
                    # Here, IMol is correct but there are too many things in the line
                    errors.append(f"Line {i}: IMol info may be valid but too many values present in line. Expected format '[IMol=X] x y z'")
                    continue
                else:
                    errors.append(f"Line {i}: Not enough elements in geometry line.")
                    continue
                coord_parts = [x, y, z]
            try:
                xx, yy, zz = map(float, coord_parts)
            except ValueError:
                errors.append(f"Line {i}: Coordinates must be numeric: {coord_parts}")
                continue
        
            atomtypes = atomtype.capitalize() if len(atomtype) <= 2 else atomtype.upper()
            processed_geometry.append((atomtypes, imol, xx, yy, zz))
            unique_atomtypes.add(atomtypes)
            nmol_set.add(imol)
        
        nmol = len(nmol_set)
    
        # Check if all molecule indices are present
        max_imol = max(max_imol, imol)
        missing_mols = sorted(set(range(1, max_imol + 1)) - nmol_set)
        if missing_mols:
            errors.append(f"The following molecule indices are missing in the input geometry: {missing_mols}")        
        string_input = "read from yaml file"

    ## CREATING GEOMETRY WITH GEOM
    elif isinstance(input_geometry, dict):
        processed_geometry, unique_atomtypes, natoms, nmol, string_input, errors = geom_generation(input_geometry, yaml_file)

    else:
        errors.append("Invalid format for input_geometry section.")

    return processed_geometry, unique_atomtypes, natoms, nmol, string_input, errors

def initialize_output_file(start_cpu_time, start_wall_time, output_file, configurations, errors):
    """
    Initialize the output file with formatted header information.

    Parameters:
        start_cpu_time (float): Starting CPU time.
        start_wall_time (float): Starting wall-clock time.
        output_file (str): Path to the output file.
        configurations (dict): Configuration dictionary.
        errors (list): List of errors to report, if any.
    """    
    sticks = "-" * 80  # 80 ----
    # calculate the length
    max_key_length = max(len(key) for key in configurations.keys())

    with open(output_file, "w") as f:  
        f.write(f" {sticks}\n")
        f.write("                         _                                __  __ \n")
        f.write("                   _ __ | | __ _ ___ _ __ ___   ___  _ __ \ \/ / \n")
        f.write("                  | '_ \| |/ _` / __| '_ ` _ \ / _ \| '_ \ \  /  \n")
        f.write("                  | |_) | | (_| \__ \ | | | | | (_) | | | |/  \  \n")
        f.write("                  | .__/|_|\__,_|___/_| |_| |_|\___/|_| |_/_/\_\ \n")
        f.write("                  |_|                                            ")
        f.write("\n\n")
        f.write(f" {sticks}\n")
        f.write(" Tommaso Giovannini    University of Rome Tor Vergata, Rome, Italy\n")
        f.write(" Chiara Cappelli       Scuola Normale Superiore, Pisa, Italy\n")
        f.write(" Stefano Corni         University of Padova, Padova, Italy\n")
        f.write(" Luca Bonatti          Scuola Normale Superiore, Pisa, Italy\n")
        f.write(" Pablo Grobas Illobre  Scuola Normale Superiore, Pisa, Italy\n")
        f.write(" Piero Lafiosca        Scuola Normale Superiore, Pisa, Italy\n")
        f.write(" Luca Nicoli           Scuola Normale Superiore, Pisa, Italy\n")
        f.write(f" {sticks}\n")
        for key, value in configurations.items():
            f.write(f" {key.ljust(max_key_length)} : {value}\n")
        f.write(f" {sticks}\n")
    if errors:
        print_execution_summary('plasmonX',start_cpu_time, start_wall_time, 0, False, "", output_file, errors=errors)
        sys.exit()

def yaml_to_fortran_input(start_cpu_time, start_wall_time, yaml_file, fortran_file, out_file, n_omp, memory, project_root, configurations, citations):
    """
    Process a YAML input file and generate a corresponding FORTRAN input file.

    Parameters:
        start_cpu_time (float): Starting CPU time.
        start_wall_time (float): Starting wall-clock time.
        yaml_file (str): Path to the YAML input file.
        fortran_file (str): Path to the FORTRAN output file.
        out_file (str): Path to the main output log file.
        n_omp (int): Number of OpenMP threads.
        memory (int): Available memory in MB.
        project_root (str): Root path of the project.
        configurations (dict): Configuration templates and settings.
    """

    errors = []

    data, errors_yaml = load_yaml_with_duplicate_check(yaml_file)
    errors.extend(errors_yaml)

    # Convert everything to lowercase, except for 'atom_types' and 'parameters'
    exclude_keys = {'external xyz file'}
    data = convert_to_lowercase(data, exclude_keys=exclude_keys)

    starting_keywords = create_starting_keywords()

    # Check the names of the sections
    errors.extend(check_unknown_section(data, starting_keywords))
    if errors == [] :

        #Input Geometry to get atomtypes
        if "input_geometry" in data:
            processed_geometry, atomtypes, natoms, nmol, str_geometry, errors = read_geometry(data, yaml_file, starting_keywords)
        else:  
            atomtypes = []
            natoms = 0
            nmol = 0
            processed_geometry = []

        errors_input, used_defaults = validate_input(yaml_file, data, atomtypes, project_root)
        errors.extend(errors_input)

    initialize_output_file(start_cpu_time, start_wall_time, out_file, configurations, errors)

    #change algorithm based on calculation
    errors.extend(check_best_algorithm(data, atomtypes, natoms, n_omp, memory))
    if errors:
        print_execution_summary('plasmonX',start_cpu_time, start_wall_time, 0, False, "", out_file, errors=errors)
        sys.exit()


    # These are the defaults
    fields = [
        ("output file", out_file),
        ("omp threads", n_omp),
        ("memory", memory),
        ("what", "energy"),
        ("algorithm method", "inversion"),
        ("algorithm parallel execution", "frequencies"),
        ("algorithm number of iterations", 0),
        ("algorithm gmres dimension", 0),
        ("algorithm tolerance", 0.0),
        ("algorithm adaptive tuning", "yes"),
        ("forcefield static", "none"),
        ("forcefield dynamic", "none"),
        ("forcefield kernel", "none"),
        ("field type", "none"),
        ("field rhs type", "field"),
        ("field field intensity", "0.0"),
        ("field nfreq", "0"),
        ("field min freq", "0.0"),
        ("field max freq", "0.0"),
        ("field step freq", "0.0"),
        ("field external freq", "0.0"),
        ("field polarization", "none"),
        ("control no info file", False),
        ("control principal axes", False),
        ("output verbose", "0"),
        ("output maxima analysis", "none"),
        ("bem mesh file", "none"),
        ("bem normal scalar factor", "0.0"),
        ("bem permittivity file", "none"),
        ("bem permittivity", "none"),
        ("bem green function", "none"),
        ("bem sphere radius", "0.0d0"),
        ("bem solvent", "none"),
        ("bem epsilon solvent", "0.0d0"),
        ("bem variant", "dpcm"),
        ("used defaults for atom_types or parameters", False),
        ("atom_types number", "0"),
        ("atom_types name", ""),
        ("atom_types chi", "0.0"),
        ("atom_types eta", "0.0"),
        ("atom_types alpha", "0.0"),
        ("atom_types rq", "0.0"),
        ("atom_types rmu", "0.0"),
        ("parameters atomtype name", "none"),
        ("parameters atomtype tau", "0.0"),
        ("parameters atomtype sigma0", "0.0"),
        ("parameters atomtype scaling sigma0-tau", "0.0"),
        ("parameters atomtype a_ij", "0.0"),
        ("parameters atomtype fermi function d", "0.0"),
        ("parameters atomtype fermi function s", "0.0"),
        ("parameters atomtype fermi energy", "0.0"),
        ("parameters atomtype wfqfmu file", "none"),
        ("parameters atomtype permittivity", "none"),
        ("parameters interaction fermi function d", "0.0"),
        ("parameters interaction fermi function s", "0.0"),
        ("number of atoms", "0"),
        ("number of molecules", "0"),
    ]
    # updated "used defaults for atom_types or parameters"
    fields = [
        (field, used_defaults if field == "used defaults for atom_types or parameters" else default_value)
        for field, default_value in fields
    ]

    citations.extend(define_references_to_be_printed(data, atomtypes))

    with open(fortran_file, "w") as file:
        # Write generic fields (not atom_types e parameters)
        for field, default_value in fields:
            if field.startswith("atom_types") or field.startswith("parameters") or field.startswith("number"):
                continue
            if field.startswith("control"):
                key = field.split(" ", 1)[1]  # Key without "control"
                control_list = data.get("control", [])  # List "control" from YAML
                value = True if key in control_list else default_value  # True if present
                file.write(f"{field}: {value}\n")
            else:
                section, key = field.split(" ", 1) if " " in field else (field, None)
                if section == "what":  
                    what_list = data.get("what", [])
                    if isinstance(what_list, list) and len(what_list) == 1:
                        value = what_list[0]  
                    else:
                        value = what_list if what_list else default_value
                elif key:
                    value = data.get(section, {}).get(key, default_value)
                else:
                    value = data.get(section, default_value)
                
                #write
                if isinstance(value, list):
                    value_str = ", ".join(map(str, value))
                    file.write(f"{field}: {value_str}\n")
                else:
                    file.write(f"{field}: {value}\n")

        # Atom_types
        atom_types = data.get("atom_types", None)
        if atom_types is None:
            num_atom_types = 0
            file.write(f"atom_types number: {num_atom_types}\n")
        else:
            num_atom_types = atom_types.get("number", 1)
            file.write(f"atom_types number: {num_atom_types}\n")
        
            # reordered atomtypes
            ordered_atomtypes = sorted(atomtypes)  
            atom_index = 1
            for atom in ordered_atomtypes:  
                properties = atom_types.get(atom, {})
                file.write(f"#{atom_index}\n")
                file.write(f"atom_types name: {atom}\n")
                for key, default_value in fields:
                    if key.startswith("atom_types") and key != "atom_types number" and key != "atom_types name":
                        property_key = key.split(" ", 1)[1]
                        value = properties.get(property_key, default_value)
                        file.write(f"{key}: {value}\n")
                atom_index += 1
        
        # parameters
        parameters = data.get("parameters", {})
        atom_index = 1
        
        # atomtypes ---> parameters order
        ordered_atomtypes = sorted(atomtypes)  
        for atomtype in ordered_atomtypes:  
            atomtype_key = f"atomtype {atomtype}"
            properties = parameters.get(atomtype_key, {})
            file.write(f"#{atom_index}\n")
            file.write(f"parameters atomtype name: {atomtype}\n")
        
            fermi_function_written = False
        
            for key, default_value in fields:
                if key.startswith("parameters atomtype") and key != "parameters atomtype name":
                    property_key = key.split(" ", 2)[2]
        
                    if property_key == "fermi function d" or property_key == "fermi function s":
                        # write fermi function once
                        if not fermi_function_written:
                            fermi_function = properties.get("fermi_function", {})
                            for subkey in ["d", "s"]:
                                value = fermi_function.get(subkey, default_value)
                                file.write(f"parameters atomtype fermi function {subkey}: {value}\n")
                            fermi_function_written = True  
                    elif property_key not in properties:  
                        file.write(f"{key}: {default_value}\n")
                    else:  
                        value = properties[property_key]
                        file.write(f"{key}: {value}\n")
            atom_index += 1
        
        # Order the interactions
        interaction_keys = sorted({
            f"{a}->{b}" for a in atomtypes for b in atomtypes if a != b
        })  
        for interaction in interaction_keys:
            interaction_name = f"interaction {interaction}"
            indices = interaction.split("->")
            properties = parameters.get(interaction_name, {})
            file.write(f"#{indices[0]}-{indices[1]}\n")
            file.write(f"parameters interaction name: {interaction}\n")
        
            
            fermi_function_written = False
        
            for key, default_value in fields:
                if key.startswith("parameters interaction") and key != "parameters interaction name":
                    property_key = key.split(" ", 2)[2]
        
                    if property_key == "fermi function d" or property_key == "fermi function s":
                        if not fermi_function_written:
                            fermi_function = properties.get("fermi_function", {})
                            for subkey in ["d", "s"]:
                                value = fermi_function.get(subkey, default_value)
                                file.write(f"parameters interaction fermi function {subkey}: {value}\n")
                            fermi_function_written = True  
                    elif property_key not in properties: 
                        file.write(f"{key}: {default_value}\n")
                    else:  
                        value = properties[property_key]
                        file.write(f"{key}: {value}\n")

        for field, default_value in fields:
            if field == "number of atoms":
                file.write(f"{field}: {natoms}\n")
            elif field == "number of molecules":
                file.write(f"{field}: {nmol}\n")
    
        #print the geometry at the end
        if "input_geometry" in data:
            # Write processed geometry
            file.write(f"# Input geometry\n")
            file.write(f"input_geometry: "+str_geometry+"\n")
            for atomtype, imol, xx, yy, zz in processed_geometry:
                file.write(f"{atomtype:6s} {imol:10d} {xx:25.16f} {yy:25.16f} {zz:25.16f}\n")

def check_input_file(in_file): 
    input_file = Path(in_file)
    if not input_file.exists():
        print(f"Error: Input file '{input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    if not input_file.is_file():
        print(f"Error: The path '{input_file}' is not a file.", file=sys.stderr)
        sys.exit(1)
    if not os.access(input_file, os.R_OK):
        print(f"Error: The file '{input_file}' is not readable.", file=sys.stderr)
        sys.exit(1)
    if input_file.suffix.lower() != ".yaml":
        print(f"Error: The input file '{input_file}' does not have a '.yaml' extension.", file=sys.stderr)
        sys.exit(1)

def run_fortran_code(executable, input_file):
    """
    Execute a Fortran program and check if it completes successfully.

    Parameters:
        executable (str): Path to the Fortran executable.
        input_file (str): Input file to pass to the executable.

    Returns:
        tuple: (stdout, stderr, success) — Standard output, standard error, and success status as a boolean.
    """    
    try:
        # Run fortran
        fortran_process = subprocess.Popen(
            [executable, input_file], 
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Usa psutil for Fortran
        ps_process = psutil.Process(fortran_process.pid)

        # CPU init
        cpu_time_fortran = 0.0

        # CPU time until not stopped
        while fortran_process.poll() is None: 
            try:
                # Sum CPU time from all threads
                cpu_time_fortran = ps_process.cpu_times().user + ps_process.cpu_times().system
            except psutil.NoSuchProcess:
                break  

            time.sleep(0.1)  # Pause to be sure

        # Take stdout and stderr
        stdout, stderr = fortran_process.communicate()

        if fortran_process.returncode == 0:
            return stdout, stderr, True, cpu_time_fortran

        # Errors
        error_message = f"Program failed with exit code {fortran_process.returncode}"
        if fortran_process.returncode == 139:
            error_message = "Segmentation fault detected. Check memory usage or allocations."
        elif fortran_process.returncode == 134:
            error_message = "Abort signal received. Possible memory allocation issue."

        return stdout, error_message + f"\nStderr from process:\n{stderr}", False, cpu_time_fortran
    
    except subprocess.CalledProcessError as e:
        stderr = f"Program failed with exit code {e.returncode}"
        stderr += f"\nStderr from process:\n{e.stderr}"
        return e.stdout, stderr, False, 0 

    except FileNotFoundError:
        return None, "Executable not found.", False, 0 
