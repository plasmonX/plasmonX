#!/usr/bin/env python

import sys
import os
import subprocess
import shutil
from ase.io import read
from .validation import parse_value_with_unit, convert_length_to_angstrom

VALID_KEYS = {
    "type": str,  # implicit or atomistic
    "shape": str,
    "atomtype": str,
    "radius": (float,str),  # sphere
    "main_axis": str,  # rod main axis
    "length": (float,str),  # length [rod, pyramid]
    "width": (float,str),  # width rod
    "core_atomtype": str,  # alloy
    "core_radius": (float,str),  # alloy sphere
    "core_length": (float,str),  # alloy rod
    "core_width": (float,str),  # alloy rod
    "shell_atomtype": str,  # alloy
    "shell_radius": (float,str),  # alloy sphere ; shell radius for rings
    "shell_length": (float,str),  # alloy sphere
    "shell_width": (float,str),  # alloy sphere
    "base_length": (float,str),  # pyramid
    "side_length": (float,str),  # triangle graphene
    "height": (float,str),  # height [pyramid, cone, tip, microscope]
    "mesh_size": (float,str),  # bem mesh
    "distance": (float,str),  # dimer -- this means that this is dimer
    "edge_type": str,  # graphene edge
    "a": float,  # parameters for paraboloid controlla lungo X [tip]
    "b": float,  # parameters for paraboloid controlla lungo Y [tip]
    "dimer_axis": str,    # axis for the dimer creation
    "radius_in": (float,str),   # internal radius for rings
    "height_tip": (float,str),  # only used if microscope -- height of the small pyramid on top of paraboloid
    "dimer": str,         # geometry for dimer "bowtie" or "default"
}

SHAPE_REQUIRED_PARAMS = {
    "sphere": {"radius"},  # Can be either single atomtype or core-shell
    "rod": {"main_axis", "length", "width"},  # Can be either single atomtype or core-shell
    "pyramid": {"base_length", "height", "atomtype"},
    "icosahedron": {"radius", "atomtype"},
    "ico": {"radius", "atomtype"},
    "ih": {"radius", "atomtype"},
    "cuboctahedron": {"radius", "atomtype"},
    "cto": {"radius", "atomtype"},
    "decahedron": {"radius", "atomtype"},
    "idh": {"radius", "atomtype"},
    "cone": {"radius", "height", "atomtype"},
    "microscope": {"a", "b", "height", "base_length", "height_tip", "atomtype"},
    "paraboloid": {"a", "b", "height", "atomtype"},
    "tip": {"a", "b", "height", "atomtype"},
    "ribbon": {"atomtype", "length", "width", "edge_type"},
    "disk": {"radius", "atomtype"},
    "ring": {"radius", "radius_in", "atomtype"},  
    "triangle": {"side_length", "edge_type", "atomtype"},
}

# Definition of valid params by shape
SHAPE_VALID_PARAMS = {
     
    "sphere": {"radius"     , "atomtype" ,         # Simple atomistic sphere
               "type"       , "mesh_size",         # Simple continuum sphere
               "core_radius", "core_atomtype", "shell_radius", "shell_atomtype", # Atomistic core-shell sphere
               "dimer", "dimer_axis", "distance"}, # Dimer
    
    "rod": {"main_axis", "length", "width", "atomtype",  # Simple atomistic rod
            "type"     , "mesh_size",                    # Simple continuum rod
            "core_length", "core_width", "core_atomtype", "shell_length", "shell_width", "shell_atomtype", # Atomistic core-shell rod 
            "dimer", "dimer_axis", "distance"},    # Dimer

    "pyramid": {"base_length", "height", "atomtype", 
                "dimer", "dimer_axis", "distance"},  # Dimer + bowtie

    "icosahedron": {"radius", "atomtype", 
                    "dimer", "dimer_axis", "distance"},  # Dimer

    "ico": {"radius", "atomtype",
            "dimer", "dimer_axis", "distance"}, # Dimer

    "ih": {"radius", "atomtype",
           "dimer", "dimer_axis", "distance"},  # Dimer

    "cuboctahedron": {"radius", "atomtype",
                      "dimer", "dimer_axis", "distance"},  # Dimer

    "cto": {"radius", "atomtype",
            "dimer", "dimer_axis", "distance"},  # Dimer

    "decahedron": {"radius", "atomtype",
                   "dimer", "dimer_axis", "distance"},  # Dimer

    "idh": {"radius", "atomtype",
            "dimer", "dimer_axis", "distance"},  # Dimer

    "cone": {"radius", "height", "atomtype",
             "dimer", "dimer_axis", "distance"},  # Dimer + bowtie

    "microscope": {"a", "b", "height", "base_length", "height_tip", "atomtype",
                   "dimer", "dimer_axis", "distance"},  # Dimer + bowtie

    "paraboloid": {"a", "b", "height", "atomtype",
                   "dimer", "dimer_axis", "distance"},  # Dimer + bowtie

    "tip": {"a", "b", "height", "atomtype",
            "dimer", "dimer_axis", "distance"},  # Dimer + bowtie

    "ribbon": {"atomtype", "length", "width", "edge_type"},

    "disk": {"atomtype", "radius"},

    "ring": {"atomtype", "radius", "radius_in"}, 

    "triangle": {"atomtype", "side_length", "edge_type"},
}

SHAPE_MAPPING = {
    "sphere": "-sphere",
    "rod": "-rod",
    "pyramid": "-pyramid",
    "icosahedron": "-ico",
    "ico": "-ico",
    "ih": "-ico",
    "cuboctahedron": "-cto",
    "cto": "-cto",
    "decahedron": "-idh",
    "idh": "-idh",
    "paraboloid": "-tip",
    "tip": "-tip",
    "cone": "-cone",
    "ribbon": "-graphene rib",
    "disk": "-graphene disk",
    "ring": "-graphene ring",
    "triangle": "-graphene triangle",
}

def geom_generation(data, yaml_file):
    """
    Generate geometry using GEOM from user-defined parameters.

    Parameters:
        data (dict): Dictionary of user-defined geometry parameters.
        yaml_file (str): Path to the input YAML file.
    """    

    errors = []

    # Check keys
    for key in data:
        if key not in VALID_KEYS:
            errors.append(f"Invalid key found: '{key}' is not a recognized parameter.")

    # Check types
    for key, value in data.items():
        expected_type = VALID_KEYS.get(key)
        if expected_type and not isinstance(value, expected_type):
            errors.append(f"Invalid type for '{key}': expected {expected_type.__name__}, got {type(value).__name__}.")

    # Conversions of the unit
    length_keys = [ #this can have units
        "radius",
        "length",
        "width",
        "core_radius",
        "core_length",
        "core_width",
        "shell_radius",
        "shell_length",
        "shell_width",
        "base_length",
        "side_length",
        "height",
        "mesh_size",
        "distance",
        "radius_in",
        "height_tip",
    ]
    for key in data:
        if key in length_keys and isinstance(data[key], str):  # Key = length
            value, unit, errors_parse = parse_value_with_unit(data[key]) # Default angstrom
            if errors_parse == "": 
                try:
                    data[key] = convert_length_to_angstrom(value, unit)
                except ValueError as e:
                    errors.append(str(e))  # Error conversion
            else:
                errors.append(errors_parse)  # Error parsing

    # Check `type` (if present)
    geom_type = data.get("type", "atomistic")  # Default is atomistic
    if geom_type not in ["atomistic", "implicit"]:
        errors.append(f"Invalid type '{geom_type}'. Must be 'atomistic' or 'implicit'.")

    # if `implicit`, shape: sphere or rod and must have mesh_size (default 10.0)
    if geom_type == "implicit":
        if "shape" not in data or data["shape"] not in ["sphere", "rod"]:
            errors.append("Implicit type requires 'shape' to be either 'sphere' or 'rod'.")
        if "mesh_size" not in data:
            data["mesh_size"] = 10.0  # Default mesh size
    else:
        # if `atomistic`, not  mesh_size
        if "mesh_size" in data:
            errors.append("'mesh_size' is not allowed for atomistic type.")

    # Check `shape`
    shape = data.get("shape")
    if not shape:
        errors.append("Missing required key: 'shape'.")
    elif shape not in SHAPE_REQUIRED_PARAMS:
        errors.append(f"Invalid shape '{shape}'. Must be one of {list(SHAPE_REQUIRED_PARAMS.keys())}.")
    
    # Check atomtype and core-shell
    has_atomtype = "atomtype" in data
    has_core = "core_atomtype" in data
    has_shell = "shell_atomtype" in data
    
    if has_atomtype:
        if has_core or has_shell:
            errors.append("Cannot specify 'atomtype' alongside 'core_atomtype' and 'shell_atomtype'.")
    elif has_core and has_shell:
        if shape == "sphere":
           if "core_radius" not in data or "shell_radius" not in data:
               errors.append("Core-shell spheres require both 'core_radius' and 'shell_radius'.")
        elif shape == "rod":
           if "core_length" not in data or "core_width" not in data or "shell_length" not in data or "shell_width" not in data:
               errors.append("Core-shell rods require 'core/shell_length' and 'core/shell_width'.")
    else:
        if geom_type == "atomistic" :
            errors.append("Either 'atomtype' or both 'core_atomtype' and 'shell_atomtype' must be specified.")
    
    # Check all requested parameters in `shape`
    if shape in SHAPE_REQUIRED_PARAMS:
        required_params = SHAPE_REQUIRED_PARAMS[shape]
        missing_keys = required_params - set(data.keys())
    
        # Exception sphere e rod : core-shell
        if shape == "sphere":
            if has_core and has_shell:
                missing_keys -= {"radius"}  # Core-shell uses core_radius and shell_radius instead of radius
            elif not has_core and not has_shell:
                missing_keys -= {"core_radius", "shell_radius"}  # Standard sphere
        elif shape == "rod":
            if has_core and has_shell:
                missing_keys -= {"length", "width"}  

        if missing_keys:
            errors.append(f"Missing required keys for shape '{shape}': {missing_keys}")

        # Check for unexpected parameters
        if shape in SHAPE_VALID_PARAMS:
            valid_keys = SHAPE_VALID_PARAMS[shape]  
            provided_keys = set(data.keys()) - {"shape"}  # Ignore 'shape' key
            extra_keys = provided_keys - valid_keys

            if extra_keys:
                errors.append(f"Unexpected keys for shape '{shape}': {extra_keys}")

    # Grafene
    if shape in {"disk", "ring", "ribbon", "triangle"} and data.get("atomtype") != "c":
        errors.append(f"Invalid 'atomtype' for {shape}. Must be 'C' (Graphene).")

    # Ring
    if shape == "ring":
        if "shell_radius" in data or "radius_in" in data:
            if "shell_radius" in data:
                data["radius_in"] = data["radius"] - data["shell_radius"]
            if data["radius_in"] <= 0:
                errors.append(f"Invalid 'shell_radius' ({data['shell_radius']}): must be smaller than 'radius' ({data['radius']}).")
        else:
            errors.append("Ring must specify either 'radius_in' or 'shell_radius'.")

    # Ribbon
    if shape == "ribbon":
        edge_type = data.get("edge_type", "armchair")  # Default: armchair
        if edge_type not in ["armchair", "zigzag", "ac", "zz"]:
            errors.append(f"Invalid 'edge_type' ({edge_type}). Must be 'armchair/ac' or 'zigzag/zz'.")

        if "width" in data and "length" in data:
            if edge_type in ["armchair", "ac"]:
                if data["width"] > data["length"]:
                    data["width"], data["length"] = data["length"], data["width"]  # Swap
            elif edge_type in ["zigzag", "zz"]:
                if data["length"] > data["width"]:
                    data["width"], data["length"] = data["length"], data["width"]  # Swap
        else:
            errors.append("Ribbon must specify both 'width' and 'length'.")

    # Errors
    if errors:
        return [], set(), 0, 0, "GEOM generation failed", errors

    # GEOM command
    if shape not in SHAPE_MAPPING:
        errors.append(f"Invalid shape '{shape}' provided!")

    geom_command = ["python3", "-m", "geom", "-create", SHAPE_MAPPING[shape]]


    # Standard parameters
    if "atomtype" in data and data["atomtype"] != "c":
        geom_command.append(data["atomtype"])
    if geom_type == "implicit":
        geom_command.append("-continuum")
    if shape in ["triangle"] and "edge_type" in data:
        geom_command.append(data["edge_type"])
    if shape in ["rod"] and "main_axis" in data:
        geom_command.append(data["main_axis"])
    if "height" in data:
        geom_command.append(str(data["height"]))
    if "radius" in data:
        geom_command.append(str(data["radius"]))
    if "length" in data:
        geom_command.append(str(data["length"]))
    if "width" in data:
        geom_command.append(str(data["width"]))
    if "base_length" in data:
        geom_command.append(str(data["base_length"]))
    if "side_length" in data:
        geom_command.append(str(data["side_length"]))
    if "radius_in" in data:
        geom_command.append(str(data["radius_in"]))

    # Core-shell parameters
    if has_core and has_shell:
        if shape == "sphere":
            geom_command.extend(["-core", data["core_atomtype"], str(data["core_radius"]),
                                 "-shell", data["shell_atomtype"], str(data["shell_radius"])])
        elif shape == "rod":
            if "core_length" in data and "core_width" in data and "shell_length" in data and "shell_width" in data:
                geom_command.extend(["-core", data["core_atomtype"], str(data["core_length"]), str(data["core_width"]),
                                     "-shell", data["shell_atomtype"], str(data["shell_length"]), str(data["shell_width"])])
            else:
                errors.append("Core-shell rods require 'core_length', 'core_width', 'shell_length', and 'shell_width'.")    

    if geom_type == "implicit":
        if "mesh_size" in data:
            geom_command.append(str(data["mesh_size"]))

    if shape in ["microscope", "paraboloid"]:
        if "a" in data and "b" in data:
            geom_command.extend([str(data["a"]), str(data["b"])])
        if "height" in data:
            geom_command.append(str(data["height"]))

    if "distance" in data: 
        dimer_type = data.get("dimer","default")
        if dimer_type not in ["default","bowtie"]:
            error.append(f"Invalid dimer type {dimer_type}. Only: default or bowtie.")

        if dimer_type == "bowtie":
            geom_command.extend(["-bowtie", str(data["distance"])])
        elif dimer_type == "default":
            geom_command.extend(["-dimer", str(data["distance"])])
            if "dimer_axis" in data:
                geom_command.append(data["dimer_axis"])
            

    # now call GEOM 
    geom_command_string = " ".join(geom_command)

    # Find GEOM and add to sys.path
    GEOM_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../submodules/geom"))
    if not GEOM_PATH in os.environ.get("PYTHONPATH",""):
        os.environ["PYTHONPATH"] = GEOM_PATH + ":" + os.environ.get("PYTHONPATH", "")

    try:
        result = subprocess.run(geom_command_string, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        error_message = e.stdout.strip() if e.stdout else "Unknown error from GEOM"
        errors.append(f"GEOM execution failed: {error_message}")
    
    # Check GEOM errors
    error_output = result.stdout.strip()
    if "ERROR" in error_output or "Error" in error_output or "error" in error_output:
        errors.append(f"GEOM execution failed: {error_output}")
        return [], set(), 0, 0, "GEOM execution failed", errors

    results_dir = "results_geom"
    
    # check GEOM folder
    if not os.path.exists(results_dir):
        errors.append("GEOM execution did not create the 'results_geom' directory.")
        return [], set(), 0, 0, "GEOM execution failed", errors
    
    if geom_type == "implicit":
        #msh
        msh_files = [f for f in os.listdir(results_dir) if f.endswith(".msh")]

        # only 1 .msh
        if len(msh_files) == 0:
            errors.append("No .msh file was found in 'results_geom'.")
            return [], set(), 0, 0, "GEOM execution failed", errors
        elif len(msh_files) > 1:
            errors.append("Multiple .msh files found in 'results_geom'. Cannot determine the correct one.")
            return [], set(), 0, 0, "GEOM execution failed", errors
        
        msh_file = msh_files[0]
        
        msh_new_name = os.path.splitext(yaml_file)[0] + ".msh"
        
        shutil.move(os.path.join(results_dir, msh_file), msh_new_name)
        shutil.rmtree(results_dir, ignore_errors=True)

        return [], set(), 0, 0, f"Generated msh by using GEOM: {geom_command_string}", errors

    elif geom_type == "atomistic": 
        # xyz in the folder
        xyz_files = [f for f in os.listdir(results_dir) if f.endswith(".xyz")]
        
        # Check 1 .xyz
        if len(xyz_files) == 0:
            errors.append("No .xyz file was found in 'results_geom'.")
            return [], set(), 0, 0, "GEOM execution failed", errors
        elif len(xyz_files) > 1:
            errors.append("Multiple .xyz files found in 'results_geom'. Cannot determine the correct one.")
            return [], set(), 0, 0, "GEOM execution failed", errors
        
        xyz_file = xyz_files[0]
        
        xyz_new_name = os.path.splitext(yaml_file)[0] + ".xyz"
        
        shutil.move(os.path.join(results_dir, xyz_file), xyz_new_name)
        shutil.rmtree(results_dir, ignore_errors=True)
        
        atoms = read(xyz_new_name)
        natoms = len(atoms)
        nmol = 1  # Default for single-molecule systems
        
        processed_geometry = []
        unique_atomtypes = set()
        
        for atom in atoms:
            processed_geometry.append((atom.symbol, 1, atom.position[0], atom.position[1], atom.position[2]))
            unique_atomtypes.add(atom.symbol)

        return processed_geometry, unique_atomtypes, natoms, nmol, f"Generated xyz by using GEOM: {geom_command_string}", errors
