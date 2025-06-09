#!/usr/bin/env python

# Global configurations
configurations = {
    "what": {
        "valid_keywords": [
            "energy",
            "static response",
            "dynamic response",
            "restart",
        ],
        "allowed_values": {},
        "type_spec": {}
    },
    "algorithm": {
        "valid_keywords": [
            "method",
            "parallel execution",
            "number of iterations",
            "gmres dimension",
            "tolerance",
            "adaptive tuning",
        ],
        "allowed_values": {
            "method": {"inversion", 
                       "iterative", 
                       "iterative on the fly"},
            "parallel execution": {"matrix", 
                                   "frequencies"},
            "adaptive tuning": {"yes", 
                                "no"}
        },
        "type_spec": {
            "method": str, 
            "parallel execution": str,
            "number of iterations": int,
            "gmres dimension": int,
            "tolerance": (float,str),
            "adaptive tuning": str
        }
    },
    "forcefield": {
        "valid_keywords": [
            "static",
            "dynamic",
            "kernel",
        ],
        "allowed_values": {
            "static": {"fq", 
                       "fqfmu"},
            "dynamic": {"wfq", 
                        "wfqfmu"},
            "kernel": {"ohno", 
                       "gaussian",
                       "coulomb"}
         },
        "type_spec": {
            "static": str, 
            "dynamic": str,
            "kernel": str,
        }
    },
    "field": {
        "valid_keywords": [
            "type", 
            "rhs type", 
            "field intensity", 
            "nfreq",
            "min freq", 
            "max freq", 
            "step freq", 
            "external freq", 
            "polarization"
        ],
        "type_spec": {
            "type": str,
            "rhs type": str,
            "field intensity": float,
            "nfreq": int,
            "min freq": (float, str),
            "max freq": (float, str),
            "step freq": (float, str),
            "external freq": (float, str),  
            "polarization": str
        },
        "allowed_values": {
            "type": {"static", "dynamic"},
            "rhs type": {"field", "potential"},
            "polarization": {"none", "x", "y", "z", "all"}
        },
    },
    "control": {
        "valid_keywords": [
            "no info file", 
            "principal axes"
        ],
        "allowed_values": {},
        "type_spec": {}
    },
    "output": {
        "valid_keywords": [
            "verbose", 
            "maxima analysis"
        ],
        "allowed_values": {
            "maxima analysis": {"absorption",
                                "scattering",
                                "excitinction"}
        },
        "type_spec": {
            "verbose": int,
            "maxima analysis": str
        }
    },
    "bem": {
        "valid_keywords": [
            "mesh file",
            "normal scalar factor",
            "permittivity file",
            "permittivity",
            "green function",
            "sphere radius",
            "solvent",
            "solvent epsilon",
            "variant"
        ],
        "allowed_values": {
            "normal scalar factor": {1.0,
                                     -1.0},
            "permittivity" : {"silver etchegoin", 
                              "silver jc",
                              "silver johnson-christy",
                              "silver johnson christy",
                              "silver bb",
                              "silver brendel-bormann",
                              "silver brendel bormann",
                              "silver palik",
                              "gold etchegoin",
                              "gold jc",
                              "gold johnson-christy",
                              "gold johnson christy",
                              "gold bb",
                              "gold brendel-bormann",
                              "gold brendel bormann",
                              "gold palik"},
            "green function": {"exact",
                               "approximate"},
            "variant" : {"dpcm",
                         "iefpcm"}
        },
        "type_spec": {
            "mesh file": str,
            "normal scalar factor": float,
            "permittivity file": str,
            "permittivity": str,
            "green function": str,
            "sphere radius": (float,str),
            "solvent": str,
            "solvent epsilon": float,
            "variant": str,
        }
    },
    "input_geometry": {
        "valid_keywords": [
            "external xyz file", 
            "type", 
            "shape", 
            "atomtype",
            "radius",
            "main_axis", 
            "length", "width",
            "core_atomtype", 
            "core_radius",
            "core_length",
            "core_width",
            "shell_atomtype",
            "shell_radius",
            "shell_length",
            "shell_width",
            "base_length",
            "side_length", 
            "height",
            "mesh_size", 
            "distance", 
            "edge_type",
            "a", 
            "b", 
            "dimer_axis",
            "radius_in", 
            "height_tip",
            "dimer"
        ],
        "allowed_values": {},
        "type_spec": {
            "external xyz file": str,
            "type": str,
            "shape": str,
            "atomtype": str,
            "radius": (float,str),
            "main_axis": str,
            "length": (float,str),
            "width": (float,str),
            "core_atomtype": str,
            "core_radius": (float,str),
            "core_length": (float,str),
            "core_width": (float,str),
            "shell_atomtype": str,
            "shell_radius": (float,str),
            "shell_length": (float,str),
            "shell_width": (float,str),
            "base_length": (float,str),
            "side_length": (float,str),
            "height": (float,str),
            "mesh_size": (float,str),
            "distance": (float,str),
            "edge_type": str,
            "a": float,
            "b": float,
            "dimer_axis": str,
            "radius_in": (float,str),
            "height_tip": (float,str),
            "dimer": str
        }
    },
    # Placeholder for dynamic sections
    "atom_types": {
        "valid_keywords": [
            "chi", 
            "eta", 
            "alpha", 
            "rq", 
            "rmu"
        ],
        "allowed_values": {},
        "type_spec": {
            "chi": (float,str),
            "eta": (float,str),
            "alpha": (float,str),
            "rq": (float,str),
            "rmu": (float,str)
        }
    },
    "parameters": {
        "valid_keywords": [
            "atomtype tau", 
            "atomtype sigma0", 
            "atomtype w_p", 
            "atomtype gamma",
            "atomtype scaling sigma0-tau", 
            "atomtype a_ij", 
            "atomtype fermi function d",
            "atomtype fermi function s", 
            "atomtype fermi energy", 
            "atomtype wfqfmu file",
            "atomtype permittivity", 
            "interaction fermi function d", 
            "interaction fermi function s"
        ],
        "allowed_values": { 
            "atomtype permittivity" : {"silver etchegoin", 
                                       "gold etchegoin"},
        },
        "type_spec": {
            "atomtype tau": (float,str),
            "atomtype sigma0": (float,str),
            "atomtype w_p": (float,str),
            "atomtype gamma": (float,str),
            "atomtype scaling sigma0-tau": (float),
            "atomtype a_ij": (float,str),
            "atomtype fermi function d": (float,str),
            "atomtype fermi function s": (float,str),
            "atomtype fermi energy": (float,str),
            "atomtype wfqfmu file": str,
            "atomtype permittivity": str,
            "interaction fermi function d": (float,str),
            "interaction fermi function s": (float,str),
        }
    },
}

def create_keywords(data, input_atomtypes):
    """
    Initialize configuration entries for atom types and parameters based on input atom types.

    Parameters:
        data (dict): Unused input data.
        input_atomtypes (list): List of atom type strings.
    """    

    # Configuration initial
    base_atomtype_config = configurations["atom_types"]

    # dictionary for atom_types
    expanded_atomtypes = {}
    for atom in input_atomtypes:
        atom_lower = atom.lower()
        expanded_atomtypes[atom_lower] = {
            "valid_keywords": base_atomtype_config["valid_keywords"].copy(),
            "allowed_values": base_atomtype_config["allowed_values"].copy(),
            "type_spec": base_atomtype_config["type_spec"].copy()
        }

    configurations["atom_types"] = expanded_atomtypes    

    # Configuration to be copied
    base_parameters_config = configurations["parameters"]

    # Nuovo dizionario da sostituire nella sezione parameters
    expanded_parameters = {}
    for atom in input_atomtypes:
        atom_lower = atom.lower()
        expanded_parameters[atom_lower] = {
            "valid_keywords": [
                k.replace("atomtype ", "")
                for k in base_parameters_config["valid_keywords"]
                if k.startswith("atomtype ")
            ],
            "allowed_values": {
                k.replace("atomtype ", ""): v
                for k, v in base_parameters_config["allowed_values"].items()
                if k.startswith("atomtype ")
            },
            "type_spec": {
                k.replace("atomtype ", ""): v
                for k, v in base_parameters_config["type_spec"].items()
                if k.startswith("atomtype ")
            }
        }

        # Interactions update parameters allowed
        for a in input_atomtypes:
            for b in input_atomtypes:
                if a != b:
                    a_lower = a.lower()
                    b_lower = b.lower()
                    interaction_key = f"interaction {a_lower}->{b_lower}"
                    expanded_parameters[interaction_key] = {
                        "valid_keywords": [
                            k.replace("interaction ", "")
                            for k in base_parameters_config["valid_keywords"]
                            if k.startswith("interaction ")
                        ],
                        "allowed_values": {
                            k.replace("interaction ", ""): v
                            for k, v in base_parameters_config["allowed_values"].items()
                            if k.startswith("interaction ")
                        },
                        "type_spec": {
                            k.replace("interaction ", ""): v
                            for k, v in base_parameters_config["type_spec"].items()
                            if k.startswith("interaction ")
                        }
                    }

    configurations["parameters"] = expanded_parameters

    return configurations

def create_starting_keywords():
    """
    Return the initial configuration dictionary.
    """
    
    conf_starting = configurations
    return conf_starting
    
def check_unknown_section(data, configurations):
    """
    Check if all sections in the input are recognized in the configuration.

    Parameters:
        data (dict): Input data with sections to validate.
        configurations (dict): Reference configuration with allowed sections.
    """    
    errors = []
    allowed_sections = list(configurations.keys())

    formatted_sections = "\n    - " + "\n    - ".join(allowed_sections)

    for section in data:
        if section not in configurations:
            errors.append(
                f"Unknown section '{section}' in input.\nAllowed sections are:{formatted_sections}"
            )

    return errors    
    
def check_section(section_name, section_data, configurations):
    """
    Validate a configuration section and its entries against allowed structure and known atom types.

    Parameters:
        section_name (str): Name of the section to validate.
        section_data (dict): Data entries for the section.
        configurations (dict): Full configuration reference.
    """
    # configuration of the section
    section_config = configurations.get(section_name)

    errors = []

    if not section_config:
        errors.append(f"Section '{section_name}' not allowed.")
        return errors

    if section_name == "atom_types" or section_name == "parameters":
        for atomtype, atom_data in section_data.items():
            # Handle "interaction ..." cases separately
            if section_name == "parameters" and atomtype.startswith("interaction "):
                try:
                    interaction = atomtype.split(" ", 1)[1]
                    a, b = interaction.split("->")
                    a = a.strip()
                    b = b.strip()
    
                    # Check that both atomtypes are in the known list
                    if a not in section_config and b not in section_config:
                        errors.append(
                            f"Invalid interaction '{interaction}' in section '{section_name}': atomtypes '{a}' and '{b}' not found in geometry."
                        )
                    elif a not in section_config:
                        errors.append(
                            f"Invalid interaction '{interaction}' in section '{section_name}': atomtype '{a}' not found in geometry."
                        )
                    elif b not in section_config:
                        errors.append(
                            f"Invalid interaction '{interaction}' in section '{section_name}': atomtype '{b}' not found in geometry."
                        )
                except Exception:
                    errors.append(
                        f"Wrong interaction name '{atomtype}' in section '{section_name}'. Expected format: 'interaction A->B'"
                    )
            # all the rests
            atom_config = section_config.get(atomtype.lower())
            if not atom_config:
                #Try to split to see if there is an errore in interaction
                parts = atomtype.strip().split()
                if len(parts) == 2:
                    # This should be an interaction
                    if parts[0] != "interaction":
                        errors.append(
                            f"Invalid interaction key '{parts[0]}'. Expected: 'interaction'"
                        )
                elif len(parts) == 1:
                    probable_atom = parts[0].strip()
                    # Heuristic: if it's short, likely an atomtype
                    if len(probable_atom) <= 5:
                        errors.append(
                            f"Atom type '{probable_atom}' in section '{section_name}' is not defined in input geometry"
                        )
                    else:
                        errors.append(
                            f"Invalid interaction or atomtype key '{atomtype}'. Unexpected format"
                        )
                else:
                    errors.append(
                        f"Invalid key '{atomtype}' in section '{section_name}'. Unexpected format"
                    )                
                continue
            section_name_atom = "atom_types: "+atomtype
            errors.extend(check_keywords(section_name_atom, atom_data, atom_config["valid_keywords"]))
            errors.extend(check_allowed_values(section_name_atom, atom_data, atom_config["allowed_values"]))
            errors.extend(check_data_types(section_name_atom, atom_data, atom_config["type_spec"]))
    else: 
        errors.extend(check_keywords(section_name, section_data, section_config["valid_keywords"]))
        errors.extend(check_allowed_values(section_name, section_data, section_config["allowed_values"]))
        errors.extend(check_data_types(section_name, section_data, section_config["type_spec"]))

    return errors


def check_keywords(section_name, section_data, valid_keywords):
    """
    Validate that all keywords in a section match the list of allowed keywords.

    Parameters:
        section_name (str): Name of the section being checked.
        section_data (str | list | dict): Data associated with the section.
        valid_keywords (list): List of allowed keywords for validation.
    """
    errors = []

    error_yaml = False
    if isinstance(section_data, list):
        for item in section_data:
            if item not in valid_keywords:
                errors.append(f"Invalid keyword '{item}' in section '{section_name}'.")
    
    elif isinstance(section_data, dict):
        for keyword, value in section_data.items():
            if isinstance(value, dict):
                # Keyword annidata, es. fermi_function: {d, s}
                for subkey in value:
                    full_key = f"{keyword.replace('_', ' ')} {subkey}"
                    if full_key not in valid_keywords:
                        errors.append(f"Invalid keyword '{full_key}' in section '{section_name}'.")
            else:
                if keyword not in valid_keywords:
                    errors.append(f"Invalid keyword '{keyword}' in section '{section_name}'.")
    elif isinstance(section_data, str):
        # Only for input geometry
        if section_name == "input_geometry":
            # Se contiene almeno un \n la trattiamo come input valido multilinea
            if "\n" in section_data.strip():
                pass  # ok, no error
            else:
                errors.append(f"Expected multiline string for '{section_name}'.\nYou probably forgot '|' after 'input_geometry:' in YAML file.")
                error_yaml = True
        else:
            errors.append(f"Unsupported string format in section '{section_name}'.")        
    else:
        errors.append(f"Unsupported type in section '{section_name}'.")

    if errors != [] and not error_yaml:
        formatted_keywords = "\n".join(f"    - {kw}" for kw in valid_keywords)
        errors.append(f"Allowed keywords for section '{section_name}' are:\n{formatted_keywords}")

    return errors

def check_allowed_values(section_name, section_data, allowed_values):
    """
    Validate that keyword values in a section are among the allowed options.

    Parameters:
        section_name (str): Name of the section being validated.
        section_data (dict): Dictionary of keyword-value pairs.
        allowed_values (dict): Dictionary mapping keywords to their allowed values.
    """
    errors = []
    
    if isinstance(section_data, dict):
        for keyword, value in section_data.items():
            if keyword in allowed_values:
                allowed = allowed_values[keyword]
                if value not in allowed:
                    formatted_allowed = "\n".join(f"    - {val}" for val in allowed)
                    errors.append(f"Invalid value '{value}' for keyword '{keyword}' in section '{section_name}'. Allowed values are:\n{formatted_allowed}.")
    return errors


def check_data_types(section_name, section_data, type_spec):
    """
    Check that each keyword in a section has a value of the expected data type.

    Parameters:
        section_name (str): Name of the section being validated.
        section_data (dict): Dictionary of keyword-value pairs.
        type_spec (dict): Dictionary mapping keywords to their expected data types.
    """
    errors = []

    if isinstance(section_data, dict):
        for keyword, value in section_data.items():
            expected_type = type_spec.get(keyword)
            if expected_type:
                if isinstance(expected_type, tuple):  # Supporto per tipi multipli (e.g., float or string)
                    if not isinstance(value, expected_type):
                        errors.append(
                            f"Invalid type for '{keyword}' in section '{section_name}'.\n"
                            f"  Expected: {expected_type}\n"
                            f"  Got:      {type(value)}"
                        )
                else:
                    if not isinstance(value, expected_type):
                        errors.append(
                            f"Invalid type for '{keyword}' in section '{section_name}'.\n"
                            f"  Expected: {expected_type}\n"
                            f"  Got:      {type(value)}"
                        )

    return errors
