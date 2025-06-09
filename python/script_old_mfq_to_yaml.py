#!/usr/bin/env python

import yaml
import argparse


def parse_input_file(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Creiamo un iteratore per scorrere le righe
    lines_iterator = iter(lines)

    #rename map
    key_rename_map = {
        "efield": "type",
        "iter": "method",
        "otf": "method",
        "onthefly": "method",
        "iterative": "method",
        "inversion": "method",
        "matrix in parallel": "parallel execution",
        "freq in parallel": "parallel execution",
        "no change parameters": "adaptive tuning",
        "n_iter": "number of iterations",
        "n_diis": "gmres dimension",
        "rhs": "rhs type",
        "step": "step freq",
        "gmsh file": "mesh file",
        "scaling": "scaling sigma0-tau",
        "ri": "A_ij"
    }
    value_rename_map = {
        "iter": "iterative",
        "otf": "iterative on the fly",
        "onthefly": "iterative on the fly",
        "matrix in parallel": "matrix",  # Nuova mappatura valore
        "freq in parallel": "frequencies",  # Nuova mappatura valore
        "no change parameters": "no",  # Nuova mappatura valore
        "silver": "silver etchegoin", 
        "gold": "gold etchegoin",
        "fqfmu_pqeq": "fqfmu"
    }

    data = {}
    current_section = None
    current_subsection = None

    for line in lines_iterator:
        line = line.strip()
        if not line or line.startswith("end"):
            current_section = None
            continue


        # Inizializza le sezioni principali
        if line in ["what", "algorithm", "forcefield", "field", "parameters", "atom types", "input geometry", "control", "output", "bem"]:
            current_section = line.replace(" ", "_")
            if current_section in ["what", "input_geometry", "control"]:
                data[current_section] = []  # Queste sezioni devono essere liste
            else:
                data[current_section] = {}  # Le altre sezioni sono dizionari
            continue

        # Aggiungi i dati nelle sezioni corrispondenti
        if current_section == "what":
            if (line == "cross section"):
                data[current_section].append("dynamic response")
            else:
                data[current_section].append(line)

        elif current_section == "forcefield":
            if ":" in line:
                key, value = map(str.strip, line.split(":", 1))  # Divide solo alla prima occorrenza di ":"

                # Applica la mappatura delle chiavi
                if key in key_rename_map:
                    key = key_rename_map[key]
                # Rinomina il valore se presente nel dizionario
                if value in value_rename_map:
                    value = value_rename_map[value]

                data[current_section][key] = float(value) if value.replace(".", "", 1).isdigit() else value
            else:
                raise ValueError(f"Formato non valido nella sezione {current_section}: {line}")

        elif current_section == "field":

            if ":" in line:
                key, value = map(str.strip, line.split(":", 1))  # Divide solo alla prima occorrenza di ":"

                # Applica la mappatura delle chiavi
                if key in key_rename_map:
                    key = key_rename_map[key]

                if(key != "nfreq") :
                    data[current_section][key] = float(value) if value.replace(".", "", 1).isdigit() else value
                else:
                    data[current_section][key] = int(value) if value.replace(".", "", 1).isdigit() else value
            else:
                raise ValueError(f"Formato non valido nella sezione {current_section}: {line}")

            # **Correzione: Aggiorna min freq se necessario**
            if "min freq" in data[current_section] and "step freq" in data[current_section]:
                data[current_section]["min freq"] += data[current_section]["step freq"]
       
        elif current_section == "output":
            if ":" in line:
                key, value = map(str.strip, line.split(":", 1))  # Divide solo alla prima occorrenza di ":"

                # Applica la mappatura delle chiavi
                if key in key_rename_map:
                    key = key_rename_map[key]
                data[current_section][key] = int(value) if value.replace(".", "", 1).isdigit() else value
            else:
                raise ValueError(f"Formato non valido nella sezione {current_section}: {line}")

        #algorithm
        elif current_section == "algorithm":
            if ":" in line:
                key, value = map(str.strip, line.split(":", 1))

                # Applica la mappatura delle chiavi
                if key in key_rename_map:
                    key = key_rename_map[key]
                # Rinomina il valore se presente nel dizionario
                if value in value_rename_map:
                    value = value_rename_map[value]
                # Conversione condizionale in float o mantieni stringa
                try:
                    value = float(value) if "E" in value.upper() or "." in value else float(value.replace("d", "e"))
                except ValueError:
                    pass
                data[current_section][key] = value
            else:
                # Usa la mappatura per aggiornare la chiave
                if line in key_rename_map:
                    mapped_key = key_rename_map[line]
                    normalized_value = value_rename_map.get(line, line.lower())
                    data[current_section][mapped_key] = normalized_value
                else:
                    data[current_section][line] = True

            # Set adaptive tuning to "no" unless "no change parameters" is explicitly provided
            if "adaptive tuning" not in data[current_section]:
                data[current_section]["adaptive tuning"] = "no"

        # control
        elif current_section == "control":
            data[current_section] = []
            while True:
                if line and not line.startswith("end"):
                    data[current_section].append(line)  # Aggiungi la riga corrente prima di leggere la prossima
                line = next(lines_iterator, "").strip()  # Leggi la prossima riga
                if not line or line.startswith("end"):
                    break        

        #bem
        elif current_section == "bem":
            if ":" in line:
                key, value = map(str.strip, line.split(":", 1))  # Divide alla prima occorrenza di ":"
        
                # Applica una mappatura delle chiavi, se necessario
                if key in key_rename_map:
                    key = key_rename_map[key]
                if value in value_rename_map:
                    value = value_rename_map[value]
        
                # Conversione condizionale in tipi appropriati
                try:
                    if "d" in value.lower():
                        value = float(value.replace("d", "e"))  # Gestione dei numeri in notazione "d"
                    elif value.replace(".", "", 1).isdigit():
                        value = float(value) if "." in value else int(value)
                except ValueError:
                    pass  # Mantiene il valore come stringa se non può essere convertito
        
                data[current_section][key] = value
            else:
                raise ValueError(f"Formato non valido nella sezione {current_section}: {line}")

        elif current_section == "parameters":
            parameters_lines = []
            parameters_lines.append(line)
        
            # Raccogli tutte le righe fino a "end parameters"
            while True:
                line = next(lines_iterator, "").strip()
                if not line or line == "end parameters":
                    break
                parameters_lines.append(line)
        
            # Analizza ogni riga raccolta
            current_subsection = None
            for param_line in parameters_lines:
                if param_line.startswith("atomtype") or param_line.startswith("interaction"):
                    current_subsection = param_line.strip()  # Usa l'intera riga come chiave
                    data[current_section][current_subsection] = {}
                elif param_line.startswith("fermi function"):
                    # Gestione delle "fermi function"
                    parts = param_line.split(":", 1)
                    if len(parts) == 2:
                        key = parts[0].strip().replace("fermi function ", "")
                        value = float(parts[1].strip())
                        data[current_section][current_subsection].setdefault("fermi_function", {})[key] = value
                    else:
                        raise ValueError(f"Formato non valido nella sezione parameters: {param_line}")
                elif ":" in param_line:
                    # Aggiunta delle altre chiavi-valore
                    key, value = map(str.strip, param_line.split(":", 1))
        
                    # Applica una mappatura delle chiavi, se necessario
                    if key.lower() in key_rename_map:
                        key = key_rename_map[key.lower()]
                    # Se il parametro è "ri", rinominiamolo in "A_ij" e calcoliamo A_ij = valore²
                    if key == "A_ij":
                        try:
                            value = float(value) ** 2  # Sovrascriviamo con A_ij = ri²
                        except ValueError:
                            errors.append(f"Invalid value for 'ri'. Must be a numerical value.")
                    
                    if value in value_rename_map:
                        value = value_rename_map[value]
                    if isinstance(value, str) and ("." in value or "E" in value):
                        try:
                            value = float(value)
                        except ValueError:
                            pass  # Se non è convertibile in float, lo lasciamo come stringa                    
                    data[current_section][current_subsection][key] = value

            # Controlla se "scaling sigma0-tau" è presente solo per gli atomtype
            for subsection, properties in data[current_section].items():
                if subsection.startswith("atomtype") and "scaling sigma0-tau" not in properties:
                    properties["scaling sigma0-tau"] = 1.0
        
        elif current_section == "atom_types":
            if "number" in line:
                # Gestisce la riga con "number"
                key, value = map(str.strip, line.split(":", 1))
                # Applica la mappatura delle chiavi
                if key in key_rename_map:
                    key = key_rename_map[key]
                data[current_section][key] = int(value)
            elif ":" in line:
                # Gestisce le righe con il formato Na: [chi=0.000000, eta=0.292000]
                atomtype, properties = line.split(":", 1)
                properties = properties.replace("[", "").replace("]", "").replace(" ", "").split(",")

                # **Inizializziamo parsed_properties per evitare problemi**
                parsed_properties = {}
        
                for prop in properties:
                    key, val = prop.split("=")
                    parsed_properties[key.strip()] = float(val.strip())
        
                # **Se rq e/o rmu sono presenti, aggiorniamo con la formula sqrt(x^2 / 0.2341)**
                if "rq" in parsed_properties:
                    parsed_properties["rq"] = (parsed_properties["rq"] ** 2 / 0.2341) ** 0.5
                if "rmu" in parsed_properties:
                    parsed_properties["rmu"] = (parsed_properties["rmu"] ** 2 / 0.2341) ** 0.5
        
                data[current_section][atomtype.strip()] = parsed_properties

        elif current_section == "input_geometry":
            if current_section not in data:
                data[current_section] = []
            # Inizializza come dizionario o lista in base al contenuto
            if line.startswith("external xyz file:"):
                if not isinstance(data.get(current_section, {}), dict):
                    data[current_section] = {}  # Converte in dizionario
                key, value = map(str.strip, line.split(":", 1))
                data[current_section][key] = value
            else:
                if not isinstance(data.get(current_section, []), list):
                    data[current_section] = []  # Converte in lista
                data[current_section].append(line)

    return data

# Funzione personalizzata per rappresentare stringhe multilinea come blocco letterale
def represent_multiline(dumper, data):
    if isinstance(data, str) and "\n" in data:
        return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")
    return dumper.represent_scalar("tag:yaml.org,2002:str", data)

class CustomDumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        # Garantisce rientri corretti senza spazi superflui
        super().increase_indent(flow, indentless=False)

def write_yaml(data, output_path):

    yaml.add_representer(str, represent_multiline)
    if "input_geometry" in data:
        if isinstance(data["input_geometry"], list):
            # Convertire in blocco letterale se è una lista
            data["input_geometry"] = "\n".join(data["input_geometry"]) + "\n"
        elif isinstance(data["input_geometry"], dict):
            # Lascia come dizionario se contiene chiavi come "external xyz file"
            pass


    # Scrivere il file YAML con righe vuote tra sezioni
    with open(output_path, "w") as file:
        yaml.dump(
            data,
            file,
            sort_keys=False,
            Dumper=CustomDumper,
            default_flow_style=False,  # Usa blocchi leggibili
            allow_unicode=True
        )    
    # Aggiungere righe vuote tra le sezioni
    with open(output_path, "r") as file:
        lines = file.readlines()

    with open(output_path, "w") as file:
        for i, line in enumerate(lines):
            # Se la riga contiene una sezione principale, aggiunge una riga vuota PRIMA di scriverla
            if line.strip() in ["forcefield:", "algorithm:", "field:", "output:", "parameters:", "atom_types:", "control:", "bem:", "input_geometry: |", "input_geometry:"]:
                file.write("\n")
            # Scrive la riga corrente
            file.write(line)
    

def main():
    parser = argparse.ArgumentParser(description="Convert an old .mfq input file to YAML format. By default a file with the same name and .yaml extension is created")
    parser.add_argument("-i", "--input_file", type=str, required=True,
                        help="input file [old .mfq file]")
    parser.add_argument("-o", "--output_file", type=str, required=False,
                        help="Optional: output file [YAML]")
    args = parser.parse_args()

    if args.output_file:
        # Scrive il file YAML se è specificato l'output_file
        output_ = args.output_file
    else: 
        output_ = args.input_file[:-4]+".yaml"

    # Parse and convert
    parsed_data = parse_input_file(args.input_file)
    write_yaml(parsed_data, output_)

    print(f"File YAML scritto in: {output_}")

if __name__ == "__main__":
    main()

