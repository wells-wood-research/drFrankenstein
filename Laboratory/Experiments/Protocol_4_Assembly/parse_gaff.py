import yaml
import io
import re

def parse_amber_force_field(file_content, constant_charge = None):
    """
    Parses the AMBER General Force Field file content into a nested dictionary.

    Args:
        file_content (str): The string content of the force field file.

    Returns:
        dict: A dictionary with the parsed force field data in the format
              {ELEMENT: {ATOM_TYPE: {CHARGE, NONBONDED_PARAMETERS, MASS, COMMENT}}}.
    """
    force_field_dict = {}
    non_bonded_radii = {}
    non_bonded_well_depths = {}
    two_letter_elements = {
        'Ac', 'ac', 'Ag', 'ag', 'Al', 'al', 'Am', 'am', 'Ar', 'ar', 'As', 'as', 'At', 'at', 'Au', 'au',
        'Ba', 'ba', 'Be', 'be', 'Bh', 'bh', 'Bi', 'bi', 'Bk', 'bk', 'Br', 'br', 'Ca', 'ca', 'Cd', 'cd',
        'Ce', 'ce', 'Cf', 'cf', 'Cl', 'cl', 'Cm', 'cm', 'Cn', 'cn', 'Co', 'co', 'Cr', 'cr', 'Cs', 'cs',
        'Cu', 'cu', 'Db', 'db', 'Ds', 'ds', 'Dy', 'dy', 'Er', 'er', 'Es', 'es', 'Eu', 'eu', 'Fe', 'fe',
        'Fl', 'fl', 'Fm', 'fm', 'Fr', 'fr', 'Ga', 'ga', 'Gd', 'gd', 'Ge', 'ge', 'He', 'he', 'Hf', 'hf',
        'Hg', 'hg', 'Ho', 'ho', 'Hs', 'hs', 'In', 'in', 'Ir', 'ir', 'Kr', 'kr', 'La', 'la', 'Li', 'li',
        'Lr', 'lr', 'Lu', 'lu', 'Lv', 'lv', 'Mc', 'mc', 'Md', 'md', 'Mg', 'mg', 'Mn', 'mn', 'Mo', 'mo',
        'Mt', 'mt', 'Na', 'na', 'Nb', 'nb', 'Nd', 'nd', 'Ne', 'ne', 'Ni', 'ni', 'No', 'no', 'Np', 'np',
        'Og', 'og', 'Os', 'os', 'Pa', 'pa', 'Pb', 'pb', 'Pd', 'pd', 'Pm', 'pm', 'Po', 'po', 'Pr', 'pr',
        'Pt', 'pt', 'Pu', 'pu', 'Ra', 'ra', 'Rb', 'rb', 'Re', 're', 'Rf', 'rf', 'Rg', 'rg', 'Rh', 'rh',
        'Rn', 'rn', 'Ru', 'ru', 'Sb', 'sb', 'Sc', 'sc', 'Se', 'se', 'Sg', 'sg', 'Si', 'si', 'Sm', 'sm',
        'Sn', 'sn', 'Sr', 'sr', 'Ta', 'ta', 'Tb', 'tb', 'Tc', 'tc', 'Te', 'te', 'Th', 'th', 'Ti', 'ti',
        'Tl', 'tl', 'Tm', 'tm', 'Ts', 'ts', 'Yb', 'yb', 'Zn', 'zn', 'Zr', 'zr'
    }
    file = io.StringIO(file_content)
    lines = file.readlines()

    parsing_atom_types = False
    parsing_mod4_re = False
    
    # This regex helps identify where the atom type definitions end and other sections begin
    section_break_marker = re.compile(r'^[a-z0-9\+]{1,2}-[a-z0-9\+]{1,2}\s+')

    # --- Part 1: Parse Atom Types, Mass, and Comments ---
    for line in lines:
        line = line.strip()
        print(line)
        if not line:
            parsing_atom_types = False
            continue
        if "MOD4      RE" in line or "NONBON" in line:
            parsing_mod4_re = True
            parsing_atom_types = False
            continue
        elif "MASS" in line or line.startswith("AMBER "):
            parsing_atom_types = True
            continue
    
        
        if section_break_marker.match(line):
            parsing_atom_types = False

        if parsing_atom_types:
            parts = line.split(None, 3)
            if parts[0][0].isalpha():
                atom_type = parts[0]
                mass = float(parts[1])
                if not constant_charge:
                    charge = float(parts[2])
                else:
                    chargeSign = re.search(r'[+-]', atom_type)[0]
                    chargeValue = re.search(r'[0-9]', atom_type)
                   
                    if chargeValue:
                        chargeValue = chargeValue[0]
                    else:
                        chargeValue = 1

                    if chargeSign == '-':
                        charge = float(chargeValue) * -1
                    else:
                        charge = float(chargeValue)
                comment = ' '.join(parts[3:])

                # Determine the element from the atom type
                element = ''
                if atom_type[:2].lower() in two_letter_elements:
                    element = atom_type[:2].capitalize()
                else:
                    element = atom_type[0].upper()

                if element not in force_field_dict:
                    force_field_dict[element] = {}

                force_field_dict[element][atom_type] = {
                    'MASS': mass,
                    'COMMENT': comment.strip(),
                    'CHARGE': charge,
                    'RADIUS': None,
                    'WELL_DEPTH': None
                }
        
        if parsing_mod4_re:
            if "END" in line:
                parsing_mod4_re = False
                continue
            parts = line.split()
            if parts[0][0].isalpha():
                try:
                    atom_type = parts[0]
                    radius = float(parts[1])
                    well_depth = float(parts[2])
                    non_bonded_radii[atom_type] = radius
                    non_bonded_well_depths[atom_type] = well_depth
                except (ValueError, IndexError):
                    continue

    # --- Part 2: Merge Nonbonded Parameters ---
    for element in force_field_dict:
        for atom_type in force_field_dict[element]:
            if atom_type in non_bonded_radii:
                force_field_dict[element][atom_type]['RADIUS'] = non_bonded_radii[atom_type]
            if atom_type in non_bonded_well_depths:
                force_field_dict[element][atom_type]['WELL_DEPTH'] = non_bonded_well_depths[atom_type]

    return force_field_dict

def write_dict_to_yaml(data_dict, file_path):
    """
    Writes a dictionary to a YAML file.

    Args:
        data_dict (dict): The dictionary to write.
        file_path (str): The path to the output YAML file.
    """
    with open(file_path, 'w') as yaml_file:
        yaml.dump(data_dict, yaml_file, sort_keys=False, default_flow_style=False, indent=2)

def merge_and_update_dicts(dicts):
    """
    Merges multiple dictionaries, and for shared keys, it merges their
    sub-dictionaries.
    """
    megaDict = {}
    for paramDict in dicts:
        for key, value in paramDict.items():
            if not key in megaDict:
                megaDict[key] = value
            else:
                megaDict[key].update(value)
    return megaDict












if __name__ == '__main__':
    # Read the AMBER General Force Field file content
    with open('/home/eugene/.conda/envs/Igor/dat/leap/parm/gaff2.dat', 'r') as f:
        file_content = f.read()

    # Parse the file content into a dictionary
    gaffDict = parse_amber_force_field(file_content, constant_charge = False)


    monovalentIons = "/home/eugene/.conda/envs/Igor/dat/leap/parm/frcmod.ions1lm_126_tip3p"
    with open(monovalentIons, 'r') as f:
        monovalentIons = f.read()

    monovalentIonDict = parse_amber_force_field(monovalentIons, constant_charge = True)


    divalentIons = "/home/eugene/.conda/envs/Igor/dat/leap/parm/frcmod.ions234lm_1264_tip3p"
    with open(divalentIons, 'r') as f:
        divalentIons = f.read()

    divalentIonDict = parse_amber_force_field(divalentIons, constant_charge = True)


    force_field_dict = merge_and_update_dicts([gaffDict, monovalentIonDict, divalentIonDict])
    

    write_dict_to_yaml(force_field_dict, "/home/eugene/drFrankenstein/Laboratory/Experiments/Protocol_4_Assembly/gaff2.yaml")