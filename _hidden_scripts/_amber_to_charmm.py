import os
from os import path as p
import pandas as pd
import yaml

#####################################################
def inputs():
    creatorDir = "/home/esp/scriptDevelopment/drFrankenstein/Benchmarking/03_GLY_XTB_MP2_XTB_MP2/06_final_creation"
    moleculeName = "GLY"

    configFile = "/home/esp/scriptDevelopment/drFrankenstein/Benchmarking/03_GLY_XTB_MP2_XTB_MP2/drFrankenstein.yaml"
    with open(configFile,"r") as f:
        config = yaml.safe_load(f)

    amberFrcmod = p.join(creatorDir, f"{moleculeName}.frcmod")
    amberLib = p.join(creatorDir, f"{moleculeName}.lib")
    amberMol2 = p.join(creatorDir, f"{moleculeName}.mol2")

    return amberFrcmod, amberLib, amberMol2, creatorDir, moleculeName, config
#####################################################

def main():
    amberFrcmod, amberLib, amberMol2, creatorDir, moleculeName, config = inputs()

    charmmParamDir = p.join(creatorDir, "charmm_paramerters")
    os.makedirs(charmmParamDir,exist_ok=True)


    charmmMol2 = p.join(charmmParamDir,f"{moleculeName}.mol2")
    amber_to_charmm_mol2(amberMol2, charmmMol2)

    charmmStr = p.join(charmmParamDir, f"{moleculeName}.str")
    amber_to_charmm_str(amberFrcmod, charmmMol2, charmmStr, moleculeName, config)

#####################################################


def amber_to_charmm_str(amberFrcmod, charmmMol2, charmmStr, moleculeName, config):
    charmmMol2Df, bonds = parse_mol2(charmmMol2)
    parsedFrcmod = parse_frcmod(amberFrcmod)
    write_charmm_str_file(charmmStr, moleculeName, charmmMol2Df, bonds, parsedFrcmod, config)
#####################################################
def write_str_header_section(f, moleculeName):
        f.write("* CHARMM Stream file generated from AMBER frcmod and mol2\n")
        f.write(f"* For residue: {moleculeName}\n")
        f.write("* Date: April 09, 2025\n\n")

#####################################################
def write_str_topology_section(f, moleculeName, charmmMol2Df, config):
    # Topology section
    totalCharge = config["moleculeInfo"]["charge"]
    f.write(f"RESI {moleculeName} {totalCharge:.6f}\n")

    chargeGroups = config["runtimeInfo"]["madeByCharges"]["chargeGroups"]

    print(chargeGroups)

    for chargeGroupName, chargeGroup in chargeGroups.items():
        if chargeGroupName.startswith("ACE_cap") or chargeGroupName.startswith("NME_cap"): 
            continue
        print(chargeGroupName)
        f.write("GROUP\n")
        groupDf = charmmMol2Df[charmmMol2Df["ATOM_ID"].isin(chargeGroup["indexes"])]
        for _, row in groupDf.iterrows():
            f.write(f"ATOM {row['ATOM_NAME']:<6} {row['ELEMENT']:<4} {row['CHARGE']:.6f}\n")
    f.write("\n")

#####################################################
def write_str_connect_section(f, charmmMol2Df, bonds, config):
        f.write("! Bonds from MOL2\n")
        for bond in bonds:
            atom1 = charmmMol2Df.loc[charmmMol2Df['ATOM_ID'] == bond[0], 'ATOM_NAME'].iloc[0]
            atom2 = charmmMol2Df.loc[charmmMol2Df['ATOM_ID'] == bond[1], 'ATOM_NAME'].iloc[0]
            f.write(f"BOND {atom1:<6} {atom2:<6}\n")

        nTerminalAtoms = config["moleculeInfo"]["nTermini"]
        cTerminalAtoms = config["moleculeInfo"]["cTermini"]
        for nTermini in nTerminalAtoms:
            f.write(f"BOND {nTermini:<6} {'-C':<6}\n")
        for cTermini in cTerminalAtoms:
            f.write(f"BOND {cTermini:<6} {'+N':<6}\n")
#####################################################

def write_str_bonds_section(f, parsedFrcmod, amberToCharmmTypes):
    # Parameter sections
    f.write("\nBONDS\n")
    for bond in parsedFrcmod['BONDS']:
        a1, a2 = bond['atoms']
        k = float(bond['k'])
        r0 = float(bond['r0'])
        charmmA1 = amberToCharmmTypes.get(a1.lower(), a1.upper())
        charmmA2 = amberToCharmmTypes.get(a2.lower(), a2.upper())
        f.write(f"{charmmA1:<4} {charmmA2:<4} {k:>8.2f} {r0:>6.4f}\n")


def write_str_angles_section(f, parsedFrcmod, amberToCharmmTypes):
    f.write("\nANGLES\n")
    for angle in parsedFrcmod['ANGLES']:
        a1, a2, a3 = angle['atoms']
        k = float(angle['k'])
        theta0 = float(angle['theta0'])
        charmmA1 = amberToCharmmTypes.get(a1.lower(), a1.upper())
        charmmA2 = amberToCharmmTypes.get(a2.lower(), a2.upper())
        charmmA3 = amberToCharmmTypes.get(a3.lower(), a3.upper())
        f.write(f"{charmmA1:<4} {charmmA2:<4} {charmmA3:<4} {k:>8.2f} {theta0:>6.2f}\n")

def write_str_dihedrals_section(f, parsedFrcmod, amberToCharmmTypes):
    f.write("\nDIHEDRALS\n")
    for dihedral in parsedFrcmod['DIHEDRALS']:
        a1, a2, a3, a4 = dihedral['atoms']
        k = float(dihedral['k'])
        n = abs(int(float(dihedral['periodicity'])))  # Convert to integer
        phi = float(dihedral['phase'])
        charmmA1 = amberToCharmmTypes.get(a1.lower(), a1.upper())
        charmmA2 = amberToCharmmTypes.get(a2.lower(), a2.upper())
        charmmA3 = amberToCharmmTypes.get(a3.lower(), a3.upper())
        charmmA4 = amberToCharmmTypes.get(a4.lower(), a4.upper())
        f.write(f"{charmmA1:<4} {charmmA2:<4} {charmmA3:<4} {charmmA4:<4} {k:>8.2f} {n:>2d} {phi:>6.2f}\n")


def write_str_impropers_section(f, parsedFrcmod, amberToCharmmTypes):
    f.write("\nIMPROPER\n")
    for improper in parsedFrcmod['IMPROPERS']:
        a1, a2, a3, a4 = improper['atoms']
        k = float(improper['k'])
        phi0 = float(improper['phi0'])
        n = int(float(improper['periodicity']))
        charmmA1 = amberToCharmmTypes.get(a1.lower(), a1.upper())
        charmmA2 = amberToCharmmTypes.get(a2.lower(), a2.upper())
        charmmA3 = amberToCharmmTypes.get(a3.lower(), a3.upper())
        charmmA4 = amberToCharmmTypes.get(a4.lower(), a4.upper())
        f.write(f"{charmmA1:<4} {charmmA2:<4} {charmmA3:<4} {charmmA4:<4} {k:>8.2f} {n:>2d} {phi0:>6.2f}\n")



#####################################################

def write_charmm_str_file(charmmStr, moleculeName, charmmMol2Df, bonds, parsedFrcmod, config):
    """Write CHARMM stream file with topology and parameters."""
    amberToCharmmTypes = init_amber_to_charmm_types()
    
    with open(charmmStr, 'w') as f:
        write_str_header_section(f, moleculeName)
        write_str_topology_section(f, moleculeName, charmmMol2Df, config)
        write_str_connect_section(f, charmmMol2Df, bonds, config)
        write_str_bonds_section(f, parsedFrcmod, amberToCharmmTypes)
        write_str_angles_section(f, parsedFrcmod, amberToCharmmTypes)
        write_str_dihedrals_section(f, parsedFrcmod, amberToCharmmTypes)
        write_str_impropers_section(f, parsedFrcmod, amberToCharmmTypes)
        f.write("END\n")
        f.write("RETURN\n")




def parse_mol2(mol2File):
    """Parse MOL2 file into a DataFrame for atoms and a list for bonds."""
    atomData = []
    bonds = []
    with open(mol2File, 'r') as f:
        inAtomSection = False
        inBondSection = False
        for line in f:
            if line.startswith('@<TRIPOS>ATOM'):
                inAtomSection = True
                continue
            elif line.startswith('@<TRIPOS>BOND'):
                inAtomSection = False
                inBondSection = True
                continue
            elif line.startswith('@<TRIPOS>'):
                inBondSection = False
                continue
            
            if inAtomSection and line.strip():
                parts = line.split()
                atomData.append({
                    'ATOM_ID': int(parts[0]),
                    'ATOM_NAME': parts[1],
                    'X': float(parts[2]),
                    'Y': float(parts[3]),
                    'Z': float(parts[4]),
                    'ELEMENT': parts[5],
                    'RES_ID': int(parts[6]),
                    'RES_NAME': parts[7],
                    'CHARGE': float(parts[8]) if len(parts) > 8 else 0.0
                })
            elif inBondSection and line.strip():
                parts = line.split()
                bonds.append((int(parts[1]), int(parts[2])))  # 1-based atom IDs
    
    # Convert atom data to DataFrame
    charmmMol2Df = pd.DataFrame(atomData)
    return charmmMol2Df, bonds


def amber_to_charmm_mol2(amberMol2, charmmMol2):

    amberMol2Df, bonds = parse_mol2(amberMol2)

    amberToCharmmTypes = init_amber_to_charmm_types()
    # Map AMBER atom types to elemental symbols
    amberMol2Df['ELEMENT'] = amberMol2Df['ELEMENT'].map(
        lambda x: amberToCharmmTypes.get(x, x[0].upper() if x[0].upper() in ['H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'] else x.upper())
    )

    # Open files for reading and writing
    with open(amberMol2, 'r') as fIn, open(charmmMol2, 'w') as fOut:
        inAtomSection = False
        atomIndex = 0  # To track which atom we're writing from the DataFrame
        
        for line in fIn:
            # Detect section start/end
            if line.startswith('@<TRIPOS>ATOM'):
                inAtomSection = True
                fOut.write(line)
                continue
            elif line.startswith('@<TRIPOS>'):
                inAtomSection = False
                fOut.write(line)
                continue

            if inAtomSection:
                # Write atom data from the DataFrame
                if atomIndex < len(amberMol2Df):
                    row = amberMol2Df.iloc[atomIndex]
                    fOut.write(f"{row['ATOM_ID']:<8} {row['ATOM_NAME']:<8} {row['X']:>9.6f} {row['Y']:>9.6f} {row['Z']:>9.6f} "
                               f"{row['ELEMENT']:<5} {row['RES_ID']:<5} {row['RES_NAME']:<8} {row['CHARGE']:>9.6f}\n")
                    atomIndex += 1
            else:
                fOut.write(line)  # Copy non-atom sections unchanged



def parse_frcmod(molFrcmod) -> dict:
    """
    Reads through a FRCMOD file and returns a dict with the contents of the file

    Args:
        molFrcmod (FilePath): frcmod file

    Returns:
        parsedFrcmod (dict): dict with the contents of the frcmod file
    
    """
    ## init a bunch of bools
    readingMass: bool       = False
    readingBonds: bool      = False
    readingAngles: bool     = False
    readingDihedrals: bool  = False
    readingImpropers: bool  = False
    readingNonbonded: bool  = False

    ## init empty parsedFrcmod dict
    parsedFrcmod = {"ATOMS": [], "BONDS": [], "ANGLES": [], "DIHEDRALS": [], "IMPROPERS": [], "NONBONDED": []}
    ## open molFrcmod for reading
    with open(molFrcmod, 'r') as f:
        ## loop through lines
        for line in f:
            ## skip empty lines
            if line.strip() == "":
                continue
            ## use bools to determine which section of the frcmod file we are in
            elif line.startswith("MASS"):
                readingMass = True
            elif line.startswith("BOND"):
                readingBonds = True
                readingMass = False
            elif line.startswith("ANGLE"):
                readingAngles = True
                readingBonds = False
            elif line.startswith("DIHE"):
                readingDihedrals = True
                readingAngles = False
            elif line.startswith("IMPROPER"):
                readingImpropers = True
                readingDihedrals = False
            elif line.startswith("NONBON"):
                readingNonbonded = True
                readingImpropers = False
            ## if we are in a data section, parse the line
            else:
                ## process mass data
                if readingMass:
                    lineData = line.split()
                    lineParsed = {"atomType": lineData[0], "mass": lineData[1], "vdw-radius": lineData[2]}
                    parsedFrcmod["ATOMS"].append(lineParsed)
                ## process bond data
                elif readingBonds:
                    atomData = get_frcmod_atom_data(line, [[0, 2], [4, 6]])
                    paramData = "".join(line[6:]).split()[0:2]
                    lineParsed = {"atoms": atomData, "r0": paramData[0], "k": paramData[1]}
                    parsedFrcmod["BONDS"].append(lineParsed)
                ## process angle data
                elif readingAngles:
                    atomData = get_frcmod_atom_data(line, [[0, 2], [3, 5], [6, 8]])
                    paramData = "".join(line[9:]).split()[0:2]
                    lineParsed = {"atoms": atomData, "theta0": paramData[0], "k": paramData[1]}
                    parsedFrcmod["ANGLES"].append(lineParsed)
                ## process dihedral data
                elif readingDihedrals:
                    atomData = get_frcmod_atom_data(line, [[0, 2], [3, 5], [6, 8], [9, 11]])
                    paramData = "".join(line[12:]).split()[0:4]
                    lineParsed = {"atoms": atomData, "multiplicity": paramData[0], "k": paramData[1], "phase": paramData[2], "periodicity": paramData[3]}  
                    parsedFrcmod["DIHEDRALS"].append(lineParsed)  
                ## process improper data
                elif readingImpropers:
                    atomData = get_frcmod_atom_data(line, [[0, 2], [3, 5], [6, 8], [9, 11]])
                    paramData = "".join(line[12:]).split()[0:3]
                    lineParsed = {"atoms": atomData, "k": paramData[0], "phi0": paramData[1], "periodicity": paramData[2]}
                    parsedFrcmod["IMPROPERS"].append(lineParsed)
                ## process nonbonded data
                elif readingNonbonded:
                    atomData = get_frcmod_atom_data(line, [[0, 6]])
                    paramData = "".join(line[6:]).split()[0:2]
                    lineParsed = {"atoms": atomData, "vdw-radius": paramData[0], "well-depth": paramData[1]}
                    parsedFrcmod["NONBONDED"].append(lineParsed)
                    
    return parsedFrcmod
def get_frcmod_atom_data(line: str, indexes: list) -> list:
    """
    Small function for getting atom names from FRCMOD file lines

    Args:
        line (str): line from a FRCMOD file
        indexes (List[int,int]): location of atom types in FRCMOD line

    Returns:
        atomData (List[str]): list of atom types 
    """

    return [line[start:end].strip() for start, end in indexes]



def  init_amber_to_charmm_types():
    amberToCharmmTypes = {
    # Carbons
    'c': 'C',    # Carbonyl sp²
    'c1': 'CY',  # Alkyne sp
    'c2': 'CA',  # Aromatic sp²
    'c3': 'CT3', # sp³ (adjust to CT1/CT2 based on H count)
    'ca': 'CA',  # Aromatic
    'cb': 'CA',  # Fused ring aromatic
    'cc': 'CPH1',# Imidazole-like
    'cd': 'CPH1',# Imidazole-like
    'cp': 'CA',  # Nucleic acid aromatic
    'cq': 'CA',  # Nucleic acid aromatic
    'cx': 'CT1', # sp³ cyclic (your data)
    'cy': 'CT2', # sp³ cyclic
    'CX': 'CT3' ## wildcard sp3
    # Hydrogens
    'h1': 'HA',  # Aliphatic CH
    'h2': 'HA',  # Aliphatic CH₂
    'h3': 'HA',  # Aliphatic CH₃
    'h4': 'HR1', # Aromatic (ortho)
    'h5': 'HR1', # Aromatic (meta)
    'ha': 'HP',  # Aromatic H
    'hc': 'HA',  # Aliphatic H (your data)
    'hn': 'H',   # H on N (your data)
    'ho': 'OH1', # Hydroxyl H
    'hp': 'HA',  # Polar H (adjust if needed)
    'hs': 'HS',  # Thiol H
    'hw': 'HT',  # Water H

    # Nitrogens
    'n': 'NH1',  # Amide sp²
    'n1': 'NH3', # Protonated amine
    'n2': 'NH2', # Guanidine NH₂
    'n3': 'NH3', # sp³ charged
    'na': 'NC2', # Aromatic/guadine N
    'nb': 'NR1', # Imidazole N
    'nc': 'NC2', # Guanidine N
    'nd': 'NR3', # Imidazole unprotonated
    'ns': 'NH1', # Peptide-like N (your data)

    # Oxygens
    'o': 'O',    # Carbonyl O (your data)
    'oh': 'OH1', # Hydroxyl O
    'os': 'OS',  # Ether O
    'ow': 'OT',  # Water O

    # Sulfur
    's': 'S',    # Thioether S (your data)
    'sh': 'SH1', # Thiol S
    'ss': 'S',   # Disulfide S

    # Others (add as needed)
    'f': 'F',    # Fluorine
    'cl': 'CL',  # Chlorine
    'br': 'BR',  # Bromine
    'i': 'I',    # Iodine
    'p2': 'P',   # Phosphorus
    'p3': 'P',   # Phosphorus
}
    return amberToCharmmTypes

if __name__ == "__main__":
    main()