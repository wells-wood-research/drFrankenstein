import subprocess

import numpy as np
import pandas as pd
import parmed
import rdkit
from rdkit import Chem


def generateTemplate(forcefield, topology, molecule):
    """
    Method to make OpenMM Templates for small molecules that do not
    have a template yet. This is useful for manually assigning charges
    and atom types.
    forcefield : OpenMM Forcefield
    topology : OpenMM Topology


    Returns
    -------
    OpenMM Forcefield object with new template registered and the topology
    with the atom types and charges added.
    """
    [templates, residues] = forcefield.generateTemplatesForUnmatchedResidues(topology)
    if len(templates) > 1:
        raise Exception(
            "Adding more than one new molecule class at a time. Not supported yet!!!"
        )
    for t in templates:
        if t is None:
            continue
        for index, atom in enumerate(t.atoms):
            atom.type = molecule.atom_types[index]
            atom.parameters = {"charge": molecule.charges[index]}
        forcefield.registerResidueTemplate(t)


def generate_amberparm(
    mol,
    new_charges=None,
    out_file="temp.prmtop",
    temp_filename="TEMP_FILENAME",
    antechamber_path="antechamber",
):
    Chem.MolToPDBFile(mol, temp_filename + ".pdb", flavor=2 | 8)
    antechamber_call = [
        antechamber_path,
        "-i",
        temp_filename + ".pdb",
        "-at",
        "amber",
        "-fi",
        "pdb",
        "-o",
        temp_filename + ".mol2",
        "-fo",
        "mol2",
        "-pf",
        "y",
    ]

    if new_charges:
        subprocess.check_output(antechamber_call)
        parm = parmed.load_file(temp_filename + ".mol2")
        for atom, charge in zip(parm.atoms, new_charges):
            atom.charge = charge
        parm.save(temp_filename + ".mol2")

    else:
        print("No charges provided !! Using antechamber's mopac charges mono-RESP ... ")
        antechamber_call.extend(["-c", "bcc"])
        subprocess.check_output(antechamber_call)

    # Use LEaP to create AMBER topology and coordinate files
    gen_unknown_dih = [
        "parmchk2",
        "-f",
        "mol2",
        "-i",
        f"{temp_filename}.mol2",
        "-o",
        f"{temp_filename}.frcmod",
    ]
    print(" ".join(gen_unknown_dih))
    subprocess.check_output(gen_unknown_dih)
    tleap_script = f"""source leaprc.protein.ff14SB
source leaprc.gaff
mol = loadmol2 {temp_filename}.mol2
loadamberparams {temp_filename}.frcmod
saveamberparm mol {out_file}.prmtop {temp_filename}.rst7
quit
"""
    with open("tleap.in", "w") as f:
        f.write(tleap_script)
    subprocess.check_output(["tleap", "-f", "tleap.in"])
    # remove files generated during the process
    remove_list = [
        f"{temp_filename}.frcmod",
        f"{temp_filename}.pdb",
        f"{temp_filename}.mol2",
        f"{temp_filename}.rst7",
        "tleap.in",
        "leap.log",
    ]
    for rm in remove_list:
        subprocess.call(["rm", rm])

    """
    # parse the atom types from the Antechamber output file; not simple, thus a helper method is called
    atom_types = align_antechamber(temp_filename + ".pdb", temp_filename + ".prepc")

    return atom_types
    """


def align_antechamber(pdb_path, antechamber_output):
    # !!!! code borrowed & slightly adapted from  psiomm
    # https://github.com/mzott/Psi4-OpenMM-Interface
    # Initial method written by Daniel Smith; @dgasmith
    # Read in both files
    # Skip first two lines where number of atoms and comment are
    with open(pdb_path, "r") as inp_o:
        inp_data = inp_o.readlines()
    bad_rows = [
        x for x in range(len(inp_data)) if len(inp_data[x].strip().split()) != 11
    ]
    print([len(inp_data[x].strip().split()) for x in range(len(inp_data))])

    inp = pd.read_csv(
        pdb_path,
        skiprows=bad_rows,
        sep=" +",
        usecols=[2, 5, 6, 7],
        engine="python",
        header=None,
    )
    inp.columns = ["Atom", "X", "Y", "Z"]
    print(inp)
    # Find which lines to skip - only lines with 8 columns are desired
    ant_o = open(antechamber_output, "r")
    ant_data = ant_o.readlines()
    ant_o.close()
    bad_rows = [
        x
        for x in range(len(ant_data))
        if len(ant_data[x].strip().split()) != 8
        or "DUMM" in ant_data[x].strip().split()
    ]
    out = pd.read_csv(
        antechamber_output, skiprows=bad_rows, sep=" +", engine="python", header=None
    )
    print(out)
    out.columns = ["#", "Sym", "GAFF", "??", "X", "Y", "Z", "Charge?"]

    # Parse out dummy atoms
    out = out[out["Sym"] != "DUMM"]

    ### Figure out mapping, easy way:
    inp_xyz = np.array(inp[["X", "Y", "Z"]])
    out_xyz = np.array(out[["X", "Y", "Z"]])

    inp_size = inp_xyz.shape[0]
    out_size = out_xyz.shape[0]

    # Sizes must be equal
    if inp_size != out_size:
        raise ValueError("Warning input and output sizes are not the same")

    translation_indices = np.zeros((inp_size), dtype=np.int)

    for irow in range(inp_size):
        found = 0
        idx = None
        for orow in range(out_size):
            # Careful! It looks like they round antechamer digits
            norm = np.linalg.norm(inp_xyz[irow] - out_xyz[orow])
            if norm < 1.0e-3:
                found += 1
                idx = orow

        # Should only be one
        if found > 1:
            raise ValueError("More than one output index matches a input index")
        if idx is None:
            raise ValueError("Did not find a matching output index for row %d" % irow)

        translation_indices[irow] = idx

    # Reorganize antechamber output
    out = out.iloc[translation_indices]
    out.reset_index()

    # The norm of this should now be small
    inp_xyz = np.array(inp[["X", "Y", "Z"]])
    out_xyz = np.array(out[["X", "Y", "Z"]])

    # Maximum norm per vector
    if np.max(np.linalg.norm(inp_xyz - out_xyz, axis=1)) > 0.002:
        raise ValueError("Warning alignment mismatch")

    inp["GAFF"] = out["GAFF"].values

    # inp.to_csv('parsed_output.csv')
    return np.asarray(inp["GAFF"])


def makeamberlibfile(frag, parm, to_removefromlib, resname="MOL"):
    tleapinput = open("tleap.in", "w")

    subprocess.check_output(
        "sed '/^DIHE/,/^NONBON/{/^NONBON/!d}'  parmchk2_all.frcmod > parmchk2.frcmod  "
    )
    tleapinput.write("loadAmberParams parmchk2.frcmod  \n")
    tleapinput.write("source leaprc.protein.ff14SB \n")
    tleapinput.write("%s = loadmol2 newcharges.mol2 \n" % (resname))
    for atom in parm.atoms:
        if atom.idx not in frag or atom.idx in to_removefromlib:
            tleapinput.write("remove %s %s.1.%s\n" % (resname, resname, atom.idx + 1))
    # for atom  in to_removefromlib :
    # tleapinput.write('remove %s %s.1.%s\n' %(resname, resname, atom) )
    tleapinput.write("set %s head %s.1.N \n" % (resname, resname))
    tleapinput.write("set %s tail %s.1.C \n" % (resname, resname))
    tleapinput.write(
        'impose %s { 1 2 3 } { { "N" "CA" "C" "O" 0.0 } { "H" "N" "CA" "C"  0.0 } { "CA" "CB" "CC" "CD" 180 } { "CB" "CC" "CD" "CE" 180 } } \n'
        % (resname)
    )
    tleapinput.write("saveoff  %s %s.lib \n" % (resname, resname))
    # os.system('grep -v ATOM  1.pdb > model.pdb')

    tleapinput.write("seq = sequence { ACE %s NME }\n" % (resname))
    tleapinput.write("savemol2 seq ACE-%s-NME.mol2 1 \n" % (resname))
    tleapinput.write(
        "loadoff %s/ethyl.lib  \n" % (os.path.dirname(os.path.abspath(__file__)))
    )
    tleapinput.write("seq = combine { ETH  seq } \n")
    tleapinput.write("bond seq.1.CY seq.3.CY \n")
    tleapinput.write("select seq.1 \n")
    tleapinput.write("relax seq.1 \n")
    tleapinput.write("select seq \n")
    tleapinput.write("relax seq\n")
    tleapinput.write("savemol2 seq ACE-%s-NME.mol2 1 \n" % (resname))
    tleapinput.write("quit")
    tleapinput.close()
    subprocess.check_output("tleap -f tleap.in")


if __name__ == "__main__":
    from prep import build_molecule_from_smiles

    smile = "CNC(=O)[C@H](CCC[C@@H](NC(C)=O)C(=O)NC)NC(C)=O"
    mol, res_type, backbone_list, capping_list = build_molecule_from_smiles(smile)
    charges = [
        -0.241323277603356,
        -0.3384318625324934,
        0.4204739692537924,
        -0.331553903434658,
        -0.0375018121722635,
        -0.0942023259998395,
        -0.028188063035473,
        -0.1826194690198855,
        -0.0483418380146683,
        -0.3994391944215769,
        0.4393446639761912,
        -0.3122127770065763,
        -0.3466226914287156,
        0.4921296243590897,
        -0.337861573928724,
        -0.3934136303195055,
        -0.3026992360245263,
        -0.3420313072472226,
        0.431309189076135,
        -0.3886477015866761,
        -0.3382307877072235,
        0.1071957446607865,
        0.1071957446607865,
        0.1071957446607865,
        0.2425899886438213,
        0.0430281549836637,
        0.0440602506950658,
        0.0440602506950658,
        0.0614350902329205,
        0.0614350902329205,
        0.0678996017381434,
        0.0678996017381434,
        0.0505736181011751,
        0.2481620217012617,
        0.0801380384110159,
        0.0801380384110159,
        0.0801380384110159,
        0.2839121814722846,
        0.127256632404196,
        0.127256632404196,
        0.127256632404196,
        0.222687057752786,
        0.0995166168009762,
        0.0995166168009762,
        0.0995166168009762,
    ]
    generate_amberparm(mol)
    generate_amberparm(
        mol,
        new_charges=charges,
        out_file="temp_newcharges",
    )
