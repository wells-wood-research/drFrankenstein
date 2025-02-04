import random
from itertools import chain

from rdkit import Chem


def random_renumber_atoms(mol):
    print(list(range(mol.GetNumAtoms())))
    order = list(range(mol.GetNumAtoms()))
    random.shuffle(order)
    print(order)
    renumbered_mol = Chem.rdmolops.RenumberAtoms(mol, order)

    for at in renumbered_mol.GetAtoms():
        print(at.GetIdx(), at.GetSymbol())
    return renumbered_mol


def renumber_frags(mol, fragments, frag_type, n_backbone):
    order_fragment = np.zeros(len(frag_type))
    outputs = []
    for i in range(n_backbone):
        flag = 0
        for i, type_frag in zip(order_fragment, frag_type):
            match type_frag:
                case "backbone":
                    if flag == i:
                        order_fragment[i] = 0
                        flag = None
                        # TODO write funstions to reoder etoms in part
                    else:
                        order_fragment[i] = 4
                        flag = i
                case "sidechain":
                    order_fragment[i] = 1
                case "capping":
                    order_fragment[i] = 2

        _, frag_type_reordered, fragments_reordered = (
            list(t) for t in zip(*sorted(zip(order_fragment, frag_type, fragments)))
        )

        new_order = list(chain(fragments_reordered))
        renumbered_mol = Chem.rdmolops.RenumberAtoms(mol, new_order)

        outputs.append((renumbered_mol, frag_type_reordered, fragments_reordered))
    return outputs


def rename_atoms(n_conf, Calpha, extin, extout=".mol2"):
    import string
    import sys

    import parmed

    ABC = list(string.ascii_uppercase)
    for i in range(1, n_conf + 1):
        # os.system('antechamber -i %s%s -fi pdb -o %sout%s  -fo pdb' %(i,extin,i,extin))

        parm = parmed.rdkit.load_rdkit(extin)
        AmberOrder = ["N", "H", "C", "O", "CA"]
        k = 1
        for atom in parm.atoms:
            if atom.name in AmberOrder:
                atom.name = "R" + str(k)
                atom.number = 55
            k += 1
        prev_order = [-1, -1, -1, -1, -1]
        current = []
        index = 1
        for atom in parm.atoms:
            if atom.name == "CA" or atom.idx == Calpha:
                atom.name = AmberOrder[4]
                # atom.type='CX'
                prev_order[4] = atom.idx
                atom.number = 4
                atom.type = "CX"
                for atm in atom.bond_partners:
                    if atm.element == 7:
                        prev_order[0] = atm.idx
                        atm.name = AmberOrder[0]
                        atom.number = 0
                        for at in atm.bond_partners:
                            if at.element == 1:
                                prev_order[1] = at.idx
                                at.name = AmberOrder[1]
                                at.number = 1
                    elif atm.element == 6:
                        if 8 in [at.element for at in atm.bond_partners]:
                            prev_order[2] = atm.idx
                            atm.name = AmberOrder[2]
                            atm.number = 2
                            for at in atm.bond_partners:
                                if at.element == 8:
                                    prev_order[3] = at.idx
                                    at.name = AmberOrder[3]
                                    at.number = 3
                        else:
                            current.append(atm)
                            t = 0

                            for at in current:
                                # print(at.name , at.idx , t,  len(current))
                                if len(current) == 1:
                                    t = ""
                                else:
                                    t += 1
                                at.name = "C" + str(t) + "B"
                                prev_order.append(at.idx)
                                at.number = 6
                                # print(at.name , at.idx)

                    elif atm.element == 1:
                        atm.number = 4
                        atm.name = "HA"

        if -1 in prev_order:
            parm.save("%s-debug.mol2" % i)
            sys.exit("debug: missing atoms in backbone conformation : %s " % i)
        if len(current) == 0:
            sys.exit("debug: no next atom")


def makeamberlibfile(frag, parm, to_removefromlib, resname="MOL"):
    import subprocess

    tleapinput = open("tleap.in", "w")
    subprocess.call(
        "parmchk2 -i %s  -o  parmchk2_all.frcmod -f  mol2 " % ("newcharges.mol2")
    )
    subprocess.call(
        "sed '/^DIHE/,/^NONBON/{/^NONBON/!d}'  parmchk2_all.frcmod > parmchk2.frcmod  "
    )
    tleapinput.write("loadAmberParams parmchk2.frcmod  \n")
    tleapinput.write("source leaprc.protein.ff14SB \n")
    tleapinput.write("%s = loadmol2 newcharges.mol2 \n" % (resname))
    for atom in parm.atoms:
        if atom.idx not in frag or atom.idx in to_removefromlib:
            tleapinput.write("remove %s %s.1.%s\n" % (resname, resname, atom.idx + 1))
    tleapinput.write("set %s head %s.1.N \n" % (resname, resname))
    tleapinput.write("set %s tail %s.1.C \n" % (resname, resname))
    tleapinput.write(
        'impose %s { 1 2 3 } { { "N" "CA" "C" "O" 0.0 } { "H" "N" "CA" "C"  0.0 } { "CA" "CB" "CC" "CD" 180 } { "CB" "CC" "CD" "CE" 180 } } \n'
        % (resname)
    )
    tleapinput.write("saveoff  %s %s.lib \n" % (resname, resname))
    tleapinput.write("quit")
    tleapinput.close()
    os.system("tleap -f tleap.in")
