
from pdbUtils import pdbUtils
import numpy as np
import networkx as nx   
import pandas as pd



## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass


def find_dihedrals(graph: nx.Graph) -> list[tuple[str, str, str, str]]:
    """Find all unique dihedrals in the graph and return them as a list of atom name tuples.

    Args:
        graph: NetworkX graph representing molecular structure.

    Returns:
        List of tuples, each containing four atom names defining a dihedral.
    """
    uniqueDihedralNodes = set()

    for u, v in graph.edges():
        neighborsU = [n for n in graph.neighbors(u) if n != v]
        neighborsV = [n for n in graph.neighbors(v) if n != u]

        for nU in neighborsU:
            for nV in neighborsV:
                if nU != nV:
                    if nU < nV:
                        uniqueDihedralNodes.add((nU, u, v, nV))
                    else:
                        uniqueDihedralNodes.add((nV, v, u, nU))
    
    dihedralAtomNames = []
    for n1, n2, n3, n4 in uniqueDihedralNodes:
        nameTuple = (
            graph.nodes[n1]['ATOM_NAME'],
            graph.nodes[n2]['ATOM_NAME'],
            graph.nodes[n3]['ATOM_NAME'],
            graph.nodes[n4]['ATOM_NAME']
        )
        dihedralAtomNames.append(nameTuple)

    return dihedralAtomNames

def pdb_to_graph(pdbFile: str) -> nx.Graph:
    """Convert a PDB file to a NetworkX graph.

    Args:
        pdbFile: Path to the input PDB file.

    Returns:
        NetworkX graph with nodes as atoms and edges as bonds within 1.6 Å.
    """
    df = pdbUtils.pdb2df(pdbFile)
    graph = nx.Graph()

    for _, row in df.iterrows():
        nodeId = f"{row['ATOM_NAME']}"
        graph.add_node(
            nodeId,
            ELEMENT=row["ELEMENT"],
            COORDS=np.array([row["X"], row["Y"], row["Z"]]),
            ATOM_NAME=row["ATOM_NAME"],
            ATOM_ID=row["ATOM_ID"]
        )

    for n1 in graph.nodes():
        for n2 in graph.nodes():
            if n1 < n2:
                dist = ((graph.nodes[n1]["COORDS"] - graph.nodes[n2]["COORDS"]) ** 2).sum() ** 0.5
                if dist < 1.6:
                    graph.add_edge(n1, n2, distance=dist)

    return graph

def calculate_dihedral(df: pd.DataFrame, atomNames: list[str]) -> float:
    """Calculate the dihedral angle for a sequence of four atoms.

    Args:
        df: DataFrame with atomic coordinates.
        atomNames: List or tuple of four atom names defining the dihedral.

    Returns:
        Dihedral angle in degrees.

    Raises:
        ValueError: If atomNames does not contain exactly four atom names.
    """
    if len(atomNames) != 4:
        raise ValueError("Please provide a list of four atomNames.")

    COORDS = get_coordinates(df, atomNames)
    
    b1 = COORDS[1] - COORDS[0]
    b2 = COORDS[2] - COORDS[1]
    b3 = COORDS[3] - COORDS[2]

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    
    n1Norm = n1 / np.linalg.norm(n1)
    n2Norm = n2 / np.linalg.norm(n2)

    dotProduct = np.dot(n1Norm, n2Norm)
    dotProduct = np.clip(dotProduct, -1.0, 1.0)

    angleRad = np.arccos(dotProduct)
    
    if np.dot(n1, b3) < 0:
        angleRad = -angleRad

    return np.degrees(angleRad)

def calculate_angle(df: pd.DataFrame, atomNames: list[str]) -> float:
    """Calculate the angle formed by three atoms.

    Args:
        df: DataFrame with atomic coordinates.
        atomNames: List or tuple of three atom names defining the angle.

    Returns:
        Angle in degrees.

    Raises:
        ValueError: If atomNames does not contain exactly three atom names.
    """
    if len(atomNames) != 3:
        raise ValueError("Please provide a list of three atomNames.")

    COORDS = get_coordinates(df, atomNames)
    
    vec1 = COORDS[0] - COORDS[1]
    vec2 = COORDS[2] - COORDS[1]

    vec1Norm = vec1 / np.linalg.norm(vec1)
    vec2Norm = vec2 / np.linalg.norm(vec2)
    
    dotProduct = np.dot(vec1Norm, vec2Norm)
    dotProduct = np.clip(dotProduct, -1.0, 1.0)
    
    angleRad = np.arccos(dotProduct)
    
    return np.degrees(angleRad)

def calculate_bond_length(df: pd.DataFrame, atomNames: list[str]) -> float:
    """Calculate the bond length between two atoms.

    Args:
        df: DataFrame with atomic coordinates.
        atomNames: List or tuple containing two atom names.

    Returns:
        Bond length as a float.

    Raises:
        ValueError: If atomNames does not contain exactly two atom names.
    """
    if len(atomNames) != 2:
        raise ValueError("Please provide a list of two atomNames.")

    COORDS = get_coordinates(df, atomNames)
    return np.linalg.norm(COORDS[0] - COORDS[1])

def get_coordinates(df: pd.DataFrame, atomNames: list[str]) -> np.ndarray:
    """Retrieve XYZ coordinates for a list of atom names from the DataFrame.

    Args:
        df: DataFrame containing atomic information.
        atomNames: List of atom names.

    Returns:
        NumPy array of XYZ coordinates for the specified atoms.
    """
    COORDS = []
    for atomId in atomNames:
        COORDS.append(df[df['ATOM_NAME'] == atomId][['X', 'Y', 'Z']].values[0])
    return np.array(COORDS)

def gen_internal_coords(inPdb: FilePath) -> str:
    """Generate internal coordinate (IC) data from a PDB file and print it.

    The output format is:
    IC atom1 atom2 atom3 atom4 bond12 angle123 dihedral1234 angle234 bond34
    """
    graph = pdb_to_graph(inPdb)
    dihedrals = find_dihedrals(graph)
    pdbDf = pdbUtils.pdb2df(inPdb)

    capsToBackbone = {
        "NN": "+N",
        "HNN1": "+NH",
        "CN": "+CA",
        "CC1": "-C",
        "OC": "-O",
        "CC2": "-CA"
    }

    ignoreAtoms = ["HC1", "HC2", "HC3", "HCN1", "HCN2", "HCN3"]

    icData = []
    for dihedral in dihedrals:
        if any(atom in ignoreAtoms for atom in dihedral):
            continue
        bond12 = calculate_bond_length(pdbDf, dihedral[:2])
        angle123 = calculate_angle(pdbDf, dihedral[:3])
        dihedral1234 = calculate_dihedral(pdbDf, dihedral)
        angle234 = calculate_angle(pdbDf, dihedral[1:])
        bond34 = calculate_bond_length(pdbDf, dihedral[2:])
        atomsRenamed = [capsToBackbone.get(atom, atom) for atom in dihedral]
        icLine = f"IC {atomsRenamed[0]}\t{atomsRenamed[1]}\t{atomsRenamed[2]}\t{atomsRenamed[3]}\t{bond12:.4f}  {angle123:.4f}  {dihedral1234:.4f}\t{angle234:.4f}  {bond34:.4f}"
        icData.append(icLine)
    icBlock = "\n".join(icData)
    return icBlock

