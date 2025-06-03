## BASIC IMPORTS ##
from pdbUtils import pdbUtils
import numpy as np

## GRAPH-BASED IMPORTS ##
import networkx as nx
from networkx.algorithms import isomorphism

## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass
from typing import List, Tuple


def symmetry_protocol(pdbFile: FilePath) -> List[List[str]]:
    """
    Main protocol for symmetry analysis
    1. Convert PDB to graph
    2. Find symmetrically equivalent atoms
    """

    moleculeGraph = pdb2graph(pdbFile)
    _, symmetricIds = find_symmetric_atoms(moleculeGraph)

    return symmetricIds


def pdb2graph(pdbFile: FilePath) -> nx.Graph:
    """
    Convert a PDB file to a NetworkX graph.
    """
    # Manual construction from PDB file using DataFrame
    df = pdbUtils.pdb2df(pdbFile)
    graph = nx.Graph()

    # Add all nodes (including hydrogens)
    for _, row in df.iterrows():
        nodeId = f"{row['CHAIN_ID']}:{row['ATOM_ID']-1}"  # 0-based indexing
        graph.add_node(
            nodeId,
            element=row["ELEMENT"],
            coords=np.array([row["X"], row["Y"], row["Z"]]),
            atom_name=row["ATOM_NAME"],
            atom_index= row["ATOM_ID"]
        )

    # Add edges with 1.6 Å threshold
    for n1 in graph.nodes():
        for n2 in graph.nodes():
            if n1 < n2:
                dist = ((graph.nodes[n1]["coords"] - graph.nodes[n2]["coords"]) ** 2).sum() ** 0.5
                if dist < 1.6:
                    graph.add_edge(n1, n2, distance=dist)

    return graph

def find_symmetric_atoms(graph: nx.Graph) -> tuple[List[List[str]], List[List[str]]]:
    """
    Identify symmetrically equivalent atoms in a molecular graph using NetworkX.
    Returns a list of lists, where each sublist contains the ATOM_NAME of equivalent atoms.
    """
    # Define a node match function to ensure atoms have the same element
    nodeMatch = lambda n1, n2: n1["element"] == n2["element"]

    # Find all automorphisms (symmetry operations) of the graph
    gm = isomorphism.GraphMatcher(graph, graph, node_match=nodeMatch)
    automorphisms = list(gm.isomorphisms_iter())

    if not automorphisms:
        print("No symmetry found in the graph.")
        # Return all atom names as a single group if no symmetry
        return [[graph.nodes[n]["atom_name"] for n in graph.nodes()]]

    # Group atoms by their equivalence under automorphisms
    symmetryGroups = {}
    for node in graph.nodes():
        equivSet = set()
        for mapping in automorphisms:
            equivSet.add(mapping[node])
        symmetryGroups[frozenset(equivSet)] = list(equivSet)

    # Convert to a list of ATOM_NAMEs
    symmetricAtoms = []
    symmetricIds = []
    for group in symmetryGroups.values():
        # Extract ATOM_NAME directly from node attributes
        atomNames = [graph.nodes[n]["atom_name"] for n in group]
        atomIds = [graph.nodes[n]["atom_index"] for n in group]
        symmetricAtoms.append(atomNames)
        symmetricIds.append(atomIds)


    return symmetricAtoms, symmetricIds

if __name__ == "__main__":
    raise NotImplementedError