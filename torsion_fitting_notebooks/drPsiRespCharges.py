import psiresp
from rdkit import Chem
from rdkit.Chem import AllChem
from subprocess import call
from os import path as p
import os

'''
1. Load Pdb file from inputDir
2. Read Pdb as a rdkit molecule
3. read rdkit molecule as psi4 molecule
4. set up 

'''
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
## dummy classes
class FilePath:
    pass
class DirPath:
    pass
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def dummy_input():
    cappedDir = "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs/capped_amino_acids"
    molPdb = p.join(cappedDir, "ALA_capped.pdb")

    topDir = "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs"
    outDir = "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs/charge_calculations"

    return cappedDir, molPdb, topDir, outDir
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def main():
    cappedDir, molPdb, topDir, outDir = dummy_input()
    os.makedirs(outDir, exist_ok=True)
    psiMol = load_molecule(molPdb)
    cappingGroupIndexes = get_capping_group_indexes(psiMol)

    psiRespJob = set_up_psiresp_job(psiMol, cappingGroupIndexes, outDir)
    charges = run_psiresp_job(psiRespJob, topDir)
    print(charges)
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_capping_group_indexes(psiMol):
    rdkitMol = psiMol.to_rdkit()
    backBoneAndCapsSmarts = "[#6]([H])([H])([H])[#7][#6](=[#8])-[#6@@](-[H])-[#7]([H])-[#6](=[#8])-[#6]([H])([H])([H])"

    print(Chem.MolToSmarts(rdkitMol))
    # Find the substructure match
    backBoneAndCapsPattern = Chem.MolFromSmarts(backBoneAndCapsSmarts)
    backBoneAndCapsMatches = rdkitMol.GetSubstructMatches(backBoneAndCapsPattern)


    print(len(backBoneAndCapsMatches[0]))
    matches = [list(match) for match in backBoneAndCapsMatches]
    nMethylCapIndexes = matches[0][0:6]
    acylCapIndexes = matches[0][-6:]

    cappingGroupIndexes = nMethylCapIndexes + acylCapIndexes

    return cappingGroupIndexes
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_psiresp_job(psiRespJob, topDir, outDir):
    os.chdir(topDir)

    psiRespJob.generate_conformers()
    psiRespJob.generate_orientations()
    ###
    print("Making Optimization Job!")
    run_psiresp_job_and_catch_exit(psiRespJob)
    ###
    print("Running Optimization Step!")
    optDir = p.join(outDir, "optimization")
    os.chdir(topDir)
    call(["chmod", "+x", "run_optimization.sh"])
    call(["bash", "run_optimization.sh"])
    ###
    ##
    print("Making Single Point Job!")
    os.chdir(outDir)
    run_psiresp_job_and_catch_exit(psiRespJob)
    stageOneMolecule = psiRespJob.molecules[0]
    print(stageOneMolecule)
    ###
    print("Running Single Point Step!")
    singlePointDir = p.join(outDir, "single_point")
    os.chdir(singlePointDir)
    call(["chmod", "+x", "run_single_point.sh"])
    call(["bash", "run_single_point.sh"])
    ###
    print("making charge calculation job!")
    os.chdir(topDir)
    run_psiresp_job_and_catch_exit(psiRespJob)
    psiRespJob.compute_charges()

    charges = psiRespJob.molecules[0].stage_2_restrained_charges
    return charges

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_psiresp_job_and_catch_exit(psiRespJob):
    """
    Runs a psiresp job and catches the SystemExit error if it occurs.
    Allows us to run the psiresp job without it exiting
    """
    try:
        psiRespJob.run()
    except SystemExit as e:
        pass
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def set_up_psiresp_job(psiMol, cappingGroupIndexes, outDir):
    """
    Set up a psiresp job from a psi4 molecule.

    Parameters
    ----------
    psiMol : psiresp.Molecule
        The psi4 molecule to set up the job for.
    cappingGroupIndexes : list of int
        The indices of the atoms in the capping group.
    outDir : DirPath
        The directory to save the job files in.

    Returns
    -------
    psiRespJob : psiresp.Job
        The set up job.
    """
    constraints = psiresp.ChargeConstraintOptions(symmetric_atoms_are_equivalent=True)
    constraints.add_charge_sum_constraint_for_molecule(psiMol, charge=0, indices=cappingGroupIndexes)

    psiMol.optimize_geometry=True

    geometry_options = psiresp.QMGeometryOptimizationOptions(
    method="hf",
    basis="6-31g*")

    esp_options = psiresp.QMEnergyOptions(
    method="hf",
    basis="6-31g*")

    conformerGenerationOptions = psiresp.ConformerGenerationOptions(
        n_max_conformers=10
    )


    psiRespJob = psiresp.Job(
        molecules=[psiMol],
        charge_constraints=constraints,
        qm_optimization_options=geometry_options,
        qm_esp_options=esp_options,
        working_directory=outDir,
        n_processes = 30,
        conformer_generation_options=conformerGenerationOptions
    )
    
    return psiRespJob

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def load_molecule(molPdb: FilePath) -> psiresp.Molecule:
    rdkitMol = Chem.MolFromPDBFile(molPdb, removeHs=True)
    rdkitMol = enforce_carboxyl_bonding(rdkitMol)


    inputSmiles = Chem.MolToSmiles(rdkitMol)

    rdkitMol = Chem.MolFromSmiles(inputSmiles)

    Chem.SanitizeMol(rdkitMol)
    rdkitMol = Chem.AddHs(rdkitMol)

    AllChem.EmbedMultipleConfs(rdkitMol, numConfs=10, randomSeed=42)
    AllChem.UFFOptimizeMolecule(rdkitMol, maxIters=1000)
    psiMol = psiresp.Molecule.from_rdkit(rdkitMol)

    return psiMol

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def enforce_carboxyl_bonding(mol):
    # Define the initial and target SMARTS patterns
    amideSmarts = '[C]([O])([N])([C])'
    print(Chem.MolToSmarts(mol))

    # Find the substructure match
    amidePattern = Chem.MolFromSmarts(amideSmarts)
    matches = mol.GetSubstructMatches(amidePattern)
    
    if not matches:
        raise ValueError("No matching substructure found.")
    
    # Modify the first match found
    emol = Chem.EditableMol(mol)
    for match in matches:
        # Get the indices of the atoms in the match
        carbon_idx = match[0]
        oxygen_idx = match[1]        
        # Remove the single bond and add a double bond
        emol.RemoveBond(carbon_idx, oxygen_idx)
        emol.AddBond(carbon_idx, oxygen_idx, Chem.BondType.DOUBLE)
    mol = emol.GetMol()
    print(Chem.MolToSmarts(mol))

    # Return the modified molecule
    return emol.GetMol()    
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

if __name__ == "__main__":
    main()
