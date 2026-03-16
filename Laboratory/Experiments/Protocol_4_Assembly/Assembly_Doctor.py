import os
from os import path as p


from . import Assembly_Monster
# from . import Assembly_Assistant
# Placeholder classes (extend if needed)
class FilePath:
    pass

class DirectoryPath:
    pass


def find_hess_xyz_files(config:dict) -> list:
    solvatedDir = config["runtimeInfo"]["madeByCharges"]["solvatedDir"]

    orcaDir = p.join(solvatedDir, "01_orca_calculations")

    hessianFiles = [p.join(orcaDir, dirname, "orca_opt_freq.hess") for dirname in os.listdir(orcaDir) if os.path.isdir(p.join(orcaDir, dirname))]

    xyzFiles = [p.join(orcaDir, dirname, "orca_opt_freq.xyz") for dirname in os.listdir(orcaDir) if os.path.isdir(p.join(orcaDir, dirname))]
    return hessianFiles, xyzFiles

def agnostic_assembly_protocol(config:dict) -> dict:
    """
    Main protocol for assembly of starting parameters for agnostic parameters
    This protocol:
    [] Find .hess files from charges calculation (might depend on protocol used for charges)
    [] loop through .hess files
        [] construct hessian matrix from .hess file
        [] invert to get compliance matrix
        [] construct adjacancy and angle matrix frpm cappedPdb
        [] convert adjacency and angle matrix to bond and angle lists
        [] calculate distance and angle strengths (k values)
        [] extract optimised bond and angle values from cappedPdb

    [] average over r0 and k0 values for each bond and angle
    [] create a dict with atomIds as keys and:
        [] Element
        [] nBonds
        [] Charge
        [] r0 values
        [] k0 values
    
    """



    ## create an empty dict in config
    config["runtimeInfo"]["madeByAssembly"] = {}

    ## make a dir for assembly of initial parameters
    assemblyDir = p.join(config["pathInfo"]["outputDir"], "02_parameter_assembly")
    os.makedirs(assemblyDir, exist_ok=True)
    config["runtimeInfo"]["madeByAssembly"]["assemblyDir"] = assemblyDir

    ## find hessian files
    hessianFiles, xyzFiles = find_hess_xyz_files(config)

    bonds, angles, dihedrals, impropers = Assembly_Monster.get_molecule_connectivity(config)

    allKbValues = []
    allKThetaValues = []
    allr0values = []
    allTheta0values = []
    for hessianFile, xyzFile in zip(hessianFiles, xyzFiles):
        kbValues, r0values, kThetaValues,  theta0values = Assembly_Monster.process_hessian(hessianFile=hessianFile,
                                                                                                bonds=bonds,
                                                                                                angles=angles,
                                                                                                xyzFile=xyzFile,
                                                                                                config=config)
        allKbValues.append(kbValues)
        allr0values.append(r0values)
        allKThetaValues.append(kThetaValues)
        allTheta0values.append(theta0values)

    kbValues = Assembly_Monster.average_dict_values(allKbValues)
    r0values = Assembly_Monster.average_dict_values(allr0values)
    kThetaValues = Assembly_Monster.average_dict_values(allKThetaValues)
    theta0values = Assembly_Monster.average_dict_values(allTheta0values)

    atomFeaturesDf, atomTypesDict = Assembly_Monster.construct_atom_features(kbValues, r0values, kThetaValues, theta0values, config=config)

    atomTypes, atomTypesDf = Assembly_Monster.assign_atom_types_by_clustering(atomFeaturesDf)



    (bondParamsByType,
      angleParamsByType,
        groupedDihedrals,
          groupedImpropers) = Assembly_Monster.group_params_by_atom_type(atomTypes,
                                                                          kbValues, r0values,
                                                                            kThetaValues, theta0values, 
                                                                                dihedrals, impropers)
    
    nonBondedParams = Assembly_Monster.assign_non_bonded_by_analogy(atomTypesDf)

    config = Assembly_Monster.construct_frcmod(atomTypes,
                                        bondParamsByType, 
                                        angleParamsByType,
                                           groupedDihedrals, groupedImpropers,
                                           nonBondedParams,
                                                 config=config)

    config = Assembly_Monster.construct_mol2(atomTypesDf, config=config)


    return config