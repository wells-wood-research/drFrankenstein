from os import path as p
import pandas as pd
import numpy as np
from pdbUtils import pdbUtils
from shutil import move, copy
import parmed
from parmed.charmm import CharmmParameterSet
from parmed.topologyobjects import Cmap
import openmm.app as app
import openmm as openmm
import openmm.unit as unit
from psfgen import PsfGen

class FilePath:
    pass

class DirectoryPath:
    pass

def add_cmap_term(config: dict) -> dict:
    """
    Adds CMAP term to a CHARMM RTF file

    Args: 
        config (dict): contains all run info
    Returns:
        config (dict): updated config
    
    """
    moleculeRtf = config['runtimeInfo']['madeByStitching']['moleculeRtf']
    moleculeName = config['moleculeInfo']['moleculeName']
    parmedRtf = CharmmParameterSet(moleculeRtf)

def copy_assembled_parameters(config: dict) -> dict:
    assembledRtf = config['runtimeInfo']['madeByAssembly']['assembledRtf']
    assembledPrm = config['runtimeInfo']['madeByAssembly']['assembledPrm']
    assembledPsf = config['runtimeInfo']['madeByAssembly']['assembledPsf']
    moleculeParameterDir = config['runtimeInfo']['madeByStitching']['moleculeParameterDir']
    moleculeName = config['moleculeInfo']['moleculeName']
    destRtf = p.join(moleculeParameterDir, f'DEFAULT_PARAMS.rtf')
    copy(assembledRtf, destRtf)
    destPrm = p.join(moleculeParameterDir, f'DEFAULT_PARAMS.prm')
    copy(assembledPrm, destPrm)
    destPsf = p.join(moleculeParameterDir, f'DEFAULT_PARAMS.psf')
    copy(assembledPsf, destPsf)
    config['runtimeInfo']['madeByStitching']['moleculeRtf'] = destRtf
    config['runtimeInfo']['madeByStitching']['moleculePrm'] = destPrm
    config['runtimeInfo']['madeByStitching']['moleculePsf'] = destPsf
    return config

def update_prm(config: dict, torsionTag: str, torsionParamDf: pd.DataFrame, shuffleIndex: int) -> dict:
    """
    Uses Parmed to update a CHARMM PRM file with a torsion parameter block 
    
    """
    moleculePrm = config['runtimeInfo']['madeByStitching']['moleculePrm']
    moleculeRtf = config['runtimeInfo']['madeByStitching']['moleculeRtf']
    rotatableDihedrals = config['runtimeInfo']['madeByTwisting']['rotatableDihedrals']
    moleculeParameterDir = config['runtimeInfo']['madeByStitching']['moleculeParameterDir']
    moleculeName = config['moleculeInfo']['moleculeName']
    targetAtomTypes = tuple(rotatableDihedrals[torsionTag]['ATOM_TYPES'])
    newDihedralTypes = [parmed.DihedralType(phi_k=row['Amplitude'], per=int(row['Period']), phase=row['Phase']) for (_, row) in torsionParamDf.iterrows()]
    parmedPrm = CharmmParameterSet(moleculePrm, moleculeRtf)
    dihedralTypesToRemove = []
    for prmDihedralTypes in parmedPrm.dihedral_types:
        if prmDihedralTypes == targetAtomTypes or prmDihedralTypes == targetAtomTypes[::-1]:
            dihedralTypesToRemove.append(prmDihedralTypes)
    for dihedralTypeToRemove in dihedralTypesToRemove:
        del parmedPrm.dihedral_types[dihedralTypeToRemove]
    parmedPrm.dihedral_types[targetAtomTypes] = newDihedralTypes
    proposedPrm = p.join(moleculeParameterDir, f'{moleculeName}_capped_{shuffleIndex + 1}.prm')
    parmedPrm.write(par=proposedPrm)
    config['runtimeInfo']['madeByStitching']['proposedPrm'] = proposedPrm
    return config

def parse_prm(moleculePrm: FilePath) -> dict:
    """
    Reads through a CHARMM 
    """
    readingBonds: bool = False
    readingAngles: bool = False
    readingDihedrals: bool = False
    readingImpropers: bool = False
    parsedPrm = {'BONDS': [], 'ANGLES': [], 'DIHEDRALS': [], 'IMPROPERS': []}
    with open(moleculePrm, 'r') as f:
        for line in f.readlines():
            if line.strip() == '':
                continue
            if line.startswith('BONDS'):
                readingBonds = True
            elif line.startswith('ANGLES'):
                readingBonds = False
                readingAngles = True
            elif line.startswith('DIHEDRALS'):
                readingAngles = False
                readingDihedrals = True
            elif line.startswith('IMPROPERS'):
                readingDihedrals = False
                readingImpropers = True
            elif line.startswith('END'):
                readingImpropers = False
            else:
                lineData = line.split()
                if '!' in line:
                    comment = '!' + line.split('!')[1].strip()
                else:
                    comment = ''
                if readingBonds:
                    lineParsed = {'atoms': [lineData[0], lineData[1]], 'k': lineData[2], 'r0': lineData[3], 'comment': comment}
                    parsedPrm['BONDS'].append(lineParsed)
                elif readingAngles:
                    lineParsed = {'atoms': [lineData[0], lineData[1], lineData[2]], 'k': lineData[3], 'theta0': lineData[3], 'comment': comment}
                    parsedPrm['ANGLES'].append(lineParsed)
                elif readingDihedrals:
                    lineParsed = {'atoms': [lineData[0], lineData[1], lineData[2], lineData[3]], 'k': lineData[4], 'period': lineData[5], 'phase': lineData[6], 'comment': comment}
                    parsedPrm['DIHEDRALS'].append(lineParsed)
                elif readingImpropers:
                    lineParsed = {'atoms': [lineData[0], lineData[1], lineData[2], lineData[3]], 'k': lineData[4], 'period': lineData[5], 'phase': lineData[6], 'comment': comment}
                    parsedPrm['IMPROPERS'].append(lineParsed)
    return parsedPrm

def write_prm(parsedPrm: dict, outPrm: FilePath) -> None:
    header = '\nread param card flex append\n* Parameters generated by analogy by\n* CHARMM General Force Field (CGenFF) program version 4.0\n* And modified by drFRANKENSTEIN\n*\n\n! Penalties lower than 10 indicate the analogy is fair; penalties between 10\n! and 50 mean some basic validation is recommended; penalties higher than\n! 50 indicate poor analogy and mandate extensive validation/optimization.\n\n\n'
    bondsSection = ''
    for param in parsedPrm['BONDS']:
        bondsSection += f"{param['atoms'][0]:<7}{param['atoms'][1]:<8}{param['k']:<11}{param['r0']:<7} {param['comment']}\n"
    anglesSection = ''
    for param in parsedPrm['ANGLES']:
        anglesSection += f"{param['atoms'][0]:<7}{param['atoms'][1]:<7}{param['atoms'][2]:<7}"
        anglesSection += f"{param['k']:>7}{param['theta0']:>10} {param['comment']}\n"
    dihedralSection = ''
    for param in parsedPrm['DIHEDRALS']:
        dihedralSection += f"{param['atoms'][0]:<7}{param['atoms'][1]:<7}{param['atoms'][2]:<7}{param['atoms'][3]:<7}"
        dihedralSection += f"{param['k']:>10}{param['period']:>3}{param['phase']:>9} {param['comment']}\n"
    improperSection = ''
    for param in parsedPrm['IMPROPERS']:
        improperSection += f"{param['atoms'][0]:<7}{param['atoms'][1]:<7}{param['atoms'][2]:<7}{param['atoms'][3]:<7}"
        improperSection += f"{param['k']:>10}{param['period']:>3}{param['phase']:>9} {param['comment']}\n"
    with open(outPrm, 'w') as f:
        f.write(header)
        f.write('BONDS\n')
        f.write(bondsSection + '\n')
        f.write('ANGLES\n')
        f.write(anglesSection + '\n')
        f.write('DIHEDRALS\n')
        f.write(dihedralSection + '\n')
        f.write('IMPROPERS\n')
        f.write(improperSection + '\n')
        f.write('END\nRETURN')

def set_capping_types_to_charmm_defaults(config: dict) -> None:
    """
    Sets the capping groups to CHARMM defaults for backbone atoms

    Args:
        config (dict): contains all run info
    Returns:
        None
    """
    mapCgenffTypeToCharmmType = set_capping_types_to_charmm_defaults_rtf(config)
    set_capping_types_to_charmm_defaults_prm(config, mapCgenffTypeToCharmmType)

def set_capping_types_to_charmm_defaults_prm(config: dict, mapCgenffTypeToCharmmType: dict) -> None:
    moleculePrm = config['runtimeInfo']['madeByStitching']['moleculePrm']
    mmTorsionCalculationDir = config['runtimeInfo']['madeByStitching']['mmTorsionCalculationDir']
    tmpPrm = p.join(mmTorsionCalculationDir, 'tmp.prm')
    parsedPrm = parse_prm(moleculePrm)
    for (cgenffType, charmmType) in mapCgenffTypeToCharmmType.items():
        for paramEntry in ['BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS']:
            for param in parsedPrm[paramEntry]:
                try:
                    param['atoms'][param['atoms'].index(cgenffType)] = charmmType
                except ValueError:
                    pass
    write_prm(parsedPrm, tmpPrm)
    move(tmpPrm, moleculePrm)

def set_capping_types_to_charmm_defaults_rtf(config: dict) -> dict:
    moleculeRtf = config['runtimeInfo']['madeByStitching']['moleculeRtf']
    mmTorsionCalculationDir = config['runtimeInfo']['madeByStitching']['mmTorsionCalculationDir']
    mapAtomNameToDefaultCharmmType = {'NN': 'NH1', 'CN': 'CT1', 'HNN1': 'H', 'OC': 'O', 'CC1': 'C', 'CC2': 'CT1'}
    mapCgenffTypeToCharmmType = {}
    tmpRtf = p.join(mmTorsionCalculationDir, 'tmp.rtf')
    with open(moleculeRtf, 'r') as inRtf, open(tmpRtf, 'w') as outRtf:
        for line in inRtf.readlines():
            if line.startswith('ATOM'):
                lineData = line.split()
                if lineData[1] in mapAtomNameToDefaultCharmmType:
                    mapCgenffTypeToCharmmType[lineData[2]] = mapAtomNameToDefaultCharmmType[lineData[1]]
                    lineData[2] = mapAtomNameToDefaultCharmmType[lineData[1]]
                line = '\t'.join(lineData) + '\n'
            outRtf.write(line)
    move(tmpRtf, moleculeRtf)
    return mapCgenffTypeToCharmmType

def create_atom_type_map(config: dict) -> dict:
    """
    Creates a map of atom types for CHARMM given a RTF file

    Args: 
        config (dict): contains all run info
    Returns:
        config (dict): updated config
    """
    moleculeRtf = config['runtimeInfo']['madeByStitching']['moleculeRtf']
    atomTypeMap = {}
    with open(moleculeRtf, 'r') as inRtf:
        for line in inRtf.readlines():
            if line.startswith('ATOM'):
                lineData = line.split()
                atomTypeMap[lineData[1]] = lineData[2]
    config['runtimeInfo']['madeByStitching']['atomTypeMap'] = atomTypeMap
    return config

def edit_rtf_charges(config: dict) -> dict:
    """
    Edits an RTF file to change the charges to those calculated by drFrankenstein
    Args: 
        config (dict): contains all run info
    Returns:
        config (dict): updated config 
    """
    moleculeRtf = config['runtimeInfo']['madeByStitching']['moleculeRtf']
    chargesCsv = config['runtimeInfo']['madeByCharges']['chargesCsv']
    mmTotalCalculationDir = config['runtimeInfo']['madeByStitching']['mmTotalCalculationDir']
    chargesDf = pd.read_csv(chargesCsv, index_col='Unnamed: 0')
    atomIndex = 0
    tmpRtf = p.join(mmTotalCalculationDir, 'update_charges.rtf')
    with open(moleculeRtf, 'r') as inRtf, open(tmpRtf, 'w') as outRtf:
        for line in inRtf.readlines():
            if line.startswith('ATOM'):
                lineData = line.split()
                lineData[3] = round(chargesDf.loc[atomIndex, 'Charge'], 3).astype(str)
                line = '\t'.join(lineData) + '\n'
                atomIndex += 1
            outRtf.write(line)
    move(tmpRtf, moleculeRtf)

def make_charmm_psf(config: dict) -> dict:
    """
    Uses PSFGen to generate a PSF file from RTF and PRM files

    Args: 
        config (dict): contains all run info
    Returns:
        config (dict): updated config 
    """
    moleculeRtf = config['runtimeInfo']['madeByStitching']['moleculeRtf']
    cappedPdb = config['runtimeInfo']['madeByCapping']['cappedPdb']
    mmTorsionCalculationDir = config['runtimeInfo']['madeByStitching']['mmTorsionCalculationDir']
    moleculeName = config['moleculeInfo']['moleculeName']
    cappedDf = pdbUtils.pdb2df(cappedPdb)
    cappedDf['RES_ID'] = '1'
    cappedDf['RES_NAME'] = moleculeName
    tmpPdb = p.join(mmTorsionCalculationDir, f'united_capped_{moleculeName}.pdb')
    pdbUtils.df2pdb(cappedDf, tmpPdb)
    gen = PsfGen(output='/dev/null')
    gen.read_topology(moleculeRtf)
    gen.add_segment(segid='A', pdbfile=tmpPdb)
    gen.read_coords(segid='A', filename=tmpPdb)
    psfFile = p.join(mmTorsionCalculationDir, f'{moleculeName}_capped.psf')
    gen.write_psf(psfFile)
    config['runtimeInfo']['madeByStitching']['moleculePsf'] = psfFile
    return config

def find_cgenff_params(config: dict) -> dict:
    """
    Looks inside the drFrankenstein dir to find CGenFF PRM and RTF files

    Args: 
        config (dict): contains all run info
    Returns:
        config (dict): updated config 
    """
    stitchingSrcDir = p.dirname(p.abspath(__file__))
    cgenffPrm = p.join(stitchingSrcDir, 'par_all36_cgenff.prm')
    cgenffRtf = p.join(stitchingSrcDir, 'top_all36_cgenff.rtf')
    if not p.isfile(cgenffPrm):
        raise FileNotFoundError(f'Cannot find CGenFF PRM file at {cgenffPrm}')
    if not p.isfile(cgenffRtf):
        raise FileNotFoundError(f'Cannot find CGenFF RTF file at {cgenffRtf}')
    config['runtimeInfo']['madeByStitching']['cgenffPrm'] = cgenffPrm
    config['runtimeInfo']['madeByStitching']['cgenffRtf'] = cgenffRtf
    return config