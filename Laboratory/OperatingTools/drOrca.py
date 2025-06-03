import os
from os import path as p
from subprocess import call

class FilePath:
    pass

class DirectoryPath:
    pass
from OperatingTools import file_parsers

def run_orca(orcaInput, orcaOutput, config):
    orcaExe = config['pathInfo']['orcaExe']
    orcaCommand = [orcaExe, orcaInput]
    with open(orcaOutput, 'w') as outputFile:
        try:
            call(orcaCommand, stdout=outputFile, stderr=outputFile)
        except Exception as e:
            raise e

def make_orca_input_qmmm_opt(inputXyz: FilePath, outDir: DirectoryPath, moleculeInfo: dict, qmMethod: str, qmAtoms: str, parameterFile: FilePath, nCpus: int) -> FilePath:
    charge = moleculeInfo['charge']
    multiplicity = moleculeInfo['multiplicity']
    qmmmOptOrcaInput = p.join(outDir, f'QMMM_orca_opt.inp')
    with open(qmmmOptOrcaInput, 'w') as f:
        f.write(' # --------------------------------- #\n')
        f.write(' #  QM/MM Optimisation               #\n')
        f.write(' # --------------------------------- #\n')
        f.write(f'! QMMM {qmMethod} Opt\n')
        f.write(f'%pal nprocs {nCpus}\nend\n')
        f.write('%qmmm\n')
        f.write(f"{' ' * 4}QMAtoms {qmAtoms} end\n")
        f.write(f'''{' ' * 4}ORCAFFFilename "{parameterFile}"\n''')
        f.write('Rigid_MM_Water TRUE\n')
        f.write('end\n')
        f.write(f'*xyzfile {charge} {multiplicity} {inputXyz}\n')
    return qmmmOptOrcaInput

def make_orca_input_qmmm_singlepoint(inputXyz: FilePath, outDir: DirectoryPath, moleculeInfo: dict, qmMethod: str, qmAtoms: str, parameterFile: FilePath) -> FilePath:
    charge = moleculeInfo['charge']
    multiplicity = moleculeInfo['multiplicity']
    qmmmSinglepointOrcaInput = p.join(outDir, f'QMMM_orca_sp.inp')
    with open(qmmmSinglepointOrcaInput, 'w') as f:
        f.write(' # --------------------------------- #\n')
        f.write(' #  QM/MM Single-Point Energy        #\n')
        f.write(' # --------------------------------- #\n')
        f.write(f'! QMMM {qmMethod} SP\n')
        f.write('%qmmm\n')
        f.write(f"{' ' * 4}QMAtoms {qmAtoms} end\n")
        f.write(f'''{' ' * 4}ORCAFFFilename "{parameterFile}"\n''')
        f.write('end\n')
        f.write(f'*xyzfile {charge} {multiplicity} {inputXyz}\n')
    return qmmmSinglepointOrcaInput

def make_orca_input_for_solvator(inputXyz: FilePath, outDir: DirectoryPath, moleculeInfo: dict, qmMethod: str, solvationMethod: str, nWaters: int) -> FilePath:
    """
    Creates and ORCA input file to run the SOLVATOR protocol
    """
    charge = moleculeInfo['charge']
    multiplicity = moleculeInfo['multiplicity']
    solvatorOrcaInput = p.join(outDir, f'SOLVATOR_orca.inp')
    with open(solvatorOrcaInput, 'w') as f:
        f.write(' # --------------------------------- #\n')
        f.write(' #  SOLVATOR solvent addition        #\n')
        f.write(' # --------------------------------- #\n')
        f.write(f'!{qmMethod} {solvationMethod}\n')
        f.write(f'%SOLVATOR\n')
        f.write(f"{' ' * 4}NSOLV {nWaters}\n")
        f.write(f"{' ' * 4}CLUSTERMODE STOCHASTIC\n")
        f.write('END\n')
        f.write(f'*xyzfile {charge} {multiplicity} {inputXyz}\n')
    return solvatorOrcaInput

def write_goat_input(conformerDir, cappedXyz, config):
    charge = config['moleculeInfo']['charge']
    multiplicity = config['moleculeInfo']['multiplicity']
    goatOrcaInput = p.join(conformerDir, f'GOAT_orca.inp')
    with open(goatOrcaInput, 'w') as f:
        f.write(' # --------------------------------- #\n')
        f.write(' #  GOAT conformation generation     #\n')
        f.write(' # --------------------------------- #\n')
        f.write('!XTB2 GOAT\n')
        f.write('%PAL NPROCS 16 END\n')
        f.write('%GOAT\n')
        f.write('\tMAXITERMULT 1\n')
        f.write('\tFREEZEAMIDES TRUE\n')
        f.write('\tFREEZECISTRANS TRUE\n ')
        f.write('END\n')
        f.write(f'*xyzfile {charge} {multiplicity} {cappedXyz}\n')
    return goatOrcaInput

def make_orca_input_for_opt(inputXyz, outDir, moleculeInfo, qmMethod, solvationMethod):
    charge = moleculeInfo['charge']
    multiplicity = moleculeInfo['multiplicity']
    orcaInputFile = p.join(outDir, 'orca_opt.inp')
    with open(orcaInputFile, 'w') as f:
        f.write(' # --------------------------------- #\n')
        f.write(' #  Geometry Optimisation            #\n')
        f.write(' # --------------------------------- #\n')
        if solvationMethod is None:
            f.write(f'! {qmMethod} Opt\n')
        else:
            f.write(f'! {qmMethod} {solvationMethod} Opt\n')
        f.write(f'*xyzfile {charge} {multiplicity} {inputXyz}\n\n')
    return orcaInputFile

def make_orca_input_for_scan(inputXyz, outDir, moleculeInfo, qmMethod, solvationMethod, torsionIndexes, scanAngles):
    charge = moleculeInfo['charge']
    multiplicity = moleculeInfo['multiplicity']
    orcaInputFile = p.join(outDir, 'orca_scan.inp')
    with open(orcaInputFile, 'w') as f:
        f.write(' # --------------------------------- #\n')
        f.write(' #  Torsion Scan                     #\n')
        f.write(' # --------------------------------- #\n')
        if solvationMethod is None:
            f.write(f'! {qmMethod} Opt\n')
        else:
            f.write(f'! {qmMethod} {solvationMethod} Opt\n')
        if not scanAngles is None:
            f.write('%geom Scan\n')
            torsionText: str = f"D {' '.join(map(str, torsionIndexes))} = {scanAngles}\n"
            f.write(torsionText)
            f.write('end\n')
            f.write('end\n')
        f.write(f'*xyzfile {charge} {multiplicity} {inputXyz}\n\n')
    return orcaInputFile

def make_orca_input_for_singlepoint(inputXyz, outDir, moleculeInfo, qmMethod, solvationMethod):
    charge = moleculeInfo['charge']
    multiplicity = moleculeInfo['multiplicity']
    orcaInputFile = p.join(outDir, 'orca_sp.inp')
    with open(orcaInputFile, 'w') as f:
        f.write(' # --------------------------------- #\n')
        f.write(' #  Single Point Calculation         #\n')
        f.write(' # --------------------------------- #\n')
        if solvationMethod is None:
            f.write(f'! {qmMethod} SP\n')
        else:
            f.write(f'! {qmMethod} {solvationMethod} SP\n')
        f.write(f'*xyzfile {charge} {multiplicity} {inputXyz}\n\n')
    return orcaInputFile

def did_orca_finish_normallly(orcaOut):
    with open(orcaOut, 'r') as f:
        lines = f.readlines()
        for line in reversed(lines):
            if '****ORCA TERMINATED NORMALLY****' in line:
                return True
            elif 'ORCA finished with an error in the energy calculation' in line:
                return False
    return False