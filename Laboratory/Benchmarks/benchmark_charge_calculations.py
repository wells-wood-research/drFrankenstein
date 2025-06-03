import os
from os import path as p
import yaml
import time
from subprocess import call
blankConfig = {'moleculeInfo': {'charge': 0, 'multiplicity': 1, 'moleculePdb': '/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs/capped_amino_acids/ALA_capped.pdb'}, 'pathInfo': {'orcaExe': '/home/esp/bin/orca_6_0_1_linux_x86-64_shared_openmpi416/orca', 'multiWfnDir': '/home/esp/bin/Multiwfn_3.8_dev_bin_Linux_noGUI/'}, 'chargeFittingInfo': {'nCoresPerCalculation': 8, 'nConformers': 4}}
methodsToTry = {'HF-HF': {'optMethod': '! HF 6-31G(d)', 'singlePointMethod': '! HF 6-31G(d)'}, 'XTB2-HF': {'optMethod': '! XTB2 ALPB(WATER)', 'singlePointMethod': '! HF 6-31G(d)'}, 'XTB2-revPBE-SVP': {'optMethod': '! XTB2 ALPB(WATER)', 'singlePointMethod': '! revPBE def2-SVP D3BJ CPCM(water)'}, 'XTB2-revPBE-TZVP': {'optMethod': '! XTB2 ALPB(WATER)', 'singlePointMethod': '! revPBE def2-SVP D3BJ CPCM(water)'}}
benchmarkDir = '/home/esp/scriptDevelopment/drFrankenstein/Benchmarking'
methodYamls = []
for (methodName, methodOptions) in methodsToTry.items():
    methodConfig = blankConfig.copy()
    methodConfig['chargeFittingInfo'].update({'optMethod': methodOptions['optMethod'], 'singlePointMethod': methodOptions['singlePointMethod']})
    methodConfig['pathInfo']['outDir'] = p.join(benchmarkDir, methodName)
    methodYaml = p.join(benchmarkDir, f'{methodName}.yaml')
    with open(methodYaml, 'w') as f:
        yaml.dump(methodConfig, f, default_flow_style=False)
    methodYamls.append(methodYaml)
timeReport = p.join(benchmarkDir, 'time_report.txt')
for methodYaml in methodYamls:
    print(methodYaml)
    startTime = time.time()
    call(['python', '/home/esp/scriptDevelopment/drFrankenstein/src/drOrcaCharges.py', '--config', methodYaml])
    endTime = time.time()
    elapsedTime = endTime - startTime
    hours = int(elapsedTime // 3600)
    minutes = int(elapsedTime % 3600 // 60)
    seconds = int(elapsedTime % 60)
    formattedTime = f'{hours:02}:{minutes:02}:{seconds:02}'
    with open(timeReport, 'a+') as f:
        f.write(f'{methodYaml} {formattedTime}\n')