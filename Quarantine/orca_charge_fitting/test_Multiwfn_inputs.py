import subprocess
import os
import pexpect
from time import sleep 
def run_external_program(command,
                         chargeGroupRestraintFile, conformerFile, 
                         output_file):
    child = pexpect.spawn(command, encoding='utf-8')

    # Capture all output until the expected prompt
    output = ""

    ####### GET TO RESP MAIN MENU #######

    ## look for main menu, input "7" for charge calculations
    child.expect(r'.*300.*\n?\r?')
    child.sendline("7")

    ## look for charge calculation menu, input "18" for RESP calculations
    child.expect(r'.*20.*\n?\r?')
    child.sendline("18")


    ###### LOAD CONFORMERS #######

    ## look for RESP main menu, input "-1" to load conformer file
    child.expect(r'.*11.*\n?\r?')
    child.sendline("-1")

    ## look for directory input prompt, input directory
    child.expect(r'.*Input.*\n?\r?')
    child.sendline(conformerFile)


    ####### LOAD CHARGE CONSTRAINTS #######
    ## look for RESP main menu, input "6" to load charge constraint file
    child.expect(r'.*11.*\n?\r?')
    child.sendline("6")


    ## look for conformation, input "1" to proceed
    child.expect(r'.*1.*\n?\r?')
    child.sendline("1")


    ## look for directory input prompt, input directory
    child.expect(r'.*Input.*\n?\r?')
    child.sendline(chargeGroupRestraintFile)

    ####### RUN CALCULATIONS #######

    ## look for RESP main menu, input "2" to run charge calculations
    child.expect(r'.*11.*\n?\r?')
    child.sendline("2")
    output += child.after
    print(output)
    ## look for Note at the end of calculation, input "q" to quit
    child.expect(r'.*Note: Because.*\n?\r?')
    child.sendline("q")
    output += child.after
    print(output)

    # Write the output to the file
    with open(output_file, 'w') as file:
        file.write(output)





multiWfnDir = "/home/esp/bin/Multiwfn_3.8_dev_bin_Linux_noGUI"
os.chdir(multiWfnDir)
multiWfnExe = "/home/esp/bin/Multiwfn_3.8_dev_bin_Linux_noGUI/Multiwfn_noGUI"
conformerFile = "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs/charge_calculations/charge_fitting/conformer_list.txt"
chargeGroupRestraintFile = "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs/charge_calculations/charge_fitting/charge_group_constraints.txt"

inputMolden = "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs/charge_calculations/charge_fitting/ALA_capped_conformer_1.molden.input"
outFile = "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs/charge_calculations/charge_fitting/MultiWFN_output.txt"
command = f"{multiWfnExe} {inputMolden}"

run_external_program(command, conformerFile=conformerFile, chargeGroupRestraintFile=chargeGroupRestraintFile, output_file=outFile)

