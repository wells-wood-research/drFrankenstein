import os 
from os import path as p
import time

def show_need_cgenff_str(cappedMol2) -> None:
    ## using font ANSI Shadow
    greenText = "\033[32m"
    orangeText = "\033[38;5;172m"
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"

    lightningBar = "🗲 "*int((80//2))

    pausedAscii = """
            ██████╗  █████╗ ██╗   ██╗███████╗███████╗██████╗ 
            ██╔══██╗██╔══██╗██║   ██║██╔════╝██╔════╝██╔══██╗
            ██████╔╝███████║██║   ██║███████╗█████╗  ██║  ██║
            ██╔═══╝ ██╔══██║██║   ██║╚════██║██╔══╝  ██║  ██║
            ██║     ██║  ██║╚██████╔╝███████║███████╗██████╔╝
            ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚══════╝╚══════╝╚═════╝ 
                              """
    
    needStrText = f"CHARMM parameterisation procedure requires a {orangeText}stream (.STR){resetTextColor} file.\n" \
                "As you don't have CGenFF installed locally,\n" \
                "you can  use the CGenFF webserver to generate this:\n" \
                f"1. Go to {greenText}https://cgenff.com/{resetTextColor}\n" \
                "2. Create an account and log in\n" \
                "3. Select \"Start CGenFF Job\"\n" \
                "4. Upload the following MOL2 file:\n"\
                f"{greenText}{cappedMol2}{resetTextColor}\n" \
                "\nOnce CGenFF has finished running:\n" \
                f"Place the {orangeText}stream (.STR){resetTextColor} file in your specified inputDir\n" \
                "and re-run drFrankenstein with the same config file"
    
    getLocalText = f"\n{yellowText}NOTE that if you want to run this protocol in a\n"\
                    "fully automated manner, you can request the CGenFF binaries\n" \
                    f"by emailing: {orangeText}info@silcsbio.com{yellowText}\n"\
                    "NOTE: make sure you CC your PI into this email!\n\n" \
                    "Once you have the binary, supply the location in the\n" \
                    f"{greenText}pathInfo.cgenffExe{yellowText} entry in your config file" 
    
    print(f"{yellowText}{lightningBar}{resetTextColor}")
    print(f"{greenText}{pausedAscii}{resetTextColor}")
    print(f"{yellowText}{lightningBar}{resetTextColor}")

    print(f"{resetTextColor}{needStrText}{resetTextColor}")
    print(f"{getLocalText}")
    print(f"{yellowText}{lightningBar}{resetTextColor}")

def show_config_error(errors: dict) -> None:
    ## using font Bloody
    greenText = "\033[32m"
    redText = "\033[31m"
    orangeText = "\033[38;5;172m"
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"
    redText = "\033[31m"

    lightningBar = "🗲 "*int((94//2))

    configErrorAscii = f"""
 ▄████▄  ▒█████   ███▄    █   █████▒██▓ ▄████    ▓█████  ██▀███   ██▀███   ▒█████   ██▀███  
▒██▀ ▀█ ▒██▒  ██▒ ██ ▀█   █ ▓██   ▒▓██▒██▒ ▀█▒   ▓█   ▀ ▓██ ▒ ██▒▓██ ▒ ██▒▒██▒  ██▒▓██ ▒ ██▒
▒▓█    ▄▒██░  ██▒▓██  ▀█ ██▒▒████ ░▒██▒██░▄▄▄░   ▒███   ▓██ ░▄█ ▒▓██ ░▄█ ▒▒██░  ██▒▓██ ░▄█ ▒
▒▓▓▄ ▄██▒██   ██░▓██▒  ▐▌██▒░▓█▒  ░░██░▓█  ██▓   ▒▓█  ▄ ▒██▀▀█▄  ▒██▀▀█▄  ▒██   ██░▒██▀▀█▄  
▒ ▓███▀ ░ ████▓▒░▒██░   ▓██░░▒█░   ░██░▒▓███▀▒   ░▒████▒░██▓ ▒██▒░██▓ ▒██▒░ ████▓▒░░██▓ ▒██▒
░ ░▒ ▒  ░ ▒░▒░▒░ ░ ▒░   ▒ ▒  ▒ ░   ░▓  ░▒   ▒    ░░ ▒░ ░░ ▒▓ ░▒▓░░ ▒▓ ░▒▓░░ ▒░▒░▒░ ░ ▒▓ ░▒▓░
  ░  ▒    ░ ▒ ▒░ ░ ░░   ░ ▒░ ░      ▒ ░ ░   ░     ░ ░  ░  ░▒ ░ ▒░  ░▒ ░ ▒░  ░ ▒ ▒░   ░▒ ░ ▒░
░       ░ ░ ░ ▒     ░   ░ ░  ░ ░    ▒ ░ ░   ░       ░     ░░   ░   ░░   ░ ░ ░ ░ ▒    ░░   ░ 
░ ░         ░ ░           ░         ░       ░       ░  ░   ░        ░         ░ ░     ░     
░          
"""

    print(f"{yellowText}{lightningBar}{redText}{configErrorAscii}{resetTextColor}\n")
    print(f"{yellowText}\tThe following errors were found in your config file:{resetTextColor}")
    for key, value in errors.items():
        print(f"{yellowText}🗲  {orangeText}{key}:\n\t{redText}{value}{resetTextColor}\n")
    print(f"{yellowText}\tPlease consult the README for more details.{resetTextColor}")
    print(f"{yellowText}{lightningBar}{resetTextColor}")


def show_wriggle_splash() -> None:
    ## using font Delta Corps Priest 1
    greenText = "\033[32m"
    redText = "\033[31m"
    orangeText = "\033[38;5;172m"
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"

    lightningBar = "🗲 "*int((86//2))

    aciiWriggle = """
 ▄█     █▄     ▄████████  ▄█     ▄██████▄     ▄██████▄   ▄█          ▄████████    ▄█ 
███     ███   ███    ███ ███    ███    ███   ███    ███ ███         ███    ███   ███ 
███     ███   ███    ███ ███▌   ███    █▀    ███    █▀  ███         ███    █▀    ███▌
███     ███  ▄███▄▄▄▄██▀ ███▌  ▄███         ▄███        ███        ▄███▄▄▄       ███▌
███     ███ ▀▀███▀▀▀▀▀   ███▌ ▀▀███ ████▄  ▀▀███ ████▄  ███       ▀▀███▀▀▀       ███
███     ███ ▀███████████ ███    ███    ███   ███    ███ ███         ███    █▄    █▀ 
███ ▄█▄ ███   ███    ███ ███    ███    ███   ███    ███ ███▌    ▄   ███    ███    
 ▀███▀███▀    ███    ███ █▀     ████████▀    ████████▀  █████▄▄██   ██████████   ██
              ███    ███                                                         
"""
    print(yellowText + lightningBar +
           resetTextColor + aciiWriggle  +
             yellowText + " "*24 + "GENERATING CONFORMERS WITH GOAT...\n" + 
             lightningBar + resetTextColor)


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def show_twist_splash() -> None:

    ## using font Delta Corps Priest 1
    greenText = "\033[32m"
    redText = "\033[31m"
    orangeText = "\033[38;5;172m"
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"
    tealText = "\033[38;5;37m" 
    brightGreenText = "\033[92m"

    lightningBar = "🗲 "*int((124 // 2))
    asciiTorsionScanning = """

████████▄   ▄██████▄          ███        ▄█    █▄       ▄████████         ███      ▄█     █▄   ▄█     ▄████████     ███     
███   ▀███ ███    ███     ▀█████████▄   ███    ███     ███    ███     ▀█████████▄ ███     ███ ███    ███    ███ ▀█████████▄ 
███    ███ ███    ███        ▀███▀▀██   ███    ███     ███    █▀         ▀███▀▀██ ███     ███ ███▌   ███    █▀     ▀███▀▀██ 
███    ███ ███    ███         ███   ▀  ▄███▄▄▄▄███▄▄  ▄███▄▄▄             ███   ▀ ███     ███ ███▌   ███            ███   ▀ 
███    ███ ███    ███         ███     ▀▀███▀▀▀▀███▀  ▀▀███▀▀▀             ███     ███     ███ ███▌ ▀███████████     ███     
███    ███ ███    ███         ███       ███    ███     ███    █▄          ███     ███     ███ ███           ███     ███     
███   ▄███ ███    ███         ███       ███    ███     ███    ███         ███     ███ ▄█▄ ███ ███     ▄█    ███     ███     
████████▀   ▀██████▀         ▄████▀     ███    █▀      ██████████        ▄████▀    ▀███▀███▀  █▀    ▄████████▀     ▄████▀   
                                                                                                                                                                                                                                                         
"""

    print(yellowText + lightningBar +
           resetTextColor + asciiTorsionScanning  +
             yellowText + " "*42 + "SCANNING TORSION ANGLES...\n" + 
             lightningBar + resetTextColor)



def show_charge_splash() -> None:

    ## using font Delta Corps Priest 1
    greenText = "\033[32m"
    redText = "\033[31m"
    orangeText = "\033[38;5;172m"
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"
    tealText = "\033[38;5;37m" 
    brightGreenText = "\033[92m"

    lightningBar = "🗲 "*int((127 // 2))
    asciiCharge = """
 ▄████████    ▄█    █▄       ▄████████    ▄████████    ▄██████▄     ▄████████      ▄█      ███         ███    █▄     ▄███████▄ 
███    ███   ███    ███     ███    ███   ███    ███   ███    ███   ███    ███     ███  ▀█████████▄     ███    ███   ███    ███ 
███    █▀    ███    ███     ███    ███   ███    ███   ███    █▀    ███    █▀      ███▌    ▀███▀▀██     ███    ███   ███    ███ 
███         ▄███▄▄▄▄███▄▄   ███    ███  ▄███▄▄▄▄██▀  ▄███         ▄███▄▄▄         ███▌     ███   ▀     ███    ███   ███    ███ 
███        ▀▀███▀▀▀▀███▀  ▀███████████ ▀▀███▀▀▀▀▀   ▀▀███ ████▄  ▀▀███▀▀▀         ███▌     ███         ███    ███ ▀█████████▀  
███    █▄    ███    ███     ███    ███ ▀███████████   ███    ███   ███    █▄      ███      ███         ███    ███   ███        
███    ███   ███    ███     ███    ███   ███    ███   ███    ███   ███    ███     ███      ███         ███    ███   ███        
████████▀    ███    █▀      ███    █▀    ███    ███   ████████▀    ██████████     █▀      ▄████▀       ████████▀   ▄████▀      
                                         ███    ███                                                                                                                                                                                                          
"""

    print(yellowText + lightningBar +
           resetTextColor + asciiCharge +
             yellowText + " "*42 + "RUNNING CHARGE CALCULATIONS...\n" + 
             lightningBar + resetTextColor)

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def show_torsion_being_scanned(torsionTag, torsionIndex, nTorsions) -> None:
    greenText = "\033[32m"
    redText = "\033[31m"
    orangeText = "\033[38;5;172m"
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"
    tealText = "\033[38;5;37m"     
    print(f"{yellowText}🗲 🗲{' '*8}SCANNING TORSION:{' '*8}\
[{' '*1}{greenText}{torsionTag}{yellowText}{' '*1}]{' '*8}\
{resetTextColor}({torsionIndex+1}/{nTorsions}){yellowText}{' '*58}🗲 🗲")


def show_stitch_splash() -> None:
    ## using font Delta Corps Priest 1
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"


    lightningBar = "🗲 "*int((102 // 2))
    asciiStitch = """
        ▄████████     ███      ▄█      ███      ▄████████    ▄█    █▄         ▄█      ███             
       ███    ███ ▀█████████▄ ███  ▀█████████▄ ███    ███   ███    ███       ███  ▀█████████▄         
       ███    █▀     ▀███▀▀██ ███▌    ▀███▀▀██ ███    █▀    ███    ███       ███▌    ▀███▀▀██         
       ███            ███   ▀ ███▌     ███   ▀ ███         ▄███▄▄▄▄███▄▄     ███▌     ███   ▀         
     ▀███████████     ███     ███▌     ███     ███        ▀▀███▀▀▀▀███▀      ███▌     ███             
              ███     ███     ███      ███     ███    █▄    ███    ███       ███      ███             
        ▄█    ███     ███     ███      ███     ███    ███   ███    ███       ███      ███             
      ▄████████▀     ▄████▀   █▀      ▄████▀   ████████▀    ███    █▀        █▀      ▄████▀           
                                                                                                      
    ███      ▄██████▄     ▄██████▄     ▄████████     ███        ▄█    █▄       ▄████████    ▄████████ 
▀█████████▄ ███    ███   ███    ███   ███    ███ ▀█████████▄   ███    ███     ███    ███   ███    ███ 
   ▀███▀▀██ ███    ███   ███    █▀    ███    █▀     ▀███▀▀██   ███    ███     ███    █▀    ███    ███ 
    ███   ▀ ███    ███  ▄███         ▄███▄▄▄         ███   ▀  ▄███▄▄▄▄███▄▄  ▄███▄▄▄      ▄███▄▄▄▄██▀ 
    ███     ███    ███ ▀▀███ ████▄  ▀▀███▀▀▀         ███     ▀▀███▀▀▀▀███▀  ▀▀███▀▀▀     ▀▀███▀▀▀▀▀   
    ███     ███    ███   ███    ███   ███    █▄      ███       ███    ███     ███    █▄  ▀███████████ 
    ███     ███    ███   ███    ███   ███    ███     ███       ███    ███     ███    ███   ███    ███ 
   ▄████▀    ▀██████▀    ████████▀    ██████████    ▄████▀     ███    █▀      ██████████   ███    ███ 
                                                                                           ███    ███ 
"""
    print(yellowText + lightningBar +
          resetTextColor + asciiStitch  +
            yellowText + " "*42 + "FITTING PARAMETERS TO QM SCAN DATA...\n" +
            lightningBar + resetTextColor) 


def show_creation_splash() -> None:
    ## using font Delta Corps Priest 1
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"
    lightningBar = "🗲 "*int((76 // 2))
    asciiCreation = """
▀█████████▄     ▄████████  ▄█  ███▄▄▄▄      ▄██████▄       ▄█      ███         
  ███    ███   ███    ███ ███  ███▀▀▀██▄   ███    ███     ███  ▀█████████▄     
  ███    ███   ███    ███ ███▌ ███   ███   ███    █▀      ███▌    ▀███▀▀██     
 ▄███▄▄▄██▀   ▄███▄▄▄▄██▀ ███▌ ███   ███  ▄███            ███▌     ███   ▀     
▀▀███▀▀▀██▄  ▀▀███▀▀▀▀▀   ███▌ ███   ███ ▀▀███ ████▄      ███▌     ███         
  ███    ██▄ ▀███████████ ███  ███   ███   ███    ███     ███      ███         
  ███    ███   ███    ███ ███  ███   ███   ███    ███     ███      ███         
▄█████████▀    ███    ███ █▀    ▀█   █▀    ████████▀      █▀      ▄████▀       
               ███    ███                                                      
    ███      ▄██████▄       ▄█        ▄█     ▄████████    ▄████████            
▀█████████▄ ███    ███     ███       ███    ███    ███   ███    ███            
   ▀███▀▀██ ███    ███     ███       ███▌   ███    █▀    ███    █▀             
    ███   ▀ ███    ███     ███       ███▌  ▄███▄▄▄      ▄███▄▄▄                
    ███     ███    ███     ███       ███▌ ▀▀███▀▀▀     ▀▀███▀▀▀                
    ███     ███    ███     ███       ███    ███          ███    █▄             
    ███     ███    ███     ███▌    ▄ ███    ███          ███    ███            
   ▄████▀    ▀██████▀      █████▄▄██ █▀     ███          ██████████            
                           ▀                                                   
"""

    print(yellowText + lightningBar +
          resetTextColor + asciiCreation  +
            yellowText + " "*10 + "FITTING PARAMETERS TO QM SCAN DATA...\n" +
            lightningBar + resetTextColor) 

import time
import os
import sys

def show_what_have_we_created(moleculeName):
    # Your ASCII art stored in a list of lines
    art = [
        " █     █░ ██░ ██  ▄▄▄     ▄▄▄█████▓    ██░ ██  ▄▄▄    ██▒   █▓▓█████           ",
        "▓█░ █ ░█░▓██░ ██▒▒████▄   ▓  ██▒ ▓▒   ▓██░ ██▒▒████▄ ▓██░   █▒▓█   ▀           ",
        "▒█░ █ ░█ ▒██▀▀██░▒██  ▀█▄ ▒ ▓██░ ▒░   ▒██▀▀██░▒██  ▀█▄▓██  █▒░▒███             ",
        "░█░ █ ░█ ░▓█ ░██ ░██▄▄▄▄██░ ▓██▓ ░    ░▓█ ░██ ░██▄▄▄▄██▒██ █░░▒▓█  ▄           ",
        "░░██▒██▓ ░▓█▒░██▓ ▓█   ▓██▒ ▒██▒ ░    ░▓█▒░██▓ ▓█   ▓██▒▒▀█░  ░▒████▒          ",
        "░ ▓░▒ ▒   ▒ ░░▒░▒ ▒▒   ▓▒█░ ▒ ░░       ▒ ░░▒░▒ ▒▒   ▓▒█░░ ▐░  ░░ ▒░ ░          ",
        "  ▒ ░ ░   ▒ ░▒░ ░  ▒   ▒▒ ░   ░        ▒ ░▒░ ░  ▒   ▒▒ ░░ ░░   ░ ░  ░          ",
        "  ░   ░   ░  ░░ ░  ░   ▒    ░          ░  ░░ ░  ░   ▒     ░░     ░             ",
        "    ░     ░  ░  ░      ░  ░            ░  ░  ░      ░  ░   ░     ░  ░          ",
        "                                                                               ",
        "                            █     █░▓█████                                     ",
        "                           ▓█░ █ ░█░▓█   ▀                                     ",
        "                           ▒█░ █ ░█ ▒███                                       ",
        "                           ░█░ █ ░█ ▒▓█  ▄                                     ",
        "                           ░░██▒██▓ ░▒████▒                                    ",
        "                           ░ ▓░▒ ▒  ░░ ▒░ ░                                    ",
        "                             ▒ ░ ░   ░ ░  ░                                    ",
        "                             ░   ░     ░                                       ",
        "                               ░       ░  ░                                    ",
        " ▄████▄  ██▀███  ▓█████ ▄▄▄     ▄▄▄█████▓▓█████ ▓█████▄    ▄▄██▓     ▄▄██▓     ",
        "▒██▀ ▀█ ▓██ ▒ ██▒▓█   ▀▒████▄   ▓  ██▒ ▓▒▓█   ▀ ▒██▀ ██▌  ▓█████▒   ▓█████▒    ",
        "▒▓█    ▄▓██ ░▄█ ▒▒███  ▒██  ▀█▄ ▒ ▓██░ ▒░▒███   ░██   █▌  ██   ██   ██   ██    ",
        "▒▓▓▄ ▄██▒██▀▀█▄  ▒▓█  ▄░██▄▄▄▄██░ ▓██▓ ░ ▒▓█  ▄ ░▓█▄   ▌  ▒   █▓    ▒   █▓     ",
        "▒ ▓███▀ ░██▓ ▒██▒░▒████▒▓█   ▓██▒ ▒██▒ ░ ░▒████▒░▒████▓   ▒  ██▒    ▒  ██▒     ",
        "░ ░▒ ▒  ░ ▒▓ ░▒▓░░░ ▒░ ░▒▒   ▓▒█░ ▒ ░░   ░░ ▒░ ░ ▒▒▓  ▒     ██ ▒      ██ ▒     ",
        "  ░  ▒    ░▒ ░ ▒░ ░ ░  ░ ▒   ▒▒ ░   ░     ░ ░  ░ ░ ▒  ▒        ▒         ▒     ",
        "░         ░░   ░    ░    ░   ▒    ░         ░    ░ ░  ░     ██▓       ██▓      ",
        "░ ░        ░        ░  ░     ░  ░           ░  ░   ░        ▒         ▒        ",
        "░                                                ░          ▒         ▒        "
    ]
    
    # Clear terminal function that works cross-platform
    def clear_terminal():
        os.system('cls' if os.name == 'nt' else 'clear')
    
    # Animation parameters
    duration = 3  # seconds
    frames_per_second = 60
    total_frames = duration * frames_per_second
    max_shift = 2  # maximum spaces to shift left/right
    green_text = "\033[32m"  # ANSI escape code for green text
    reset_color = "\033[0m"  # Reset color to default
    
    # Animation loop
    for frame in range(total_frames):
        clear_terminal()
        # Calculate shift amount using sine wave for smooth motion
        import math
        shift = int(max_shift * math.sin(frame * 0.5))
        
        # Print each line with appropriate spacing and green color
        for line in art:
            padding = " " * (max_shift + shift)
            print(f"{padding}{green_text}{line}{reset_color}")
        
        # Small delay between frames
        time.sleep(1/frames_per_second)
    
    # Display final centered version
    clear_terminal()
    for line in art:
        print(f"{' ' * max_shift}{green_text}{line}{reset_color}")

    yellowText = "\033[33m"
    resetTextColor = "\033[0m"
    print(" "*16+ yellowText + "AMBER PARAMETERS FOR YOUR MOLECULE:" + " "*4 + green_text + moleculeName + resetTextColor)





#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def print_botched(errorReport) -> None:
    greenText = "\033[32m"
    redText = "\033[31m"
    orangeText = "\033[38;5;172m"
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"
    tealColor = "\033[38;5;37m" 

    # run(["clear"])
    print(redText+
          f"""
🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

                ▄▄▄▄    ▒█████  ▄▄▄█████▓ ▄████▄   ██░ ██ ▓█████ ▓█████▄  ▐██▌ 
                ▓█████▄ ▒██▒  ██▒▓  ██▒ ▓▒▒██▀ ▀█  ▓██░ ██▒▓█   ▀ ▒██▀ ██▌ ▐██▌ 
                ▒██▒ ▄██▒██░  ██▒▒ ▓██░ ▒░▒▓█    ▄ ▒██▀▀██░▒███   ░██   █▌ ▐██▌ 
                ▒██░█▀  ▒██   ██░░ ▓██▓ ░ ▒▓▓▄ ▄██▒░▓█ ░██ ▒▓█  ▄ ░▓█▄   ▌ ▓██▒ 
                ░▓█  ▀█▓░ ████▓▒░  ▒██▒ ░ ▒ ▓███▀ ░░▓█▒░██▓░▒████▒░▒████▓  ▒▄▄  
                ░▒▓███▀▒░ ▒░▒░▒░   ▒ ░░   ░ ░▒ ▒  ░ ▒ ░░▒░▒░░ ▒░ ░ ▒▒▓  ▒  ░▀▀▒ 
                ▒░▒   ░   ░ ▒ ▒░     ░      ░  ▒    ▒ ░▒░ ░ ░ ░  ░ ░ ▒  ▒  ░  ░ 
                ░    ░ ░ ░ ░ ▒    ░      ░         ░  ░░ ░   ░    ░ ░  ░     ░ 
                ░          ░ ░           ░ ░       ░  ░  ░   ░  ░   ░     ░    
                    ░                   ░                        ░          
🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
          """
          +resetTextColor)
    

    if errorReport is not None:
        print(f"{redText}{' '*7}🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲{resetTextColor}")
        print(f"{' '*7}For System:\t\t{yellowText}{errorReport['pdbName']}{resetTextColor}")
        print("")
        print(f"{tealColor}{' '*7}{'#'*4}{' '*7}Traceback{' '*7}{'#'*4}{resetTextColor}")
        print(f"{' '*7}In Script:\t\t{orangeText}{errorReport['scriptName']}{resetTextColor}")
        print(f"{' '*7}In Function:\t\t{orangeText}{errorReport['functionName']}{resetTextColor}")


        print(f"{' '*7}With Error:\t\t{redText}{errorReport['errorType']}{resetTextColor}")
        print(f"{' '*7}With Message:\t\t{redText}{errorReport['errorMessage']}{resetTextColor}")
        print(f"{' '*7}At Line Number:\t\t{redText}{errorReport['lineNumber']}{resetTextColor}")

        print(f"{tealColor}{' '*7}{'#'*4}{' '*7}Full Debug Traceback{' '*7}{'#'*4}{resetTextColor}")

        print(f"\t{orangeText}{'LINE NUMBER':<10}{yellowText}{'FUNCTION':>30}{resetTextColor}\t/path/to/crashed/{tealColor}script_name.py{resetTextColor}")
        print(f"\t{'---':<10}{'---':>30}\t{'---'}")
        for tracebackLine in errorReport["fullTraceBack"]:
            scriptPath = tracebackLine.split(":")[0]
            scriptDir = p.dirname(scriptPath)
            scriptName = p.basename(scriptPath)
            lineNumber = tracebackLine.split(":")[1].split("in")[0].strip()
            functionName = tracebackLine.split(":")[1].split("in")[1].strip()
            print(f"\t{orangeText}{lineNumber:<10}{yellowText}{functionName:>30}{resetTextColor}\t{scriptDir}/{tealColor}{scriptName}{resetTextColor}")
    print(f"{redText}🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲")
    print(resetTextColor)


if __name__ == "__main__":
    show_wriggle_splash()
    show_twist_splash()
    show_charge_splash()
    show_torsion_being_scanned("N-C-CA-CB", 1, 10)
    show_stitch_splash()
    show_creation_splash()
    show_what_have_we_created("ALANINE")