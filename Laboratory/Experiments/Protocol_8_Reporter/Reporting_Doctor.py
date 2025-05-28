from jinja2 import Environment, FileSystemLoader
import pandas as pd
from os import path as p
from shutil import copy
import py3Dmol
from pdbUtils import pdbUtils # Your custom PDB parsing utility
import os
from os import path as p

from . import Reporting_Monster
from . import Reporting_Assistant
from typing import Dict, Any # Added for type hinting

# Placeholder classes (extend if needed)
class FilePath:
    pass
class DirectoryPath:
    pass

def reporter_protocol(config: dict) -> dict:
    """
    Runs the reporter protocol

    Args:
        config (dict): Dictionary containing all information needed

    Returns:
        None
    """
    ## unpack config
    outDir: DirectoryPath = config["pathInfo"]["outputDir"] # type: ignore
    moleculeName: str = config["moleculeInfo"]["moleculeName"]
    ## make directories
    reporterDir: DirectoryPath = p.join(outDir, "08_Parameterisation_Report") # type: ignore
    os.makedirs(reporterDir, exist_ok=True)
    imagesDir: DirectoryPath = p.join(reporterDir, "Images") # type: ignore
    os.makedirs(imagesDir, exist_ok=True)

    ## update config
    config["runtimeInfo"]["madeByReporting"] = {}
    config["runtimeInfo"]["madeByReporting"]["reporterDir"] = reporterDir
    config["runtimeInfo"]["madeByReporting"]["imagesDir"] = imagesDir

    ## run protocols
    timeGanttPng: FilePath    = Reporting_Assistant.generate_gantt_chart(config) # type: ignore
    wriggleData: Dict[str, Any]     = Reporting_Monster.process_wriggle_results(config)
    twistData: Dict[str, Any]       = Reporting_Monster.process_twist_results(config)
    chargesData: Dict[str, Any]     = Reporting_Monster.process_charges_results(config)
    fittingData: Dict[str, Any]     = Reporting_Monster.process_fitting_results(config)

    reportHtml: FilePath = p.join(reporterDir, "drFrankenstein_report.html") # type: ignore
    make_html_report(time_gantt_png=timeGanttPng, 
                     wriggle_data=wriggleData, 
                     twist_data=twistData, 
                     charges_data=chargesData, 
                     fitting_data=fittingData, 
                     molecule_name=moleculeName, 
                     report_html=reportHtml)

    ## update config
    config["checkpointInfo"]["reportingComplete"] = True
    config["runtimeInfo"]["madeByReporting"]["reportHtml"] = reportHtml
    return config


def make_html_report(time_gantt_png: FilePath, 
                     wriggle_data: Dict[str, Any], 
                     twist_data: Dict[str, Any], 
                     charges_data: Dict[str, Any], 
                     fitting_data: Dict[str, Any], 
                     molecule_name: str, 
                     report_html: FilePath) -> None:
    """
    Generates an HTML report from templates and provided data.

    Args:
        time_gantt_png: Path to the Gantt chart PNG for timing information.
        wriggle_data: Dictionary containing data from conformer generation.
        twist_data: Dictionary containing data from torsion scanning.
        charges_data: Dictionary containing data from charge fitting.
        fitting_data: Dictionary containing data from parameter fitting.
        molecule_name: Name of the molecule for the report title.
        report_html: Path to save the generated HTML report.
    """
    template_dir = os.path.join(os.path.dirname(__file__), 'templates')
    if not os.path.exists(template_dir):
        # Fallback if templates are in current working dir (e.g. running from parent)
        template_dir = 'templates'
        if not os.path.exists(template_dir):
            print(f"Error: Template directory '{template_dir}' not found.")
            # Consider raising an error here instead of exit() for better library behavior
            raise FileNotFoundError(f"Template directory '{template_dir}' not found.")

    file_loader = FileSystemLoader(template_dir)
    env = Environment(loader=file_loader)

    # Load the main template
    template = env.get_template('index.html')


    # The {% include %} directive will make these variables available to the included files.
    rendered_html = template.render(
        job_name=molecule_name,
        timeGanttPng=time_gantt_png,
        conformer_data=wriggle_data,
        torsion_data=twist_data,
        charge_data=charges_data,
        fitting_data=fitting_data
    )

    # Save the rendered HTML to a file
    with open(report_html, 'w', encoding='utf-8') as f: # type: ignore
        f.write(rendered_html)
