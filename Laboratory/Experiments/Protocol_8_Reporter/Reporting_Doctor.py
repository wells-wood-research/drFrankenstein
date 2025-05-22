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

def reporter_protocol(config: dict) -> None:
    """
    Runs the reporter protocol

    Args:
        config (dict): Dictionary containing all information needed

    Returns:
        None
    """
    ## unpack config
    outDir = config["pathInfo"]["outputDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    ## make directories
    reporterDir = p.join(outDir, "08_Parameterisation_Report")
    os.makedirs(reporterDir, exist_ok=True)
    imagesDir = p.join(reporterDir, "Images")
    os.makedirs(imagesDir, exist_ok=True)

    ## update config
    config["runtimeInfo"]["madeByReporting"] = {}
    config["runtimeInfo"]["madeByReporting"]["reporterDir"] = reporterDir
    config["runtimeInfo"]["madeByReporting"]["imagesDir"] = imagesDir

    ## run protocols
    timeGanttPng    = Reporting_Assistant.generate_gantt_chart(config)
    wriggleData     = Reporting_Monster.process_wriggle_results(config)
    twistData       = Reporting_Monster.process_twist_results(config)
    chargesData     = Reporting_Monster.process_charges_results(config)
    fittingData     = Reporting_Monster.process_fitting_results(config)

    reportHtml = p.join(reporterDir, "drFrankenstein_report.html")
    make_html_report(timeGanttPng, wriggleData, twistData, chargesData, fittingData, moleculeName, reportHtml)

    config["checkpointInfo"]["reportingComplete"] = True

    return config


def make_html_report(timeGanttPng, wriggleData, twistData, chargesData, fittingData, moleculeName, reportHtml):

    template_dir = os.path.join(os.path.dirname(__file__), 'templates')
    if not os.path.exists(template_dir):
        # Fallback if templates are in current working dir (e.g. running from parent)
        template_dir = 'templates'
        if not os.path.exists(template_dir):
            print(f"Error: Template directory '{template_dir}' not found.")
            exit()

    file_loader = FileSystemLoader(template_dir)
    env = Environment(loader=file_loader)

    # Load the main template
    template = env.get_template('index.html')


    # The {% include %} directive will make these variables available to the included files.
    rendered_html = template.render(
        job_name=moleculeName,
        timeGanttPng=timeGanttPng,
        conformer_data=wriggleData,
        torsion_data=twistData,
        charge_data=chargesData,
        fitting_data=fittingData
    )

    # Save the rendered HTML to a file
    with open(reportHtml, 'w', encoding='utf-8') as f:
        f.write(rendered_html)

    print(f"Generated HTML report: {os.path.abspath(reportHtml)}")