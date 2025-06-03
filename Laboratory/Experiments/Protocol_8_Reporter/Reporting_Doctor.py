from jinja2 import Environment, FileSystemLoader
import pandas as pd
from os import path as p
from shutil import copy
import py3Dmol
from pdbUtils import pdbUtils
import os
from os import path as p
from . import Reporting_Monster
from . import Reporting_Assistant
from . import Shelly

def reporter_protocol(config: dict) -> None:
    """
    Runs the reporter protocol

    Args:
        config (dict): Dictionary containing all information needed

    Returns:
        None
    """
    outDir = config['pathInfo']['outputDir']
    moleculeName = config['moleculeInfo']['moleculeName']
    reporterDir = p.join(outDir, '08_Parameterisation_Report')
    os.makedirs(reporterDir, exist_ok=True)
    imagesDir = p.join(reporterDir, 'Images')
    os.makedirs(imagesDir, exist_ok=True)
    config['runtimeInfo']['madeByReporting'] = {}
    config['runtimeInfo']['madeByReporting']['reporterDir'] = reporterDir
    config['runtimeInfo']['madeByReporting']['imagesDir'] = imagesDir
    timeGanttPng = Reporting_Assistant.generate_gantt_chart(config)
    wriggleData = Reporting_Monster.process_wriggle_results(config)
    twistData = Reporting_Monster.process_twist_results(config)
    chargesData = Reporting_Monster.process_charges_results(config)
    fittingData = Reporting_Monster.process_fitting_results(config)
    methodsData = Shelly.methods_writer_protocol(config)
    citationsData = Shelly.gather_citations(config)
    reportHtml = p.join(reporterDir, 'drFrankenstein_report.html')
    make_html_report(timeGanttPng, wriggleData, twistData, chargesData, fittingData, moleculeName, methodsData, citationsData, reportHtml)
    config['runtimeInfo']['madeByReporting']['reportHtml'] = reportHtml
    return config

def make_html_report(timeGanttPng, wriggleData, twistData, chargesData, fittingData, moleculeName, methodsData, citationsData, reportHtml):
    templateDir = os.path.join(os.path.dirname(__file__), 'templates')
    if not os.path.exists(templateDir):
        templateDir = 'templates'
        if not os.path.exists(templateDir):
            print(f"Error: Template directory '{templateDir}' not found.")
            exit()
    fileLoader = FileSystemLoader(templateDir)
    env = Environment(loader=fileLoader)
    template = env.get_template('index.html')
    renderedHtml = template.render(job_name=moleculeName, timeGanttPng=timeGanttPng, conformer_data=wriggleData, torsion_data=twistData, charge_data=chargesData, fitting_data=fittingData, methods_data=methodsData, citation_data=citationsData)
    with open(reportHtml, 'w', encoding='utf-8') as f:
        f.write(renderedHtml)