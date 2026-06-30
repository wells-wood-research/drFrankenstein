import os
from os import path as p
import re
import textwrap
from functools import lru_cache
from textwrap import dedent
import ruamel.yaml as ruamel
from ruamel.yaml.comments import CommentedMap
class FilePath(str):
    pass


METHOD_REFERENCE_FILE = p.join(p.dirname(__file__), "methods_to_references.yaml")


def read_input_yaml(configFile: FilePath) -> dict:
    """Read a YAML file into a Python dictionary."""
    yellow = "\033[33m"
    reset = "\033[0m"
    teal = "\033[38;5;37m"
    try:
        ruamelParser = ruamel.YAML()
        ruamelParser.preserve_quotes = True  # Ensure quotes are preserved
        with open(configFile, "r") as yamlFile:
            config: dict = ruamelParser.load(yamlFile)
            return config
        

    except FileNotFoundError:
        print(f"-->{' '*4}Config file {configFile} not found.")
        exit(1)
    except ruamel.YAMLError as exc:
        print(f"-->{' '*4}{yellow}Error while parsing YAML file:{reset}")
        if hasattr(exc, 'problem_mark'):
            mark = exc.problem_mark
            print(f"{' '*7}Problem found at line {mark.line + 1}, column {mark.column + 1}:")
            if exc.context:
                print(f"{' '*7}{exc.problem} {exc.context}")
            else:
                print(f"{' '*7}{exc.problem}")
            print(f"{' '*7}Please correct the data and retry.")
        else:
            print(f"{' '*7}Something went wrong while parsing the YAML file.")
            print(f"{' '*7}Please check the file for syntax errors.")
        print(f"\n{teal}TIP:{reset} Large language models (LLMs) like GPT-4 can be helpful for debugging YAML files.")
        print(f"{' '*5}If you get stuck with the formatting, ask a LLM for help!")
        exit(1)


@lru_cache(maxsize=1)
def load_method_reference_map() -> dict:
    """Load method-to-reference mappings for report citation augmentation."""
    if not p.isfile(METHOD_REFERENCE_FILE):
        return {}

    yaml_parser = ruamel.YAML(typ="safe")
    with open(METHOD_REFERENCE_FILE, "r") as yaml_file:
        data = yaml_parser.load(yaml_file) or {}
    return data.get("methods_to_references", {})


def normalize_method_string(method_text: str) -> str:
    """Normalize a method string to improve fuzzy key matching."""
    if not method_text:
        return ""
    normalized = str(method_text).strip().lower().replace("²", "2").replace("ω", "w")
    # Preserve internal method punctuation and only normalize wrapper punctuation.
    normalized = normalized.strip("!%#.,;:()[]{}<>")
    return normalized


def tokenize_method_input(method_text: str) -> list[str]:
    """Split a method input line and return normalized tokens for exact matching."""
    if not method_text:
        return []
    raw_tokens = str(method_text).split()
    token_candidates = set()

    def _add_candidate(candidate: str) -> None:
        normalized = normalize_method_string(candidate)
        if normalized:
            token_candidates.add(normalized)

    split_separators = ["=", ",", "/", "(", ")"]
    for token in raw_tokens:
        # Recursively split tokens on common ORCA delimiters.
        queue = [token]
        visited = set()
        while queue:
            current = queue.pop()
            if current in visited:
                continue
            visited.add(current)
            _add_candidate(current)

            for separator in split_separators:
                if separator in current:
                    for part in current.split(separator):
                        if part and part not in visited:
                            queue.append(part)

            # Handle composite dispersion forms, e.g. PBE0-D3BJ, B97-D4
            if "-" in current and re.search(r"(d2|d3|d4|d5|vv10|nl)", current, re.IGNORECASE):
                for part in current.split("-"):
                    if part and part not in visited:
                        queue.append(part)

    return sorted(token_candidates)


def find_method_reference_hits(method_text: str, reference_map: dict) -> list[tuple[str, dict]]:
    """Return matched (method_key, method_data) pairs found in method_text.

    Uses regex-style full-line searching and keeps the most specific overlaps.
    """
    if not method_text:
        return []

    searchable_text = normalize_method_string(method_text)
    token_set = set(tokenize_method_input(method_text))
    if not searchable_text:
        return []

    # Expand common aliases/synonyms to canonical method keys.
    alias_map = {
        "xtb0": "gfn0-xtb",
        "xtb1": "gfn-xTB".lower(),
        "xtb2": "gfn2-xtb",
        "xtbff": "gfn-ff",
        "cpcm": "cpcm-x",
        "cpcmx": "cpcm-x",
        "ddcosmo": "ddcosmo",
    }
    expanded_tokens = set(token_set)
    for token in token_set:
        if token in alias_map:
            expanded_tokens.add(alias_map[token])

    # Make alias-expanded forms searchable in regex matching.
    if expanded_tokens:
        searchable_text = searchable_text + " " + " ".join(sorted(expanded_tokens))

    candidates: list[tuple[int, int, int, str, dict]] = []
    for method_key, method_data in reference_map.items():
        norm_key = normalize_method_string(method_key)
        if not norm_key:
            continue

        # Regex-style full-line search.
        match = re.search(re.escape(norm_key), searchable_text)
        if not match:
            continue
        start, end = match.span()
        candidates.append((-(end - start), start, end, method_key, method_data))

    # Choose the option that matches the most (longest), suppressing shorter overlaps.
    candidates.sort(key=lambda item: (item[0], item[1], item[3].lower()))
    selected: list[tuple[int, int, str, dict]] = []
    selected_spans: list[tuple[int, int]] = []
    for neg_len, start, end, method_key, method_data in candidates:
        overlaps = any(not (end <= span_start or start >= span_end) for span_start, span_end in selected_spans)
        if overlaps:
            continue
        selected.append((start, end, method_key, method_data))
        selected_spans.append((start, end))

    selected.sort(key=lambda item: (item[0], item[2].lower()))
    return [(method_key, method_data) for _, _, method_key, method_data in selected]


def extract_doi_from_citation(citation_text: str) -> str | None:
    """Extract DOI string and normalize to doi.org/<doi> when possible."""
    if not citation_text:
        return None

    url_match = re.search(r"(doi\.org/10\.\d{4,9}/\S+)", citation_text, flags=re.IGNORECASE)
    if url_match:
        doi = url_match.group(1).strip().rstrip(".,;")
        doi = re.sub(r"^https?://", "", doi, flags=re.IGNORECASE)
        return doi

    doi_match = re.search(r"(10\.\d{4,9}/\S+)", citation_text, flags=re.IGNORECASE)
    if doi_match:
        raw = doi_match.group(1).strip().rstrip(".,;")
        return f"doi.org/{raw}"

    return None


def extract_authors_from_citation(citation_text: str) -> str:
    """Best-effort extraction of the author block from a compact citation string."""
    if not citation_text:
        return ""

    # Remove DOI suffix for cleaner author extraction.
    citation_no_doi = re.sub(r"(DOI:\\s*\\S+|doi\\.org/\\S+)", "", citation_text, flags=re.IGNORECASE).strip()

    # Usually authors occupy the first sentence before the title.
    head = citation_no_doi.split(". ", 1)[0].strip()
    return head


def extract_title_from_citation(citation_text: str) -> str:
    """Best-effort extraction of paper title from compact citation text."""
    if not citation_text:
        return ""

    cleaned = re.sub(r"(DOI:\s*\S+|doi\.org/\S+)", "", citation_text, flags=re.IGNORECASE).strip()
    parts = [p.strip() for p in cleaned.split(". ") if p.strip()]
    if len(parts) >= 2:
        # Typically: Authors. Title. Journal...
        candidate = parts[1].rstrip(".")
        # If this looks like a journal fragment (no explicit title), use full citation text.
        if len(candidate.split()) <= 3 or re.search(r"\b(Phys|Chem|Rev|J\.|Comput|Lett)\b", candidate):
            return cleaned.rstrip(".")
        return candidate
    # Fallback to full citation text if title cannot be cleanly isolated.
    return cleaned.rstrip(".")


def build_qm_methods_citation_entries(method_fields: list[tuple[str, str]]) -> list[dict]:
    """Build structured QM citation entries for report rendering."""
    reference_map = load_method_reference_map()
    if not reference_map:
        return []

    entries: list[dict] = []
    seen_method_keys = set()

    for _, field_value in method_fields:
        for method_key, method_data in find_method_reference_hits(field_value, reference_map):
            if method_key in seen_method_keys:
                continue
            seen_method_keys.add(method_key)
            refs = method_data.get("references", [])
            if not refs:
                continue

            best_citation = ""
            best_doi = None
            for ref in refs:
                citation = (ref.get("citation") or "").strip()
                if not citation:
                    continue
                # Use first available reference as primary citation to keep
                # citation ordering stable with reference_ids ordering.
                if not best_citation:
                    best_citation = citation
                    best_doi = extract_doi_from_citation(citation)

            if not best_citation:
                continue

            authors = extract_authors_from_citation(best_citation)
            title = extract_title_from_citation(best_citation)
            # Reuse same style as get_last_authors: split on ';' and take the final author.
            authors_list = authors.split(";") if authors else []
            last_author = authors_list[-1].strip() if authors_list else ""

            entry = {
                "method": method_key,
                "authors": authors,
                "title": title,
                "last_author": last_author,
                "doi": best_doi,
                "citation": best_citation
            }
            entries.append(entry)

    return entries
############################################################################
def methods_writer_protocol(config: dict) -> dict:
    """Build the report methods text sections."""
    conformerMethods : str = write_conformer_method(config)


    scanningMethods: str = write_torsion_scanning_method(config)


    chargeMethods: str = write_charge_calculation_method(config)

    fittingMethods = write_fitting_methods(config)

    methods = {
        "conformerMethods": conformerMethods,
        "scanningMethods": scanningMethods,
        "chargeMethods": chargeMethods,
        "fittingMethods": fittingMethods,
        "qmMethodsCitations": {
            "conformer": build_qm_methods_citation_entries([
                ("conformerMethod", "GFN2-XTB"),
            ]),
            "scanning": build_qm_methods_citation_entries([
                ("scanMethod", config["torsionScanInfo"]["scanMethod"]),
                ("singlePointMethod (scanning)", config["torsionScanInfo"]["singlePointMethod"]),
            ]),
            "charge": build_qm_methods_citation_entries([
                ("optMethod (charge calculations)", config["chargeFittingInfo"]["optMethod"]),
                ("singlePointMethod (charge calculations)", config["chargeFittingInfo"]["singlePointMethod"]),
            ]),
            "fitting": build_qm_methods_citation_entries([]),
        }
    }

    return methods
############################################################################


def write_conformer_method(config: dict) -> str:
    """Write the conformer-generation methods text."""
    ## unpack config ##
    moleculeName: str = config["moleculeInfo"]["moleculeName"]
    backboneAliases: dict = config["moleculeInfo"]["backboneAliases"]
    nConformersGenerated = config["runtimeInfo"]["madeByConformers"]["nConformersGenerated"]
    isCapped = False
    if len(backboneAliases["N"]) != 0:
        isCapped = True

    
    conformerLabel = f"capped {moleculeName}" if isCapped else moleculeName
    conformerMethod = dedent(f"""
A library of {nConformersGenerated} low-energy conformers of {conformerLabel} was generated using ORCA's GOAT protocol. \
GOAT conformer generation was run with the GFN2-XTB method. \
These conformers were then used as starting structures for torsion scanning and charge fitting.
                            """)
    
    return conformerMethod

############################################################################
def write_torsion_scanning_method(config: dict) -> str:
    """Write the torsion-scanning methods text."""
    ## unpack config ##
    moleculeName: str = config["moleculeInfo"]["moleculeName"]
    nConformersGenerated = config["runtimeInfo"]["madeByConformers"]["nConformersGenerated"]
    nConformers = config["torsionScanInfo"]["nConformers"]
    scanMethod = config["torsionScanInfo"]["scanMethod"]
    scanSolvationMethod = config["torsionScanInfo"]["scanSolvationMethod"]
    singlePointMethod = config["torsionScanInfo"]["singlePointMethod"]
    singlePointSolvationMethod = config["torsionScanInfo"]["singlePointSolvationMethod"]

    if nConformers == -1:
        nConformersUsed = nConformersGenerated
    elif nConformers > nConformersGenerated:
        nConformersUsed = nConformersGenerated
    else:
        nConformersUsed = nConformers

    torsionMethod = dedent(f"""
Torsion scans were performed on the set of rotatable dihedrals selected by the torsion-scan settings for {moleculeName}. \
Each selected torsion was scanned over 360 degrees in both forward and backward directions using 10 degree increments (37 points per direction). \
Relaxed scans were run at the {scanMethod} level{" with the " + scanSolvationMethod + " solvation model" if scanSolvationMethod is not None else ""}. \
                            """)
    
    if singlePointMethod:
        singlePointText = dedent(f"""
The resulting scan geometries were refined with single-point calculations at the {singlePointMethod} level \
{"with the " + singlePointSolvationMethod + " solvation model" if singlePointSolvationMethod is not None else ""}. \
                                 """)
        torsionMethod = torsionMethod + singlePointText
    
    averagingText = dedent(f"""
The torsion-scanning workflow was repeated for {nConformersUsed} starting conformers. \
Only successfully completed forward/backward scan pairs were retained, and the final profile for each torsion was taken as the average over those retained scans. \
                            """)

    torsionMethod = torsionMethod + averagingText

    return torsionMethod
############################################################################
def write_charge_calculation_method(config: dict) -> str:
    """Write the charge-calculation methods text."""
    ## unpack config ##
    moleculeName: str = config["moleculeInfo"]["moleculeName"]
    nConformersGenerated = config["runtimeInfo"]["madeByConformers"]["nConformersGenerated"]
    nConformers = config["chargeFittingInfo"]["nConformers"]
    chargeFittingProtocol = config["chargeFittingInfo"]["chargeFittingProtocol"]

    optMethod  = config["chargeFittingInfo"]["optMethod"]
    optSolvationMethod = config["chargeFittingInfo"]["optSolvationMethod"]
    singlePointMethod = config["chargeFittingInfo"]["singlePointMethod"]
    singlePointSolvationMethod = config["chargeFittingInfo"]["singlePointSolvationMethod"]

    nWaters = None
    if chargeFittingProtocol == "SOLVATOR":
        nWaters = config["runtimeInfo"]["madeByCharges"]["nWaters"]


    if nConformers == -1:
        nConformersUsed = nConformersGenerated
    elif nConformers > nConformersGenerated:
        nConformersUsed = nConformersGenerated
    else:
        nConformersUsed = nConformers

    ## write charge method
    chargesMethod = dedent(f"""
Partial charges for {moleculeName} were generated using the {chargeFittingProtocol} charge-fitting protocol. \
                             """)
    
    if chargeFittingProtocol == "SOLVATOR":
        chargesMethod = chargesMethod + dedent(f"""
ORCA SOLVATOR was used to place {nWaters} explicit water molecules around each conformer, using the stochastic solvation workflow. \
The solvated systems were then evaluated with QM/MM: the solute was treated at the {optMethod} level for geometry optimization and at the {singlePointMethod} level for single-point energies, while solvent waters were modeled as TIP3P at the MM level. \
                                                   """)

    elif chargeFittingProtocol == "RESP":
        chargesMethod = chargesMethod + dedent(f"""
Each conformer was geometry-optimized at the {optMethod} level \
{"with the " + optSolvationMethod + " solvation model" if optSolvationMethod is not None else ""}, \
followed by a single-point calculation at the {singlePointMethod} level \
{"with the " + singlePointSolvationMethod + " solvation model" if singlePointSolvationMethod is not None else ""}. \
                                                   """)

    elif chargeFittingProtocol == "RESP2":
        chargesMethod = chargesMethod + dedent(f"""
Each conformer was geometry-optimized at the {optMethod} level \
{"with the " + optSolvationMethod + " solvation model" if optSolvationMethod is not None else ""}. \
Two {singlePointMethod} evaluations were then run: one with the {singlePointSolvationMethod} solvation model and one in the gas phase. \
                                                   """)

    chargesMethod = chargesMethod + dedent(f"""
This workflow was applied to {nConformersUsed} GOAT-generated conformers of {moleculeName}. \
                                               """)

    chargesMethod = chargesMethod + dedent(f"""
MultiWFN was used for RESP charge fitting, and conformer contributions were Boltzmann-weighted using their relative electronic energies. \
                                                   """)
        
    if chargeFittingProtocol == "RESP2":
        chargesMethod = chargesMethod + dedent(f"""
For RESP2, fitting was performed separately for solvated and gas-phase data, then combined as a 60:40 weighted average (solvated:gas-phase). \
                                               """)

    return chargesMethod


def write_fitting_methods(config: dict) -> str:
    """Write the torsion fitting methods text."""
    ## unpack config
    maxCosineFunctions = config["parameterFittingInfo"]["maxCosineFunctions"]
    shufflesCompleted = config["runtimeInfo"]["madeByStitching"]["shufflesCompleted"]
    fittingMethods = dedent(f"""
MM torsion parameters were optimized iteratively to reproduce the QM torsion profiles from the scanning step. \
For each torsion, the fitted target was defined as:\n\
QM(torsion) = QM(total) - MM(total) + MM(torsion)\n\
where QM(total) comes from the scan data, MM(total) is computed in OpenMM on scan geometries, and MM(torsion) is reconstructed from the current torsion terms. \
Up to {maxCosineFunctions} cosine terms were selected from a Fourier representation of QM(torsion), and the resulting amplitudes, phases, and periodicities were written back to the parameter file. \
This update cycle was repeated for {shufflesCompleted} shuffles of torsion order to reduce order-dependent bias. \
                            """)

    return fittingMethods




def extract_citations(orcaOut: FilePath) -> dict:
    """Extract citation metadata from an ORCA output file."""

    with open(orcaOut, "r") as f:
        fileContent = f.read()
    # Variable names in camelCase
    citationsDict = {}
    linesList = fileContent.splitlines()
    
    inCitationSection = False
    currentCategoryKey = None
    lineIndex = 0

    while lineIndex < len(linesList):
        line = linesList[lineIndex]

        if not inCitationSection:
            if "SUGGESTED CITATIONS FOR THIS RUN" in line:
                inCitationSection = True
            lineIndex += 1
            continue

        # Check for citation category headers
        isNewCategoryFound = False
        if "List of essential papers." in line:
            currentCategoryKey = "essential_papers"
            isNewCategoryFound = True
        elif "List of papers to cite with high priority." in line:
            currentCategoryKey = "high_priority_papers"
            isNewCategoryFound = True
        elif "List of suggested additional citations." in line:
            currentCategoryKey = "suggested_additional_citations"
            isNewCategoryFound = True
        elif "List of optional additional citations" in line: 
            currentCategoryKey = "optional_additional_citations"
            isNewCategoryFound = True
        
        if isNewCategoryFound:
            if currentCategoryKey not in citationsDict: # Initialize list for new category
                citationsDict[currentCategoryKey] = []
            lineIndex += 1
            continue

        if currentCategoryKey: 
            # Attempt to match an author line for a citation
            # Format: "  1. Author(s)"
            authorMatch = re.match(r"^\s*\d+\.\s*(.*)", line)
            if authorMatch:
                authorsText = authorMatch.group(1).strip()
                
                titleText = None
                journalInfoText = None
                doiText = None
                
                # Try to parse Title from the next line
                if lineIndex + 1 < len(linesList):
                    titleLineCandidate = linesList[lineIndex + 1]
                    # Title line must be indented (5 spaces) and non-empty
                    if titleLineCandidate.startswith("     ") and titleLineCandidate.strip():
                        titleText = titleLineCandidate.strip()
                    else:
                        # Title not found or malformed, advance past author line and re-evaluate
                        lineIndex += 1 
                        continue
                else:
                    # End of file after author line
                    lineIndex += 1
                    continue
                    
                # Try to parse Journal Info from the line after title
                if lineIndex + 2 < len(linesList):
                    journalLineCandidate = linesList[lineIndex + 2]
                    if journalLineCandidate.startswith("     ") and journalLineCandidate.strip():
                        journalInfoText = journalLineCandidate.strip()
                    else:
                        # Journal info not found or malformed, advance past title line and re-evaluate
                        lineIndex += 2 
                        continue
                else:
                    # End of file after title line
                    lineIndex += 2
                    continue
                    
                # Try to parse DOI from the line after journal info
                if lineIndex + 3 < len(linesList):
                    doiLineCandidate = linesList[lineIndex + 3]
                    if doiLineCandidate.startswith("     ") and doiLineCandidate.strip():
                        doiText = doiLineCandidate.strip()
                    else:
                        # DOI not found or malformed, advance past journal line and re-evaluate
                        lineIndex += 3 
                        continue
                else:
                    # End of file after journal line
                    lineIndex += 3
                    continue

                # If all parts (title, journal, doi) were successfully parsed for the current author
                citationsDict[currentCategoryKey].append({
                    "authors": authorsText,
                    "title": titleText,
                    "journal_info": journalInfoText,
                    "doi": doiText
                })
                # Advance past the four lines of the citation (author, title, journal, doi)
                lineIndex += 4 
                continue 

        # If the line is not a category header, not a parseable author line, 
        # or if a citation part was malformed causing a 'continue' that advanced lineIndex,
        # this ensures we move to the next line for the next iteration.
        lineIndex += 1
        
    return citationsDict

############################################################################################
def find_orca_output_files(config: dict) -> dict:
    """Locate ORCA output files for all protocol stages."""
    ## find orca output file for capping geom opt
    cappingDir = config["runtimeInfo"]["madeByCapping"]["cappingDir"]
    cappingOrcaOut = p.join(cappingDir, "geometry_optimisation", "orca_opt.out")
    if not p.isfile(cappingOrcaOut):
        cappingOrcaOut = None


    forCapping = [cappingOrcaOut]

    ## find orca output file for GOAT calculations
    conformerDir = config["runtimeInfo"]["madeByConformers"]["conformerDir"]
    conformerOrcaOut = p.join(conformerDir, "GOAT_orca.out")
    if not p.isfile(conformerOrcaOut):
        conformerOrcaOut = None

    forConformers = [conformerOrcaOut]


    ## find orca output files for torsion scans
    torsionDir = config["runtimeInfo"]["madeByTwisting"]["torsionDir"]
    torsionTagDir = [p.join(torsionDir, dirName) for dirName in os.listdir(torsionDir)][0]
    conformerDir = [p.join(torsionTagDir, dirName) for dirName in os.listdir(torsionTagDir) if "conformer" in dirName][0]
    optDir = [p.join(conformerDir, dirName) for dirName in os.listdir(conformerDir) if "opt" in dirName][0]
    torsionScanOut = p.join(optDir, "orca_opt.out")
    if not p.isfile(torsionScanOut):
        torsionScanOut = None

    if not config["torsionScanInfo"]["singlePointMethod"] is None:
        for conformerDir in [p.join(torsionTagDir, dirName) for dirName in os.listdir(torsionTagDir) if "conformer" in dirName]:
            spTopDirs = [p.join(conformerDir, dirName) for dirName in os.listdir(conformerDir) if "SP" in dirName]
            if len(spTopDirs) == 0:
                continue
            spTopDir = spTopDirs[0]
            spDir = [p.join(spTopDir, dirName) for dirName in os.listdir(spTopDir) if "SP" in dirName][0]
            torsionScanSpOut = p.join(spDir, "orca_sp.out")
            break
        if not p.isfile(torsionScanSpOut):
            torsionScanSpOut = None
    else:
        torsionScanSpOut = None

    forTorsions = [torsionScanOut, torsionScanSpOut]

    ## find orca output files for charge calculations
    chargeDir = config["runtimeInfo"]["madeByCharges"]["chargeDir"]
    chargeCalculationProtocol = config["chargeFittingInfo"]["chargeFittingProtocol"]
    #### FOR SOLVATOR PROTOCOL ####
    if chargeCalculationProtocol == "SOLVATOR":
        solvatorDir = p.join(chargeDir, "01_solvator_calculations")
        conformerDir = [p.join(solvatorDir, dirName) for dirName in os.listdir(solvatorDir) if "conformer" in dirName][0]
        solvatorOut = p.join(conformerDir, "SOLVATOR_orca.out")
        if not p.isfile(solvatorOut):
            solvatorOut = None
        ##
        chargeOptDir = p.join(chargeDir, "02_QMMM_optimisations")
        conformerDir = [p.join(chargeOptDir, dirName) for dirName in os.listdir(chargeOptDir) if "conformer" in dirName][0]
        chargeOptOut = p.join(conformerDir, "QMMM_orca_opt.out")
        if not p.isfile(chargeOptOut):
            chargeOptOut = None
        ##
        chargesSpDir = p.join(chargeDir, "03_QMMM_singlepoints")
        conformerDir = [p.join(chargesSpDir, dirName) for dirName in os.listdir(chargesSpDir) if "conformer" in dirName][0]
        chargeSpOut = p.join(conformerDir, "QMMM_orca_sp.out")
        if not p.isfile(chargeSpOut):
            chargeSpOut = None
        
        forCharges = [solvatorOut, chargeOptOut, chargeSpOut]
    #### FOR RESP PROTOCOL ####
    elif chargeCalculationProtocol == "RESP":
        chargeOrcaDir = p.join(chargeDir, "01_orca_calculations")
        conformerDir = [p.join(chargeOrcaDir, dirName) for dirName in os.listdir(chargeOrcaDir) if "conformer" in dirName][0]
        chargeOptOut = p.join(conformerDir, "orca_opt.out")
        if not p.isfile(chargeOptOut):
            chargeOptOut = None
        chargeSpOut = p.join(conformerDir, "orca_sp.out")
        if not p.isfile(chargeSpOut):
            chargeSpOut = None
        forCharges = [chargeOptOut, chargeSpOut]
    #### FOR RESP2 PROTOCOL ####
    elif chargeCalculationProtocol == "RESP2":
        chargeOrcaDir = p.join(chargeDir, "RESP2_gas_phase", "01_orca_calculations")
        conformerDir = [p.join(chargeOrcaDir, dirName) for dirName in os.listdir(chargeOrcaDir) if "conformer" in dirName][0]

        chargeOptOut = p.join(conformerDir, "orca_opt.out")
        if not p.isfile(chargeOptOut):
            chargeOptOut = None
        chargeSpOut = p.join(conformerDir, "orca_sp.out")
        if not p.isfile(chargeSpOut):
            chargeSpOut = None
        forCharges = [chargeOptOut, chargeSpOut]


    orcaOutFiles = {
        "forCapping": forCapping,
        "forConformers": forConformers,
        "forTorsions": forTorsions,
        "forCharges": forCharges
    }

    return orcaOutFiles



def gather_citations(config: dict) -> dict:
    """Collect citation metadata for the final report."""
    ##unpack config
    chargeFittingProtocol = config["chargeFittingInfo"]["chargeFittingProtocol"]
    forcefield = config["parameterFittingInfo"]["forceField"]

    orcaOutFiles = find_orca_output_files(config)

    citations = {}
    for tag in orcaOutFiles:
        for orcaOut in orcaOutFiles[tag]:
            if orcaOut is None:
                continue
            citationsDict = extract_citations(orcaOut)
            citations[tag] = citationsDict
            if tag == "forCharges":
                citations[tag]["essential_papers"].extend(manual_charges_citations(chargeFittingProtocol))
    citations["forStitching"] = {}
    citations["forStitching"] = manual_stitching_citations(forcefield)


    citations = get_last_authors(citations)

    return citations


def get_last_authors(citations: dict) -> dict:
    """Return a reduced citation structure containing only last authors."""
    for tag in citations:
        for level in citations[tag]:
            for paper in citations[tag][level]:
                authors_list = paper["authors"].split(";")
                last_author = authors_list[-1].strip()
                paper["last_author"] = last_author

    return citations



def manual_stitching_citations(forcefield: str) -> list[dict]:
    """Return manual citations for the stitching force field."""
    essential_papers = [{
            "authors": "P. Eastman; R. Galvelis; R. P. Peláez; C. R. A. Abreu; S. E. Farr; E. Gallicchio; A. Gorenko; M. M. Henry; F. Hu; J. Huang; A. Krämer; J. Michel; J. A. Mitchell; V. S. Pande; \
                  J. PGLM Rodrigues; J. Rodriguez-Guerra; A. C. Simmonett; S. Singh; J. Swails; P. Turner; Y. Wang; I. Zhang; J. D. Chodera; G. De Fabritiis; T. E. Markland",
            "title": "OpenMM 8: Molecular Dynamics Simulation with Machine Learning Potentials",
            "journal_info": "J. Phys. Chem. B",
            "doi": "doi.org/10.1021/acs.jpcb.3c06662"      

    }]

    suggested_additional_citations = [
          {         
            "authors": "J A Izaguirre; C R Sweet; V S Pande",
            "title": "Multiscale Dynamics of Macromolecules Using Normal Mode Langevin",
            "journal_info": " Pac Symp Biocomput.",
            "doi": "doi.org/10.1142/9789814295291"     
            },
            {
            "authors": "Charles R. Harris; K. Jarrod Millman; Stéfan J. van der Walt; Ralf Gommers; Pauli Virtanen; David Cournapeau; \
                Eric Wieser; Julian Taylor; Sebastian Berg; Nathaniel J. Smith; Robert Kern; Matti Picus; Stephan Hoyer; Marten H. van Kerkwijk; \
                    Matthew Brett; Allan Haldane; Jaime Fernández del Río; Mark Wiebe; Pearu Peterson; Pierre Gérard-Marchant; \
                        Kevin Sheppard; Tyler Reddy; Warren Weckesser; Hameer Abbasi; Christoph Gohlke; Travis E. Oliphant",
            "title": "Array programming with NumPy",
            "journal_info": "Nature",
            "doi": "doi.org/10.1038/s41586-020-2649-2"     
            }
                ]

    if forcefield == "AMBER":
        essential_papers.append({
            "authors": "David A. Case; Hasan Metin Aktulg; Kellon Belfon; David S. Cerutti; G. Andrés Cisneros; Vinícius Wilian D. Cruzeiro; Negin Forouzesh; Timothy J. Giese; \
                Andreas W. Götz; Holger Gohlke; Saeed Izadi; Koushik Kasavajhala; Mehmet C. Kaymak; Edward King; Tom Kurtzman; Tai-Sung Lee; Pengfei Li; Jian Liu; Tyler Luchko; \
                    Ray Luo; Madushanka Manathunga; Matias R. Machado; Hai Minh Nguyen; Kurt A. O’Hearn; Alexey V. Onufriev; Feng Pan; Sergio Pantano; Ruxi Qi; Ali Rahnamoun; \
                        Ali Risheh; Stephan Schott-Verdugo; Akhil Shajan; Jason Swails; Junmei Wang; Haixin Wei; Xiongwu Wu; Yongxian Wu; Shi Zhang; Shiji Zhao; Qiang Zhu; Thomas E. Cheatham III; \
                            Daniel R. Roe; Adrian Roitberg; Carlos Simmerling; Darrin M. York; Maria C. Nagan; Kenneth M. Merz Jr.",
            "title": "AmberTools",
            "journal_info": "J. Chem. Inf. Model.",
            "doi": "doi.org/10.1021/acs.jcim.3c01153"        
            })
    elif forcefield == "CHARMM":
        essential_papers.append({
            "authors": "K Vanommeslaeghe; E Hatcher; C Acharya; S Kundu; S Zhong; J Shim; E Darian; O Guvench; P Lopes; I Vorobyov; A D MacKerell Jr",
            "title": "CHARMM General Force Field (CGenFF): A force field for drug-like molecules compatible with the CHARMM all-atom additive biological force fields",
            "journal_info": "J. Comput. Chem.",
            "doi": "doi.org/10.1002/jcc.21367"            
        })

    manualStitchingCitations = {
        "essential_papers": essential_papers,
        "suggested_additional_citations": suggested_additional_citations
    }

    return manualStitchingCitations

def manual_charges_citations(chargeFittingProtocol: str) -> list[dict]:
    """Return manual citations for the selected charge fitting protocol."""
    manualChargeCitations = [{
                    "authors": "Tian Lu; Feiwu Chen",
                    "title": "Multiwfn: A Multifunctional Wavefunction Analyzer",
                    "journal_info": "J. Comput. Chem.",
                    "doi": "doi.org/10.1002/jcc.22885"
                            },
                            {
                    "authors": "Jun Zhang; Tian Lu",
                    "title": "Efficient Evaluation of Electrostatic Potential with Computerized Optimized Code",
                    "journal_info": "Phys. Chem. Chem. Phys.",
                    "doi": "doi.org/10.1039/D1CP02805G"
                            },
                            {
                    "authors": "Wendy D. Cornell; Piotr Cieplak; Christopher I. Bayly; Peter A. Kollman",
                    "title": "Application of RESP charges to calculate conformational energies, hydrogen bond energies, and free energies of solvation",
                    "journal_info": "Journal of American Chemical Society",
                    "doi": "doi.org/10.1021/ja00074a030"
                            },
                ]
    
    if chargeFittingProtocol == "RESP2":
        manualChargeCitations.append({
                    "authors": "Michael Schauperl; Paul S. Nerenberg; Hyesu Jang; Lee-Ping Wang; Christopher I. Bayly; David L. Mobley; Michael K. Gilson",
                    "title": "Non-bonded force field model with advanced restrained electrostatic potential charges (RESP2)",
                    "journal_info": "Commun Chem.",
                    "doi": "doi.org/10.1038/s42004-020-0291-4"
        })
    
    return manualChargeCitations




if __name__ == "__main__":
    configFile = "/home/esp/scriptDevelopment/drFrankenstein/04_GLY_CHARMM/drFrankenstein.yaml"

    config = read_input_yaml(configFile)

    methods_writer_protocol(config)

    orcaOutFiles = find_orca_output_files(config)
    gather_citations(config)

    exit()
