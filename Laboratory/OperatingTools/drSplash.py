import os 
from os import path as p
import time
import sys
import re
from tqdm import tqdm

_status_line = None

def set_status_line(status_line: str) -> None:
    """Set the live terminal status line buffer."""
    global _status_line
    _status_line = status_line

def clear_status_line() -> None:
    """Clear the live terminal status line buffer."""
    global _status_line
    _status_line = None

def _truncate_table_tag(tag: str, max_width: int) -> str:
    """Truncate a table label without breaking ANSI-colored output."""

    noAnsiTag = strip_ansi_codes(tag)
    if len(noAnsiTag) <= max_width:
        return tag
    if max_width <= 3:
        return tag[:max_width]
    return f"{tag[:max_width-3]}..."

def _format_torsion_table_cell(tag: str, converged: bool, cell_width: int, score: float | None = None, tol: float = 0.1) -> str:
    """Format one torsion status cell for the progress table."""
    orangeText = "\033[38;5;208m"
    greenText = "\033[32m"
    resetTextColor = "\033[0m"
    scoreText = f"{score:.3f}" if isinstance(score, (int, float)) else "--"
    color = colour_by_score(score, tol=tol) if isinstance(score, (int, float)) else orangeText

    label = f"{tag} [{color}{scoreText}{resetTextColor}]"
    trimmedTag = _truncate_table_tag(label, cell_width)
    padding = " " * (cell_width - len(trimmedTag))
    color = greenText if converged else orangeText
    return f"{trimmedTag}{padding}"

def build_torsion_status_banner() -> str:
    """Build the banner that heads the torsion status display."""
    yellowText = "\033[33m"
    greenText = "\033[32m"
    orangeText = "\033[38;5;208m"
    resetTextColor = "\033[0m"
    return (
        f"{yellowText}{'🗲 '*51}{resetTextColor}\n"
        f"{yellowText}🗲🗲{' '*19}{resetTextColor} "
        f"{resetTextColor}TORSIONS BEING PARAMETERISED: "
        f"[{greenText}converged{yellowText} | {orangeText}not-converged{yellowText}] "
        f"{' '*19}{yellowText}🗲🗲{resetTextColor}"
        )


def colour_by_score(score: float, tol: float) -> str:
    """Choose a color code for a torsion fit score."""
    greenText = "\033[32m"
    orangeText = "\033[38;5;208m"
    redText = "\033[31m"
    if score < tol:
        return greenText
    elif score > tol * 2:
        return redText
    else:
        return orangeText

def build_torsion_status_table_rows(
    torsionTags: list[str],
    convergedTags: set[str],
    scores: dict[str, float] | None = None,
    columns: int = 4,
    ncols: int | None = 102,
    tol: float = 0.1,
) -> list[str]:
    """Build rows for the torsion status table."""
    if not torsionTags:
        return []
    if ncols is None:
        ncols = 102
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"
    separator = f" {yellowText}|{resetTextColor} "
    separatorWidth = 3
    availableWidth = max(ncols - (columns - 1) * separatorWidth, columns)
    cellWidth = max(1, availableWidth // columns)
    rows = []
    for start in range(0, len(torsionTags), columns):
        rowTags = torsionTags[start:start + columns]
        rowCells = [
            _format_torsion_table_cell(tag, tag in convergedTags, cellWidth, None if scores is None else scores.get(tag), tol = tol)
            for tag in rowTags
        ]
        while len(rowCells) < columns:
            rowCells.append(" " * cellWidth)
        rows.append(separator.join(rowCells).rstrip())
    return rows

def init_torsion_status_table(
    torsionTags: list[str],
    convergedTags: set[str] | None = None,
    scores: dict[str, float] | None = None,
    columns: int = 4,
    ncols: int | None = 102,
    position: int = 1,
    tol: float = 0.1,
) -> list[tqdm]:
    """Create the live torsion status table."""
    convergedTags = convergedTags or set()
    rows = build_torsion_status_table_rows(torsionTags, convergedTags, scores, columns, ncols, tol=tol)
    tableRows: list[tqdm] = []
    for rowIndex, rowText in enumerate(rows):
        rowBar = tqdm(
            total=1,
            initial=1,
            desc=rowText,
            bar_format="{desc}",
            position=position + rowIndex,
            leave=True,
            ncols=ncols,
            mininterval=0,
        )
        rowBar.refresh()
        tableRows.append(rowBar)
    return tableRows

def update_torsion_status_table(
    tableRows: list[tqdm],
    torsionTags: list[str],
    convergedTags: set[str],
    scores: dict[str, float] | None = None,
    columns: int = 4,
    ncols: int | None = 102,
    tol: float = 0.1,
) -> None:
    """Refresh the live torsion status table."""
    rows = build_torsion_status_table_rows(torsionTags, convergedTags, scores, columns, ncols, tol=tol)
    for rowBar, rowText in zip(tableRows, rows):
        rowBar.set_description_str(rowText, refresh=True)
        rowBar.refresh()

def close_torsion_status_table(tableRows: list[tqdm]) -> None:
    """Close all torsion status progress bars."""
    for rowBar in tableRows:
        rowBar.close()



def show_mad_man() -> None:
    """Display the main splash banner."""
    # -- ANSI Color Codes --
    redText = "\033[31m"
    y = "\033[33m"
    orangeText = "\033[38;5;208m"
    g = "\033[32m"
    x = "\033[0m"
    p = "\033[35m"

    theMadDocAscii = f"""⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀  {y}  ______           _             ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀  {y}  |  _  \         | |            {g}⠀⠀⠀⠀⠀⠀⠀    ⠲⠿⣶⣀⣦⣠⣤⣠⡀⣀⢀⠀⠀⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀  {y}  | | | |___   ___| |_ ___  _ __ {g}   ⠀⠀⠀⠀⠐⣷⡀⢿⣿⣿⣿⣿⣿⣿⣿⣿⣶⣦⣰⣷⣷⣿⣷⣶⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀{y}⠀⠀| | | / _ \ / __| __/ _ \| '__|{g}⠀⠀⠀    ⢘⡧⠀⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣾⣷⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀  {y}  | |/ / (_) | (__| || (_) | |   {g}         ⡏⠠⢋⣭⣶⣶⣶⣾⣽⣿⠿⣿⣿⣿⣿⣿⡿⠿⠛⠻⠁⠀⠀
⠀⠀⠀⠀{y}⠀ |___/ \___/ \___|\__\___/|_|   {g}⠀⠀⠀⠀⠀⠀⢀⠏⣴⣿⣿⣿⣿⣿⣿⣿⣿⣿⣶⣮⣭⣥⣴⣶⣾⣿⣶⣦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⡀⠀⠀⠀⠀⡌⢸⣿⣿⠿⡛⠻⠛⠻⡿⢿⣿⣿⡞⡟⣽⣿⣿⣿⣿⣿⣿⡇⠀⠀{y}⠀______               _                  _       _       ⠀⠀⠀⠀⠀⠀⠀ {g}⠀⠀⠀⠀⠀⠀
⠀⠀     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⠋⢉⠻⣷⡀⠀⠀⡇⠾⠟⠋⠀⠀⠀{p}⣴⣶⣬⡁{g}⠂⠉⠴⠷⠈⠛⠛{p}⣋⡉{g}⠉⠻⠇{y}⠀⠀⠀|  ___|             | |                | |     (_)  {g}   
⠀⠀     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⡆⡟⠀⠈⠀⠀⠀⢿⣾⣟⣥⡄⠀⠀{p}⠈⠉⠉⠁{g}⢀⣠⣤⣤⠀{p}⠀⠘⠛⠁⠀⠀{y}⠀⠀⠀⠀| |_ _ __ __ _ _ __ | | _____ _ __  ___| |_ ___ _ _ __  ⠀⠀{g}⠀⠀⠀
⠀⠀     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠟⣄⠘⡄⠀⠡⣿⣿⣿⣿⣿⣿⣷⣦⣤⠶⢊⣈⡛⢿⣿⣇⠑⢤⣤⣤⣿⠀⠀⠀⠀{y}⠀|  _| '__/ _` | '_ \| |/ / _ \ '_ \/ __| __/ _ \ | '_ \ ⠀⠀⠀⠀⠀⠀⠀⠀{g}⠀⠀⠀
⠀⠀     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢈⣣⡄⠀⠀⣼⣿⣿⠿⠟⠛⢋⣩⣴⣦⡘⠛⠛⠿⣿⣿⠈⠀⢻⣿⡟⠀⠀⠀⠀{y}⠀| | | | | (_| | | | |   <  __/ | | \__ \ ||  __/ | | | |⠀⠀⠀⠀⠀⠀⠀⠀{g}⠀⠀⠀
⠀⠀     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠀⠀⣿⡟⣡⣄⠈⠙⠛⠿⢿⣿⣿⣿⣦⣀⣉⣠⣤⠾⠈⠋⠀⠀    {y}\_| |_|  \__,_|_| |_|_|\_\___|_| |_|___/\__\___|_|_| |_|⠀⠀⠀⠀⠀⠀⠀⠀{g}⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣰⣾⡟⠁⠘⡄⠀⢳⡄{g}⠀⢽⣿⣿⣮⡻{x}⡢⡙⠿⢶⣮⣼⣭⡏⡿{g}⠀⣼⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⣀⣠⣤⣴⡆⣸⣿⣿⠀⠀⠀⠈⠀⠈⢿⣦⡀{g}⠙⠛⢿⣿⣮⣓⠀⠀⠐⠒⠒⠐⠀⣼⠏{x}⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣤⣤⣶⣾⣿⣿⣿⣿⣿⣿⣿⢡⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀⠙⣿⣦⡀{g}⠀⠈⠛⠿⣿⣶⣶⣶⣶⣶⣿⠏{x}⣼⣌⣄⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀     ⠀⠀⠀⠀⠀⠀⠀⠀⢀⡀⢀⢲⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡏⣾⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠙⠓⠤⠀{g}⠀⠀⠈⠉⠉⠉⠉⠁{x}⠈⢿⣿⣿⣾⣿⣶⣤⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀     ⠀⠀⠀⠀⠀⢀⣴⣷⣦⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⢣⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⣿⣿⣿⢽⣿⣿⣿⣿⣷⢆⡤⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀     ⠀⠀⠀⠀⢐⣡⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⣸⣿⣿⣿⣿⡀⠀⠀⠀⠀⠀⠀⠀{y}LIGHTNING-FAST{x}⠀⠀⢹⣿⣿⣯⢻⣿⣿⣟⣵⣿⣾⣿⣮⣶⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀     ⠀⠀⠀⢀⣿⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇⣿⣿⣿⣿⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀{y}PARAMETERS{x}⠀⠀⠀⠀⠀⠀⣿⣿⣿⣷⣻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀     ⠀⠀⢠⢏⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣇⠿⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣷⣻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀     ⢀⢀⡟⣼⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡎⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣦⣝⠻⣿⣿⡆⠀⠀⠀⠀   ⠀⠀  FOR⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⣿⣿⣿⣷⣻⣿⣿⣿⣿⣿⣿⣿⢸⣿⣷⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀     ⠀⣾⢱⣿⣿⣿⣿⣿⣿⣿⣾⢿⣿⣿⣿⡘⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣌⠻⣷⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⢟⣽⣿⣿⣿⣿⣿⣿⣿⣸⣿⣿⣿⡀⠀⠀⠀⠀⠀⠀    ⠀⠀⠀⠀⠀⠀⠀⠀     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀     ⣾⢹⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣻⣿⣿⠘⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣟⣵⣿⣿⣿⣧⠀⠀⠀⠀⠀⠀ {p} MONSTROUS{x}⠀⠀⠀⠀⠀⠀⠀⣿⣿⣿⣾⣿⣿⣿⣿⣿⣿⣿⡇⣿⣿⣷⣶⣿⠀⠀⠀⠀⠀          ⠀⠀      ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀     ⣿⡉⣹⣿⣿⣿⣿⣿⣿⣿⡻⣿⣿⣧⢻⣿⠀⣿⣿⣿⣿⣿⣿⣿⣿⣿⣏⢻⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀ {p} MOLECULES{x}⠀⠀⠀⠀⠀⠀⠀⢻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⢻⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀       ⠀⠀⠀⠀      ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀     ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣎⢿⣿⡎⡿⢠⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣆⢻⣿⣿⣿⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⣼⠿⣛⣿⣿⣿⠀⠀⠀⠀⠀         ⠀⠀⠀      ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀{x}⠀    ⠀⣼⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣯⢿⣇⠁⣼⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣆⢻⣿⣿⣿⣧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⠟⣯⣱⣶⣿⣷⣯⡺⢟⠃             ⠀⠀⠀     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    """
    lightningBar = "🗲 "*int((140//2))
    print(y+lightningBar+x)
    print(theMadDocAscii)
    print(y+lightningBar+x)

    time.sleep(1)

def print_config_error(configErrors: dict) -> None:
    """Print a detailed, color-coded error report for the config file and exit."""
    # -- ANSI Color Codes --
    redText = "\033[31m"
    yellowText = "\033[33m"
    orangeText = "\033[38;5;208m"
    greenText = "\033[32m"
    resetTextColor = "\033[0m"

    configErrorAscii = """
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
    lightningBar = "🗲 "*int((80//2))


    print(f"{resetTextColor}Fatal errors were found in your config file. Please fix them to proceed.")
    print(f"{resetTextColor}Colour Key: | {greenText}Input Correct{resetTextColor} | {orangeText}Non-Fatal Issue (Default Used){resetTextColor} | {redText}Fatal Issue{resetTextColor} |")

    # -- Helper function to print a single line --
    def print_config_text(argName: str, argDisorder: str, textColor: str, indentationLevel: int = 0) -> None:
        print(f"{' ' * (indentationLevel * 3 + 2)}{yellowText}{argName}: {textColor}{argDisorder}{resetTextColor}")

    # -- Helper function to recursively print nested dictionaries --
    def loop_error_dict(argName: str, disorderDict: dict, indentationLevel: int = 0) -> None:
        print(f"{'--' * indentationLevel}--> In sub-entry {yellowText}{argName}{resetTextColor}:")
        for sub_argName, sub_argDisorder in disorderDict.items():
            if sub_argDisorder is None:
                print_config_text(sub_argName, "OK", greenText, indentationLevel)
            elif isinstance(sub_argDisorder, str):
                textColor = orangeText if "default" in sub_argDisorder.lower() else redText
                print_config_text(sub_argName, sub_argDisorder, textColor, indentationLevel)
            elif isinstance(sub_argDisorder, dict):
                loop_error_dict(sub_argName, sub_argDisorder, indentationLevel + 1)

    # -- Main loop to process the error dictionary --
    for infoName, infoErrors in configErrors.items():
        print(f"\n> For the config entry {yellowText}{infoName}{resetTextColor}:")
        for argName, argDisorder in infoErrors.items():
            if argDisorder is None:
                print_config_text(argName, "OK", greenText, 0)
            elif isinstance(argDisorder, str):
                textColor = orangeText if "default" in argDisorder.lower() else redText
                print_config_text(argName, argDisorder, textColor, 0)
            elif isinstance(argDisorder, dict):
                loop_error_dict(argName, argDisorder)

    print(f"\n{redText}⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕")
    print(resetTextColor)
    exit(1)




def strip_ansi_codes(text: str) -> str:
    """Remove ANSI color codes from a string."""
    ansiPattern = re.compile(r'\033\[[0-9;]*m')
    return ansiPattern.sub('', text)

def show_getting_mm_total(torsionTag: str, inline: bool = True) -> None:
    """Show the status message for MM torsion energy calculation."""
    greenText = "\033[32m"
    redText = "\033[31m"
    orangeText = "\033[38;5;172m"
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"
    tealText = "\033[38;5;37m" 
    brightGreenText = "\033[92m"

    text = f"{yellowText}🗲🗲{resetTextColor}{' '*8}Fitting torsion parameters for {greenText}{torsionTag}"
  # Calculate visible length (excluding ANSI codes)
    visibleText = strip_ansi_codes(text)
    nSpaces = 102 - len(visibleText) - 3  # Account for trailing spaces and 🗲🗲
    text += f"{' '*nSpaces}{yellowText}🗲🗲{resetTextColor}"

    if _status_line is not None:
        _status_line.set_description_str(text, refresh=True)
        _status_line.refresh()
        return
    if inline:
        sys.stdout.write("\r\033[K" + text)  # \r: move to start of line, \033[K: clear line
        sys.stdout.flush()
    else:
        tqdm.write(text)



def show_need_cgenff_str(cappedMol2: str) -> None:
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


def show_capping_splash() -> None:
    ## using font Delta Corps Priest 1
    greenText = "\033[32m"
    redText = "\033[31m"
    orangeText = "\033[38;5;172m"
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"

    lightningBar = "🗲 "*int((84//2))

    cappingSplash = """
 ▄████████    ▄████████    ▄███████▄    ▄███████▄  ▄█  ███▄▄▄▄      ▄██████▄    ▄█ 
███    ███   ███    ███   ███    ███   ███    ███ ███  ███▀▀▀██▄   ███    ███  ███ 
███    █▀    ███    ███   ███    ███   ███    ███ ███▌ ███   ███   ███    █▀   ███▌
███          ███    ███   ███    ███   ███    ███ ███▌ ███   ███  ▄███         ███▌
███        ▀███████████ ▀█████████▀  ▀█████████▀  ███▌ ███   ███ ▀▀███ ████▄   ███
███    █▄    ███    ███   ███          ███        ███  ███   ███   ███    ███  █▀ 
███    ███   ███    ███   ███          ███        ███  ███   ███   ███    ███   
████████▀    ███    █▀   ▄████▀       ▄████▀      █▀    ▀█   █▀    ████████▀   ██
                                                                              
"""
    print(yellowText + lightningBar +
           resetTextColor + cappingSplash  +
             yellowText + " "*24 + "ADDING CAPPING GROUPS...\n" + 
             lightningBar + resetTextColor)
    
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
def show_torsion_being_scanned(torsionTag: str, torsionIndex: int, nTorsions: int) -> None:
    greenText = "\033[32m"
    redText = "\033[31m"
    orangeText = "\033[38;5;172m"
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"
    tealText = "\033[38;5;37m"    

    textToPrint = f"{yellowText}🗲 🗲{' '*8}SCANNING TORSION:{' '*8}\
[{' '*1}{greenText}{torsionTag}{yellowText}{' '*1}]{' '*8}\
{resetTextColor}({torsionIndex+1}/{nTorsions}){yellowText}"
    
    visibleText = strip_ansi_codes(textToPrint)

    
    nSpaces = 124 - len(visibleText) - 3  # Account for trailing spaces and 🗲🗲
    textToPrint += f"{' '*nSpaces}{yellowText}🗲🗲{resetTextColor}"
    print(textToPrint)





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

def show_what_have_we_created(config: dict) -> None:
    """Display a summary of the generated outputs."""
    ## unpack config ##
    moleculeName = config["moleculeInfo"]["moleculeName"]
    reportHtml = config["runtimeInfo"]["madeByReporting"]["reportHtml"]
    forcefield = config["parameterFittingInfo"]["forceField"]


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
    def clear_terminal() -> None:
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
    print(" "*16+ yellowText + f"{forcefield} PARAMETERS FOR YOUR MOLECULE:" + " "*4 + green_text + moleculeName + resetTextColor)
    print()
    print(" "*16+ yellowText + f"VIEW REPORT AT:\n" + " "*4 + f"{green_text}{reportHtml}{resetTextColor}")




#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def print_botched(errorReport: dict) -> None:
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
    show_mad_man()
    show_capping_splash()
    show_wriggle_splash()
    show_twist_splash()
    show_charge_splash()
    show_torsion_being_scanned("N-C-CA-CB", 1, 10)
    show_stitch_splash()
    show_creation_splash()

    config = {
        "moleculeInfo": {
            "moleculeName": "ALA"
        },
        "runtimeInfo": {
            "madeByReporting": {
                "reportHtml": "/path/to/report/html"
            },
      
        },
        "parameterFittingInfo": {
            "forceField": "AMBER"
            }
    }


    show_what_have_we_created(config)
