### The4Bytes-2025
# Protein Visualizer

# Running the "get" code Locally
Requires VSCode and Python 3.x.x

## Virtual Environment
Create a Virtual Environment
```
python3 -m venv venv
```

Set up Virtual Environment (Windows)
```
./venv/Scripts/activate
```


Set up Virtual Environment (Linux / Mac)
```
source venv/bin/activate.(zsh/bash/fish/bat)
```

## Install Dependencies
```
pip install --upgarde pip
pip install flask pandas py3dmodel biopython
```
# Protein Structure Viewer & Analysis
## Overview

This project is a web-based protein structure analysis tool built with Flask and Py3Dmol. It allows users to:

 - Visualize 3D protein structures from the RCSB PDB database.

 - Display protein information including structure description and sequence analysis.

 - Render interactive 3D molecular models with customizable color schemes and transparent or semi-transparent backgrounds.

 - Generate styled tables of sequence or structural data for easy interpretation.

This tool is designed for researchers, educators, and students who want a fast and interactive way to explore protein structures.

## Features

 - Fetch Protein Data

 - Pulls protein structure information from RCSB PDB using the PDB ID.

 - Retrieves weaknesses associated to the protein and means to destroy it

 - 3D Molecular Visualization

 - Uses Py3Dmol to render protein structures in 3D.

 - Supports cartoon-style rendering with custom colors (e.g., orange, spectrum, rainbow).

 - Sequence Analysis

 - Displays sequence or structural data in interactive HTML tables.

 - Tables and 3D viewer are styled via CSS.


