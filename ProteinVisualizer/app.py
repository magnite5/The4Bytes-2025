import os
from flask import Flask, request, render_template
from Bio.PDB import PDBList, PDBParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import seq1
import py3Dmol
import pandas as pd
import requests

app = Flask(__name__)

### Helper functions ###
def download_pdb(pdb_id, out_dir='output/pdb_files'):
    """Download a PDB file by ID to the given directory."""
    pdb_id = pdb_id.lower()
    os.makedirs(out_dir, exist_ok=True)
    pdbl = PDBList()
    filepath = pdbl.retrieve_pdb_file(pdb_id, pdir=out_dir, file_format='pdb')
    return filepath

def parse_pdb_sequences(pdb_path):
    """Parse amino acid sequences from PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_path)
    sequences = {}

    for model in structure:
        for chain in model:
            seq_chars = []
            for residue in chain:
                # Skip heteroatoms or water
                if residue.id[0] != " ":
                    continue
                try:
                    aa1 = seq1(residue.get_resname())
                    if aa1 != "X":
                        seq_chars.append(aa1)
                except Exception:
                    continue
            if seq_chars:
                sequences[chain.get_id()] = ''.join(seq_chars)
        break  # use first model only
    return sequences

def analyze_sequence(seq):
    """Analyze a protein sequence and return basic statistics."""
    pa = ProteinAnalysis(seq)
    return {
        'length': len(seq),
        'molecular_weight': round(pa.molecular_weight(), 2),
        'aromaticity': round(pa.aromaticity(), 4),
        'instability_index': round(pa.instability_index(), 2),
        'pI': round(pa.isoelectric_point(), 2)
    }

def fetch_pdb_description(pdb_id):
    """Fetch the structure description and PubMed ID from RCSB PDB."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        description = data.get("struct", {}).get("title", "No description available.")

        # Try both citation sources

        return description
    except Exception as e:
        return f"Error fetching description: {e}", None
    
### Flask routes ###

@app.route('/')
def home():
    """Render home page with PDB ID input form."""
    return render_template('home.html')


@app.route('/view')
def view_pdb():
    """Render analysis and 3D structure page for given PDB ID."""
    pdb_id = request.args.get('pdb_id')
    if not pdb_id:
        return "No PDB ID provided."
    try:
        pdb_path = download_pdb(pdb_id)
    except Exception as e:
        return f"Error downloading PDB file: {e}"

    seqs = parse_pdb_sequences(pdb_path)

    if not seqs:
        return f"No valid sequences found in {pdb_id}."

    # Analyze each chain
    results = []
    for chain_id, seq in seqs.items():
        try:
            analysis = analyze_sequence(seq)
            results.append({'Chain': chain_id, **analysis})
        except Exception as e:
            print(f"Error analyzing chain {chain_id}: {e}")

    # Build data table
    df = pd.DataFrame(results)
    table_html = df.to_html(index=False, classes='table table-striped', border=1)
    # Build description
    description = fetch_pdb_description(pdb_id)

    # Generate 3Dmol HTML viewer (built in function from the library)
    with open(pdb_path, 'r') as fh:
        pdb_text = fh.read()

    view = py3Dmol.view(width=1000, height=500)
    view.addModel(pdb_text, 'pdb')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.setBackgroundColor('342c55') 
    view.zoomTo()
    viewer_html = view._make_html()
    return render_template(
        'view.html',
        pdb_id=pdb_id,
        viewer_html=viewer_html,
        table_html=table_html,
        description=description
    )

### Run Flask ###
if __name__ == '__main__':
    app.run(debug=True)
