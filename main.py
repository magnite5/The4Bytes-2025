import os
from flask import Flask, request, render_template_string
from Bio.PDB import PDBList, PDBParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import seq1
import py3Dmol
import pandas as pd

app = Flask(__name__)

# ========== Helper functions ==========

def download_pdb(pdb_id, out_dir='pdb_files'):
    pdb_id = pdb_id.lower()
    os.makedirs(out_dir, exist_ok=True)
    pdbl = PDBList()
    filepath = pdbl.retrieve_pdb_file(pdb_id, pdir=out_dir, file_format='pdb')
    return filepath

def parse_pdb_sequences(pdb_path):
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
    pa = ProteinAnalysis(seq)
    return {
        'length': len(seq),
        'molecular_weight': round(pa.molecular_weight(), 2),
        'aromaticity': round(pa.aromaticity(), 4),
        'instability_index': round(pa.instability_index(), 2),
        'pI': round(pa.isoelectric_point(), 2)
    }

# ========== Flask routes ==========

@app.route('/')
def home():
    return '''
        <h2>Protein 3D Viewer & Sequence Analyzer</h2>
        <form action="/view">
            <label>Enter PDB ID:</label>
            <input type="text" name="pdb_id" placeholder="e.g. 1A3N" required>
            <button type="submit">Analyze</button>
        </form>
    '''

@app.route('/view')
def view_pdb():
    pdb_id = request.args.get('pdb_id')
    if not pdb_id:
        return "No PDB ID provided."

    try:
        pdb_path = download_pdb(pdb_id)
    except Exception as e:
        return f"Error downloading PDB file: {e}"

    # Parse amino acid sequences
    seqs = parse_pdb_sequences(pdb_path)

    if not seqs:
        return f"No valid sequences found in {pdb_id}."

    # Analyze each chain
    results = []
    for chain_id, seq in seqs.items():
        try:
            analysis = analyze_sequence(seq)
            row = {'Chain': chain_id, **analysis}
            results.append(row)
        except Exception as e:
            print(f"Error analyzing chain {chain_id}: {e}")

    df = pd.DataFrame(results)

    # Generate 3Dmol HTML
    with open(pdb_path, 'r') as fh:
        pdb_text = fh.read()

    view = py3Dmol.view(width=800, height=500)
    view.addModel(pdb_text, 'pdb')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.setBackgroundColor('0xeeeeee')
    view.zoomTo()
    viewer_html = view._make_html()

    # Build HTML page
    table_html = df.to_html(index=False, classes='table table-striped', border=0)

    html = f"""
    <html>
    <head>
        <title>{pdb_id.upper()} Viewer</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                background-color: #fafafa;
                padding: 20px;
            }}
            .table {{
                border-collapse: collapse;
                margin-top: 20px;
            }}
            th, td {{
                border: 1px solid #ccc;
                padding: 8px 12px;
                text-align: center;
            }}
            th {{
                background-color: #f2f2f2;
            }}
        </style>
    </head>
    <body>
        <h1>Protein Structure: {pdb_id.upper()}</h1>
        {viewer_html}
        <h2>Sequence Analysis</h2>
        {table_html}
        <br><a href="/">← Back</a>
    </body>
    </html>
    """

    return render_template_string(html)

# ========== Run Flask app ==========
if __name__ == '__main__':
    app.run(debug=True)
