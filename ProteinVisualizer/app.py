import os
from flask import Flask, request, render_template
from Bio.PDB import PDBList, PDBParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import seq1
import py3Dmol
import pandas as pd
import requests

from utils.analysis_utils import AnalysisUtils

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

@app.route('/test_view')

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

    # Special block for protein analysis
    chain_length = df["length"].sum()
    molecular_weight = df["molecular_weight"].sum()

    count_inst = 0
    for i in range(df.shape[0]):
        current_inst = df.loc[i]["instability_index"]
        current_length = df.loc[i]["length"]
        count_inst += current_inst * current_length
    weighted_instability_index = count_inst/chain_length

    count_aromaticity = 0
    for i in range(df.shape[0]):
        current_aromaticity = df.loc[i]["aromaticity"]
        current_length = df.loc[i]["length"]
        count_aromaticity += current_aromaticity * current_length
    weighted_aromaticity = count_aromaticity/chain_length

    count_pI = 0
    for i in range(df.shape[0]):
        current_pI = df.loc[i]["pI"]
        current_length = df.loc[i]["length"]
        count_pI += current_pI * current_length
    weighted_pi = count_pI/chain_length

    analysis = AnalysisUtils.analyse(
        chain_length, 
        molecular_weight, 
        weighted_aromaticity, 
        weighted_instability_index, 
        weighted_pi
    )
    
    chain_length_analysis = AnalysisUtils.analyse_chain_length(chain_length)
    molecular_weight_analysis = AnalysisUtils.analyse_molecular_weight(molecular_weight)
    aromaticity_analysis = AnalysisUtils.analyse_aromaticity(weighted_aromaticity)
    instability_index_analysis = AnalysisUtils.analyse_instability_index(weighted_instability_index)
    pi_analysis = AnalysisUtils.analyse_pi(weighted_pi)
        
    #cumulative_sentence = ''
    cumulative_sentence = f'Overall, this protein is {analysis}. '
    if chain_length_analysis == 'Weak':
        cumulative_sentence1 = f'Because of its small chain length, this protein could be denatured using mild heat or chaotropes. '
    elif chain_length_analysis == 'Moderate':
        cumulative_sentence1 = f'Because of its moderate chain length, this protein could be denatured using moderate to strong denaturants such as 6M urea, or sodium dodecyl sulfate combined with heat. '
    else:
        cumulative_sentence1 = 'Because of its long chain length, this protein is harder to denature. It is possible to denature it by using a combination of chaotropes, reducing agents and heat. '
    #print(analyse_chain_length())
    if molecular_weight_analysis == 'Weak':
        cumulative_sentence2 = f'Because of its low mollecular weight, this protein could be denatured using mild heat or urea. '
    elif molecular_weight_analysis == 'Moderate':
        cumulative_sentence2 = f'Because of its moderate mollecular weight, this protein could be denatured using moderate to strong denaturants such as 6M urea, or sodium dodecyl sulfate combined with heat. '
    else:
        cumulative_sentence2 = 'Because of its large mollecular weight, this protein is harder to denature. It is possible to denature it by using guanidium alongside a reducing agent. '
    #print(analyse_molecular_weight())
    if aromaticity_analysis == 'Weak':
        cumulative_sentence3 = f'Because of its low aromaticity, this protein could be denatured using mild heat or a change in pH. '
    elif aromaticity_analysis == 'Moderate':
        cumulative_sentence3 = f'Because of its moderate aromaticity, this protein could be denatured using moderate to strong denaturants such as 6M urea, or sodium dodecyl sulfate. '
    else:
        cumulative_sentence3 = 'Because of its strong aromaticity, this protein is harder to denature. It is possible to denature it by using guanidium alongside heat. '
    #print(analyse_aromaticity())
    if instability_index_analysis == 'Weak':
        cumulative_sentence4 = f'Because of its high instability, this protein could be denatured using mild heat or a change in pH. '
    elif instability_index_analysis == 'Moderate':
        cumulative_sentence4 = f'Because of its moderate instability, this protein could be denatured using moderate to strong denaturants such as 6M urea, or sodium dodecyl sulfate. '
    else:
        cumulative_sentence4 = 'Because of its low instability, this protein is harder to denature. It is possible to denature it by using guanidium alongside heat. '
    #print(analyse_instability_index())
    if pi_analysis == 'Weak':
        cumulative_sentence5 = f'Because the isoelectric point is near the pH of 7, this protein could be denatured using mild heat or a change in pH. '
    elif pi_analysis == 'Moderate':
        cumulative_sentence5 = f'Because the isoelectric point is moderately near the pH of 7, this protein could be denatured using moderate to strong denaturants such as 6M urea, or sodium dodecyl sulfate. '
    else:
        cumulative_sentence5 = 'Because the isoelectric point is far from the pH of 7, this protein is harder to denature. It is possible to denature it by using chaotropes and sodium dodecyl sulfate. '
    #print(analyse_pI())
    if analysis == 'Weak':
        cumulative_sentence6 = f'To conclusion, this protein is weak to denaturation. It should then be easy to break it using basic techniques such as mild heat, a change in pH or a compound such as urea. '
    elif analysis == 'Moderate':
        cumulative_sentence6 = f'In conclusion, this protein is moderately resistant to denaturation. It should then be easy to break it using moderate heat combined with a buffer to avoid aggregation, using chemical denaturants such as sodium dodecyl sulfate or reducing agents if disulfide bonds are present within the protein. '
    else:
        cumulative_sentence6 = 'In conclusion, this protein is resistant to usual denaturation methods. We should then break it using more extreme techniques such as a combination of multiple methods. Exemples include combinations of strong chaotropic agents, detergents such as Sodium dodecyl sulfate, reducing agents, extreme heat or an adjustement to pH. '
    #print(analyse_pI())
    # Generate 3Dmol HTML viewer (built in function from the library)
    with open(pdb_path, 'r') as fh:
        pdb_text = fh.read()

    view = py3Dmol.view(width=1000, height=500)
    view.addModel(pdb_text, 'pdb')
    view.setStyle({'cartoon': {'color': '#FFD580'}})  # base color
    view.addStyle({'resi': list(range(1, 50))}, {'cartoon': {'color': "#FF6600"}})
    view.addStyle({'resi': list(range(51, 100))}, {'cartoon': {'color': '#E25822'}})
    view.addStyle({'resi': list(range(101, 150))}, {'cartoon': {'color': "#FF0000"}})

    view.setBackgroundColor("rgba(71, 23, 120, 0)")  # 50% transparent
    view.zoomTo()
    viewer_html = view._make_html()
    
    return render_template(
        'view.html',
        pdb_id=pdb_id,
        viewer_html=viewer_html,
        table_html=table_html,
        description=description,
        cumulative_sentence=cumulative_sentence,
        cumulative_sentence1=cumulative_sentence1,
        cumulative_sentence2=cumulative_sentence2,
        cumulative_sentence3=cumulative_sentence3,
        cumulative_sentence4=cumulative_sentence4,
        cumulative_sentence5=cumulative_sentence5,
        cumulative_sentence6=cumulative_sentence6
    )

### Run Flask ###
if __name__ == '__main__':
    app.run(debug=True)
