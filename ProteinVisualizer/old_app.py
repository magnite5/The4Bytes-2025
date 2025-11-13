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

    # Special block for protein analysis
    chain_length = df["length"].sum()
    molecular_weight = df["molecular_weight"].sum()

    count_inst = 0
    for i in range(df.shape[0]):
        current_inst = df.loc[i]["instability_index"]
        current_length = df.loc[i]["length"]
        count_inst += current_inst * current_length
    Weigthed_instability_index = count_inst/chain_length


    count_aromaticity = 0
    for i in range(df.shape[0]):
        current_aromaticity = df.loc[i]["aromaticity"]
        current_length = df.loc[i]["length"]
        count_aromaticity += current_aromaticity * current_length
    Weigthed_aromaticity = count_aromaticity/chain_length

    count_pI = 0
    for i in range(df.shape[0]):
        current_pI = df.loc[i]["pI"]
        current_length = df.loc[i]["length"]
        count_pI += current_pI * current_length
    Weigthed_pI = count_pI/chain_length


    def analyse():
        Strong = 0
        Moderate = 0
        Weak = 0

        value = analyse_chain_length()
        if value == "Strong":
            Strong += 1
        elif value == "Moderate":
            Moderate += 1
        else:
            Weak += 1

        value = analyse_molecular_weight()
        if value == "Strong":
            Strong += 1
        elif value == "Moderate":
            Moderate += 1
        else:
            Weak += 1
        
        value = analyse_aromaticity()
        if value == "Strong":
            Strong += 1
        elif value == "Moderate":
            Moderate += 1
        else:
            Weak += 1
        
        value = analyse_instability_index()
        if value == "Strong":
            Strong += 1
        elif value == "Moderate":
            Moderate += 1
        else:
            Weak += 1
        
        value = analyse_pI()
        if value == "Strong":
            Strong += 1
        elif value == "Moderate":
            Moderate += 1
        else:
            Weak += 1


        if Strong > Moderate and Strong > Weak:
            return "Strong"
        elif Moderate > Strong and Moderate > Weak:
            return "Moderate"
        elif Weak > Strong and Weak > Moderate:
            return "Weak"
        elif Strong == Moderate and Strong > Weak:
            return "Strong"
        elif Strong == Weak and Strong > Moderate:
            return "Strong"
        elif Strong == Weak and Strong == Moderate:
            return "Strong"
        elif Moderate == Weak and Moderate > Strong:
            return "Moderate"
        else:
            return "Strong"


    def analyse_chain_length():
        if chain_length > 400:
            return "Strong"
        elif chain_length >= 150 and chain_length <= 400:
            return "Moderate"
        else:
            return "Weak"

    def analyse_molecular_weight():
        if molecular_weight > 60000:
            return "Strong"
        elif molecular_weight >= 20000 and molecular_weight <= 60000:
            return "Moderate"
        else:
            return "Weak"

    def analyse_aromaticity():
        if Weigthed_aromaticity > 0.10:
            return "Strong"
        elif Weigthed_aromaticity >= 0.07 and Weigthed_aromaticity <= 0.10:
            return "Moderate"
        else:
            return "Weak"

    def analyse_instability_index():
        if Weigthed_instability_index < 30:
            return "Strong"
        elif Weigthed_instability_index >= 30 and Weigthed_instability_index <= 40:
            return "Moderate"
        else:
            return "Weak"

    def analyse_pI():
        if abs((7- Weigthed_pI)) > 2:
            return "Strong"
        elif abs((7-Weigthed_pI)) >= 1 and abs((7-Weigthed_pI)) <= 2:
            return "Moderate"
        else:
            return "Weak"
        
    cumulatif_sentence = ''
    cumulatif_sentence += f'Overall, this protein is {analyse()}. '
    if analyse_chain_length() == 'Weak':
        cumulatif_sentence += f'Because of its small chain length, this protein could be denatured using mild heat or chaotropes. '
    elif analyse_chain_length() == 'Moderate':
        cumulatif_sentence += f'Because of its moderate chain length, this protein could be denatured using moderate to strong denaturants such as 6M urea or sodium dodecyl sulfate combined with heat. '
    else:
        cumulatif_sentence += 'Because of its long chain length, this protein is harder to denature. It is possible to denature it by using a combination of chaotropes, reducing agents and heat. '
    #print(analyse_chain_length())
    if analyse_molecular_weight() == 'Weak':
        cumulatif_sentence += f'Because of its low mollecular weight, this protein could be denatured using mild heat or urea. '
    elif analyse_molecular_weight() == 'Moderate':
        cumulatif_sentence += f'Because of its moderate mollecular weight, this protein could be denatured using moderate to strong denaturants such as 6M urea or sodium dodecyl sulfate combined with heat. '
    else:
        cumulatif_sentence += 'Because of its large mollecular weight, this protein is harder to denature. It is possible to denature it by using guanidium alongside a reducing agent. '
    #print(analyse_molecular_weight())
    if analyse_aromaticity() == 'Weak':
        cumulatif_sentence += f'Because of its low aromaticity, this protein could be denatured using mild heat or a change in pH. '
    elif analyse_aromaticity() == 'Moderate':
        cumulatif_sentence += f'Because of its moderate aromaticity, this protein could be denatured using moderate to strong denaturants such as 6M urea or sodium dodecyl sulfate. '
    else:
        cumulatif_sentence += 'Because of its strong aromaticity, this protein is harder to denature. It is possible to denature it by using guanidium alongside heat. '
    #print(analyse_aromaticity())
    if analyse_instability_index() == 'Weak':
        cumulatif_sentence += f'Because of its high instability, this protein could be denatured using mild heat or a change in pH. '
    elif analyse_instability_index() == 'Moderate':
        cumulatif_sentence += f'Because of its moderate instability, this protein could be denatured using moderate to strong denaturants such as 6M urea or sodium dodecyl sulfate. '
    else:
        cumulatif_sentence += 'Because of its low instability, this protein is harder to denature. It is possible to denature it by using guanidium alongside heat. '
    #print(analyse_instability_index())
    if analyse_pI() == 'Weak':
        cumulatif_sentence += f'Because the isoelectric point is near the pH of 7, this protein could be denatured using mild heat or a change in pH. '
    elif analyse_pI() == 'Moderate':
        cumulatif_sentence += f'Because the isoelectric point is moderately near the pH of 7, this protein could be denatured using moderate to strong denaturants such as 6M urea or sodium dodecyl sulfate. '
    else:
        cumulatif_sentence += 'Because the isoelectric point is far from the pH of 7, this protein is harder to denature. It is possible to denature it by using chaotropes and sodium dodecyl sulfate. '
    #print(analyse_pI())
    if analyse() == 'Weak':
        cumulatif_sentence += f'In summary, this protein is weak to denaturation. It should then be easy to break it using basic techniques such as mild heat, a change in pH or a compound such as urea. '
    elif analyse() == 'Moderate':
        cumulatif_sentence += f'In summary, this protein is moderately resistant to denaturation. It should then be easy to break it using moderate heat combined with a buffer to avoid aggregtion, using chemical denaturants such as sodium dodecyl sulfate or reducing agents if disulfide bonds are present within the protein. '
    else:
        cumulatif_sentence += 'In summary, this protein is resistant to usual denaturation methods. We should then break it using more extreme techniques such as a combination of multiple methods. Exemples include combinations of strong chaotropic agents, detergents such as Sodium dodecyl sulfate, reducing agents, extreme heat or an adjustement to pH. '
    #print(analyse_pI())
    # Generate 3Dmol HTML viewer (built in function from the library)
    with open(pdb_path, 'r') as fh:
        pdb_text = fh.read()

    view = py3Dmol.view(width=1000, height=500)
    view.addModel(pdb_text, 'pdb')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.setBackgroundColor("transparent")

    view.zoomTo()
    viewer_html = view._make_html()
    return render_template(
        'view.html',
        pdb_id=pdb_id,
        viewer_html=viewer_html,
        table_html=table_html,
        description=description,
        cumulatif_sentence=cumulatif_sentence
    )

### Run Flask ###
if __name__ == '__main__':
    app.run(debug=True)
