import os
import requests
from Bio.PDB import PDBList, PDBParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import py3Dmol
import pandas as pd

### FOR LOCAL TESTING, USE main.ipynb FILE, SEE README ###

# Helper function to download PDB file
def download_pdb(pdb_id, out_dir='pdb_files'):
    pdb_id = pdb_id.lower()  # convert to lowercase as PDB IDs are case insensitive
    os.makedirs(out_dir, exist_ok=True)
    pdbl = PDBList()
    # Retrieve and download PDB file
    filepath = pdbl.retrieve_pdb_file(pdb_id, pdir=out_dir, file_format='pdb')
    return filepath

# Helper function to parse the PDB file and extract amino acid sequences
def parse_pdb_sequences(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_path)
    print(structure)
    sequences = {}

    for model in structure:
        for chain in model:
            seq_chars = []
            for residue in chain:
                # Only process standard amino acids (skip 'X' or unknown residues)
                #if residue.get_id()[0] == ' ':  # ' ' indicates standard residue
                    #try:
                aa = residue.get_resname()
                        # Convert three-letter code to one-letter code using Biopython's helper function
                from Bio.SeqUtils import seq1
                aa1 = seq1(aa)
                if aa1 == 'X':
                    continue
                   # except Exception:
                        #continue  # Skip unknown residues
                #print(aa)
                seq_chars.append(aa1)
            if seq_chars:  # Only add non-empty sequences
                sequences[chain.get_id()] = ''.join(seq_chars)
        break  # Usually only need the first model
    return sequences

# Helper function to compute basic sequence properties
def analyze_sequence(seq):
    #print(seq)
    pa = ProteinAnalysis(seq)
    return {
        'length': len(seq),
        'molecular_weight': pa.molecular_weight(),
        'aromaticity': pa.aromaticity(),
        'instability_index': pa.instability_index(),
        'isoelectric_point': pa.isoelectric_point(),
        'aa_composition': pa.get_amino_acids_percent()
    }

# Download PDB file for Hemoglobin (ID: 1A3N)
PDB_ID = "2ZFO"  # Human Hemoglobin Alpha Chain
pdb_path = download_pdb(PDB_ID, out_dir='pdb_files')
print("Downloaded PDB file to:", pdb_path)

# Parse the sequences from the downloaded PDB file
seqs = parse_pdb_sequences(pdb_path)
print(seqs)
print("Chains and sequences found:", list(seqs.keys()))

# Analyze each sequence (for all chains found)
results = []
for chain_id, seq in seqs.items():
    if len(seq) == 0:
        print(f"Chain {chain_id} is empty, skipping...")
        continue  # Skip empty chains

    # Perform sequence analysis only if the sequence is non-empty
    try:
        analysis = analyze_sequence(seq)
        row = {
            'chain': chain_id,
            'length': analysis['length'],
            'molecular_weight': round(analysis['molecular_weight'], 2),
            'aromaticity': round(analysis['aromaticity'], 4),
            'instability_index': round(analysis['instability_index'], 2),
            'pI': round(analysis['isoelectric_point'], 2)
        }
        results.append(row)
        print(f"Chain {chain_id}: length {row['length']}, MW {row['molecular_weight']}, pI {row['pI']}")
    except Exception as e:
        print(f"Error analyzing chain {chain_id}: {e}")

# Create a DataFrame to summarize the results
df = pd.DataFrame(results)
df

# 3D visualization of the Hemoglobin protein (you can zoom in/out interactively)
with open(pdb_path, 'r') as fh:
    pdb_text = fh.read()

view = py3Dmol.view(width=800, height=500)
view.addModel(pdb_text, 'pdb')
view.setStyle({'cartoon': {'color':'spectrum'}})   # Display in cartoon style (color-coded by structure)
view.setBackgroundColor('0xeeeeee')  # Light gray background
view.zoomTo()  # Auto zoom to fit the structure
view.show()