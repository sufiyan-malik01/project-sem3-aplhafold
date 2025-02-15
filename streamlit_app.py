# This rebuild of webapp is created by Sufiyan Malik. M.Sc Part-2 .
# The original app is created by Chanin Nantasenamat.(Data Professor) .
# This app is inspired by https://huggingface.co/spaces/osanseviero/esmfold .

import altair as alt
import pandas as pd
import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser,PPBuilder
import io

# st.set_page_config(layout = 'wide')
st.sidebar.title('PROTEINFOLD \nProtein Structure Prediction')
st.sidebar.write('[*ESMFold*](https://esmatlas.com/about) is a protein structure prediction algorithm developed by the META server. It is designed to predict the 3D structure of protein domains using an evolutionary-based approach.\n'
                 '[*ESMFold*](https://esmatlas.com/about) is an end-to-end single sequence protein structure predictor based on the ESM-2 language model. For more information, read the [research article](https://www.biorxiv.org/content/10.1101/2022.07.20.500902v2) and the [news article](https://www.nature.com/articles/d41586-022-03539-1) published in *Nature*.')

# stmol
def render_mol(pdb):
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb,'pdb')
    pdbview.setStyle({'cartoon':{'color':'spectrum'}})
    pdbview.setBackgroundColor('white')#('0xeeeeee')
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    showmol(pdbview, height = 500,width=800)

# Protein sequence input
DEFAULT_SEQ = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
txt = st.sidebar.text_area('Input sequence', DEFAULT_SEQ, height=275)

# ESMfold
def update(sequence=txt):
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
    }
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence, verify=False)
    name = sequence[:3] + sequence[-3:]
    pdb_string = response.content.decode('utf-8')

    with open('predicted.pdb', 'w') as f:
        f.write(pdb_string)

    struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
    b_value = round(struct.b_factor.mean(), 4)

    # Display protein structure
    st.subheader('Visualization of predicted protein structure')
    render_mol(pdb_string)

    # plDDT value is stored in the B-factor field
    st.subheader('plDDT')
    st.write('plDDT is a per-residue estimate of the confidence in prediction on a scale from 0-100.')
    st.info(f'plDDT: {b_value}')

    st.download_button(
        label="Download PDB",
        data=pdb_string,
        file_name='predicted.pdb',
        mime='text/plain',
    )

    ## Protein nucleotide count
    st.header('OUTPUT (Protein Amino Acid Count)')

    ### 1. Print dictionary
    st.subheader('1. Print dictionary')

    def Protein_Aminoacids_count(seq):
        d = dict([
            ('A', seq.count('A')),
            ('R', seq.count('R')),
            ('N', seq.count('N')),
            ('D', seq.count('D')),
            ('C', seq.count('C')),
            ('E', seq.count('E')),
            ('Q', seq.count('Q')),
            ('G', seq.count('G')),
            ('H', seq.count('H')),
            ('I', seq.count('I')),
            ('L', seq.count('L')),
            ('K', seq.count('K')),
            ('M', seq.count('M')),
            ('F', seq.count('F')),
            ('P', seq.count('P')),
            ('S', seq.count('S')),
            ('T', seq.count('T')),
            ('W', seq.count('W')),
            ('Y', seq.count('Y')),
            ('V', seq.count('V')),
        ])
        return d

    X = Protein_Aminoacids_count(sequence)

    # X_label = list(X)
    # X_values = list(X.values())

    X

    ### 2. Print text
    st.subheader('2. Print text')
    st.write('There are ' + str(X['A']) + ' alanine (A)')
    st.write('There are ' + str(X['R']) + ' arginine (R)')
    st.write('There are ' + str(X['N']) + ' asparagaine (N)')
    st.write('There are ' + str(X['D']) + ' cysteine (D)')
    st.write('There are ' + str(X['C']) + ' cysteine (C)')
    st.write('There are ' + str(X['E']) + ' glutamic acid (E)')
    st.write('There are ' + str(X['Q']) + ' glutamine (Q)')
    st.write('There are ' + str(X['G']) + ' glycine (G)')
    st.write('There are ' + str(X['H']) + ' histidine (H)')
    st.write('There are ' + str(X['I']) + ' isoleucine (I)')
    st.write('There are ' + str(X['L']) + ' leucine (L)')
    st.write('There are ' + str(X['K']) + ' lysine (K)')
    st.write('There are ' + str(X['M']) + ' methionine (M)')
    st.write('There are ' + str(X['F']) + ' phenylalanine (F)')
    st.write('There are ' + str(X['P']) + ' proline (P)')
    st.write('There are ' + str(X['S']) + ' serine (S)')
    st.write('There are ' + str(X['T']) + ' threonine (T)')
    st.write('There are ' + str(X['W']) + ' tryptophan (W)')
    st.write('There are ' + str(X['Y']) + ' tyrosine (Y)')
    st.write('There are ' + str(X['V']) + ' valine (V)')

    ### 3. Display DataFrame
    st.subheader('3.Data table of Protein')
    df = pd.DataFrame.from_dict(X, orient='index')
    df = df.rename({0: 'count'}, axis='columns')
    df.reset_index(inplace=True)
    df = df.rename(columns={'index': 'nucleotide'})

    # Add column for amino acid names
    amino_acid_names = {
        'A': 'Alanine (Ala)',
        'R': 'Arginine (Arg)',
        'N': 'Asparagine (Asn)',
        'D': 'Aspartic acid (Asp)',
        'C': 'Cysteine (Cys)',
        'E': 'Glutamic acid (Glu)',
        'Q': 'Glutamine (Gln)',
        'G': 'Glycine (Gly)',
        'H': 'Histidine (His)',
        'I': 'Isoleucine (Ile)',
        'L': 'Leucine (Leu)',
        'K': 'Lysine (Lys)',
        'M': 'Methionine (Met)',
        'F': 'Phenylalanine (Phe)',
        'P': 'Proline (Pro)',
        'S': 'Serine (Ser)',
        'T': 'Threonine (Thr)',
        'W': 'Tryptophan (Trp)',
        'Y': 'Tyrosine (Tyr)',
        'V': 'Valine (Val)'
    }
    df['amino_acid_name'] = df['nucleotide'].map(amino_acid_names)
    st.write(df)

    ### 4. Display Bar Chart using Altair
    st.subheader('4. Bar chart')
    p = alt.Chart(df).mark_bar().encode(
        x='nucleotide',
        y='count'
    )
    p = p.properties(
        width=alt.Step(32)  # controls width of bar.
    )
    st.write(p)

predict = st.sidebar.button('Predict', on_click=update)


if not predict:

    image = Image.open('Screenshot 2023-11-03 213448.png')

    st.image(image, use_column_width=True)

    st.warning('👈 Enter protein sequence data in the query box on left side!')

def calculate_phi_psi_coordinates(pdb_file):
    pdb_text = pdb_file.read().decode('utf-8')
    parser = PDBParser()
    structure = parser.get_structure('protein', io.StringIO(pdb_text))

    phi_angles = []
    psi_angles = []
    coordinates = []
    sequence = ""

    for model in structure:
        for chain in model:
            polypeptides = PPBuilder().build_peptides(chain)
            for poly_index, poly in enumerate(polypeptides):
                phi_psi_angles = poly.get_phi_psi_list()
                aa_sequence = poly.get_sequence()
                sequence += aa_sequence
                for angles in phi_psi_angles:
                    phi, psi = angles
                    if phi is not None and psi is not None:
                        phi_degrees = np.degrees(phi)
                        psi_degrees = np.degrees(psi)
                        phi_angles.append(phi_degrees)
                        psi_angles.append(psi_degrees)
                        coordinates.append((phi_degrees, psi_degrees))

    # Ensure all lists have the same length
    min_length = min(len(phi_angles), len(psi_angles), len(coordinates), len(sequence))
    phi_angles = phi_angles[:min_length]
    psi_angles = psi_angles[:min_length]
    coordinates = coordinates[:min_length]
    sequence = sequence[:min_length]

    return phi_angles, psi_angles, coordinates, sequence

# Streamlit UI
st.title('Ramachandran Plot and Data Table Generator')

# Upload PDB file
pdb_file = st.file_uploader('Upload PDB file', type=['pdb'])

# Button to generate Ramachandran plot and data table
if pdb_file:
    phi_angles, psi_angles, coordinates, sequence = calculate_phi_psi_coordinates(pdb_file)

    # Plot Ramachandran plot
    plt.figure(figsize=(8, 6))
    plt.scatter(phi_angles, psi_angles, s=20, c='b', marker='o', alpha=0.6)
    plt.xlabel('Phi (degrees)')
    plt.ylabel('Psi (degrees)')
    plt.title('Ramachandran Plot')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.grid(True)
    st.set_option('deprecation.showPyplotGlobalUse', False)
    st.pyplot()

    # Create DataFrame
    data = {'Phi (degrees)': phi_angles, 'Psi (degrees)': psi_angles, 'Coordinates': coordinates, 'Sequence': sequence}
    df = pd.DataFrame(data)

    # Display data table
    st.header('Data Table')
    st.table(df)