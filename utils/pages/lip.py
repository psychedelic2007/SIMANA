import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from io import StringIO, BytesIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.preprocessing import MinMaxScaler

# Function to calculate Lipinski and other properties
def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # Lipinski parameters
    mol_weight = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    num_h_donors = Descriptors.NumHDonors(mol)
    num_h_acceptors = Descriptors.NumHAcceptors(mol)
    
    # Additional properties
    ring_count = Descriptors.RingCount(mol)
    psa = Descriptors.TPSA(mol)
    num_aromatic_rings = len([ring for ring in mol.GetAromaticAtoms() if ring.IsInRing()])
    num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    
    # Atom distribution
    atom_dist = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_dist[symbol] = atom_dist.get(symbol, 0) + 1
    
    # Violations
    violations = []
    if mol_weight > 500: violations.append("MolWt > 500")
    if num_h_donors > 5: violations.append("HDonors > 5")
    if num_h_acceptors > 10: violations.append("HAcceptors > 10")
    if logp > 5: violations.append("LogP > 5")
    follows_rule = "No" if violations else "Yes"
    
    return {
        "SMILES": smiles,
        "MolWt": mol_weight,
        "LogP": logp,
        "PSA": psa,
        "RingCount": ring_count,
        "HDonors": num_h_donors,
        "HAcceptors": num_h_acceptors,
        "AromaticRings": num_aromatic_rings,
        "RotatableBonds": num_rotatable_bonds,
        "FollowsLipinski": follows_rule,
        "Violations": ", ".join(violations) if violations else "--",
        "AtomDistribution": atom_dist,
        "Mol": mol  # Save RDKit molecule object
    }

# Plot distributions for properties
def plot_distributions(data):
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    sns.histplot(data['MolWt'], ax=axs[0, 0], kde=True).set(title='Molecular Weight')
    sns.histplot(data['HDonors'], ax=axs[0, 1], kde=True).set(title='H Donors')
    sns.histplot(data['HAcceptors'], ax=axs[0, 2], kde=True).set(title='H Acceptors')
    sns.histplot(data['LogP'], ax=axs[1, 0], kde=True).set(title='LogP')
    sns.histplot(data['RingCount'], ax=axs[1, 1], kde=True).set(title='Ring Count')
    sns.histplot(data['PSA'], ax=axs[1, 2], kde=True).set(title='Polar Surface Area')
    plt.tight_layout()
    return fig

# Radar plot function
def plot_radar_normalized(selected_compound):
    properties = {
        'RingCount': selected_compound['RingCount'],
        'HDonors': selected_compound['HDonors'],
        'HAcceptors': selected_compound['HAcceptors']
    }

    # Normalize properties to a 0-1 scale
    scaler = MinMaxScaler()
    normalized_values = scaler.fit_transform(np.array(list(properties.values())).reshape(-1, 1)).flatten()

    # Radar plot setup
    labels = list(properties.keys())
    angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
    normalized_values = np.concatenate((normalized_values, [normalized_values[0]]))
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))
    ax.fill(angles, normalized_values, color='blue', alpha=0.25)
    ax.plot(angles, normalized_values, color='blue', linewidth=2)
    ax.set_yticklabels([])
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels)

    return fig

def lip():
    st.title("Lipinski's Rule of Five Calculator")
    uploaded_file = st.file_uploader("Upload a text file containing SMILES notations", type=["txt"])

    if uploaded_file:
        smiles_list = uploaded_file.read().decode("utf-8").splitlines()
        compounds = [calculate_properties(smiles.strip()) for smiles in smiles_list if smiles.strip()]
        compound_names = [f"Compound {i+1}" for i in range(len(compounds))]
        
        if st.button("Analyze Compounds"):
            data = pd.DataFrame([{
                "SMILES": c["SMILES"],
                "MolWt": c["MolWt"],
                "LogP": c["LogP"],
                "PSA": c["PSA"],
                "RingCount": c["RingCount"],
                "HDonors": c["HDonors"],
                "HAcceptors": c["HAcceptors"]
            } for c in compounds])
            
            st.session_state['processed_data'] = data
            st.session_state['plot_figure'] = plot_distributions(data)
            st.pyplot(st.session_state['plot_figure'])

        st.download_button("Download Lipinski Data CSV",
                           data=st.session_state['processed_data'].to_csv(index=False),
                           file_name="lipinski_data.csv",
                           mime="text/csv")

        st.subheader("Select a compound for detailed analysis")
        selected_compound_name = st.selectbox("Choose a compound", compound_names)
        selected_compound = compounds[compound_names.index(selected_compound_name)]

        # Display compound details using st.columns without 'with' statement
        col1, col2, col3 = st.columns([1, 1, 1])

        # Column 1: Compound details
        with col1:
            st.write(f"**Molecular Weight:** {selected_compound['MolWt']:.2f}")
            st.write(f"**LogP:** {selected_compound['LogP']:.2f}")
            st.write(f"**Polar Surface Area:** {selected_compound['PSA']:.2f}")
            st.write(f"**Ring Count:** {selected_compound['RingCount']}")
            st.write(f"**Follows Lipinski's Rule:** {selected_compound['FollowsLipinski']}")
            st.write(f"**Violations:** {selected_compound['Violations']}")
            st.write(f"**Atom Distribution:** {selected_compound['AtomDistribution']}")
            
        # Column 2: Display molecule's 2D structure
        with col2:
            img = Draw.MolToImage(selected_compound['Mol'], size=(500, 500), dpi=1200)
            st.image(img, caption="2D Structure", use_column_width=True)
            
        # Column 3: Radar plot
        with col3:
            radar_fig = plot_radar_normalized(selected_compound)
            st.pyplot(radar_fig)

if __name__ == "__main__":
    lip()
