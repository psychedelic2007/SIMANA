import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import numpy as np

def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    tpsa = Descriptors.TPSA(mol)
    wlogp = Descriptors.MolLogP(mol)
    return tpsa, wlogp

def plot_boiled_egg(tpsa_list, wlogp_list):
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Define the regions
    white_region = plt.Circle((2.5, 75), 5, color='yellow', alpha=0.3)
    yolk_region = plt.Circle((2.5, 75), 3, color='white', alpha=0.6)

    ax.add_artist(white_region)
    ax.add_artist(yolk_region)

    # Scatter plot
    ax.scatter(wlogp_list, tpsa_list, c='blue', edgecolor='k', s=100)
    
    # Labeling the molecules
    for i, (wlogp, tpsa) in enumerate(zip(wlogp_list, tpsa_list)):
        ax.annotate(f'Molecule{i+1}', (wlogp, tpsa), fontsize=9)
    
    # Threshold lines
    ax.axhline(y=140, color='r', linestyle='--', label='TPSA = 140')
    ax.axvline(x=5, color='g', linestyle='--', label='WlogP = 5')

    # Labels and title
    ax.set_xlabel('WlogP')
    ax.set_ylabel('TPSA')
    ax.set_title('Boiled Egg Plot')
    ax.legend()

    st.pyplot(fig)


def boiled_egg():
    st.title('Boiled Egg Plot')
    smiles_input = st.text_area('Enter SMILES notation (one per line):')
    if st.button('Plot'):
        smiles_list = smiles_input.split('\n')
        tpsa_list, wlogp_list = [], []
        for smiles in smiles_list:
            tpsa, wlogp = calculate_properties(smiles)
            if tpsa is not None and wlogp is not None:
                tpsa_list.append(tpsa)
                wlogp_list.append(wlogp)
        if tpsa_list and wlogp_list:
            plot_boiled_egg(tpsa_list, wlogp_list)
        else:
            st.error('Invalid SMILES notation entered. Please check the input.')

if __name__ == "__main__":
    boiled_egg()
