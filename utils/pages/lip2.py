import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import IPythonConsole
from io import StringIO, BytesIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.preprocessing import MinMaxScaler

import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import IPythonConsole
from io import StringIO, BytesIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.preprocessing import MinMaxScaler

import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import IPythonConsole
from io import StringIO, BytesIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.preprocessing import MinMaxScaler

def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # Calculate all properties with consistent naming
    properties = {
        "MW": Descriptors.ExactMolWt(mol),
        "nBonds": mol.GetNumBonds(),
        "fChar": Chem.GetFormalCharge(mol),
        "nHet": rdMolDescriptors.CalcNumHeteroatoms(mol),
        "MaxRing": max([len(ring) for ring in mol.GetRingInfo().AtomRings()]) if mol.GetRingInfo().AtomRings() else 0,
        "nRing": rdMolDescriptors.CalcNumRings(mol),
        "nRot": Descriptors.NumRotatableBonds(mol),
        "TPSA": Descriptors.TPSA(mol),
        "nHD": Descriptors.NumHDonors(mol),
        "nHA": Descriptors.NumHAcceptors(mol),
        "LogP": Descriptors.MolLogP(mol),  # Changed from logP to LogP
        "LogD": Descriptors.MolLogP(mol),  # Changed from logD to LogD
        "LogS": Descriptors.MolLogP(mol) - 0.89,  # Changed from logS to LogS
        "SC": len(Chem.FindMolChiralCenters(mol))
    }
    
    # Add violations check with consistent naming
    violations = []
    if properties["MW"] > 500: violations.append("MolWt > 500")
    if properties["nHD"] > 5: violations.append("HDonors > 5")
    if properties["nHA"] > 10: violations.append("HAcceptors > 10")
    if properties["LogP"] > 5: violations.append("LogP > 5")  # Updated to match property name
    
    properties["FollowsLipinski"] = "No" if violations else "Yes"
    properties["Violations"] = ", ".join(violations) if violations else "--"
    
    # Add atom distribution
    atom_dist = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_dist[symbol] = atom_dist.get(symbol, 0) + 1
    properties["AtomDistribution"] = atom_dist
    
    properties["Mol"] = mol  # Store molecule object
    return properties

def plot_distributions(data):
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    sns.histplot(data['MW'], ax=axs[0, 0], kde=True).set(title='Molecular Weight')
    sns.histplot(data['nHD'], ax=axs[0, 1], kde=True).set(title='H Donors')
    sns.histplot(data['nHA'], ax=axs[0, 2], kde=True).set(title='H Acceptors')
    sns.histplot(data['LogP'], ax=axs[1, 0], kde=True).set(title='LogP')  # Updated from logP to LogP
    sns.histplot(data['nRing'], ax=axs[1, 1], kde=True).set(title='Ring Count')
    sns.histplot(data['TPSA'], ax=axs[1, 2], kde=True).set(title='Polar Surface Area')
    plt.tight_layout()
    return fig

def plot_radar_normalized(selected_compound):
    # Define property ranges with consistent naming
    ranges = {
        "MW": (160, 500),
        "nBonds": (0, 50),
        "fChar": (-2, 2),
        "nHet": (0, 10),
        "MaxRing": (0, 7),
        "nRing": (0, 4),
        "nRot": (0, 10),
        "TPSA": (0, 140),
        "nHD": (0, 5),
        "nHA": (0, 10),
        "LogD": (-3, 5),  # Updated from logD to LogD
        "LogS": (-6, 1),  # Updated from logS to LogS
        "LogP": (-3, 5)   # Updated from logP to LogP
    }
    
    # Normalize values between 0 and 1
    normalized_values = {}
    for prop, (min_val, max_val) in ranges.items():
        if prop in selected_compound:
            val = selected_compound[prop]
            normalized_values[prop] = max(0, min(1, (val - min_val) / (max_val - min_val)))
    
    # Setup the radar plot
    labels = list(ranges.keys())
    num_vars = len(labels)
    
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    values = [normalized_values[key] for key in labels]
    
    # Close the plot by appending the first value
    values += values[:1]
    angles += angles[:1]
    labels += labels[:1]
    
    # Create figure and polar axis
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
    
    # Plot data
    ax.plot(angles, values, 'o-', linewidth=2, label='Compound Properties', color='orange')
    ax.fill(angles, values, alpha=0.25, color='orange')
    
    # Add lower and upper limits
    lower_limits = [0] * len(labels)
    upper_limits = [1] * len(labels)
    ax.plot(angles, lower_limits, 'o-', linewidth=1, label='Lower Limit', color='green', alpha=0.5)
    ax.fill(angles, lower_limits, alpha=0.1, color='green')
    ax.plot(angles, upper_limits, 'o-', linewidth=1, label='Upper Limit', color='blue', alpha=0.5)
    ax.fill(angles, upper_limits, alpha=0.1, color='blue')
    
    # Set chart properties
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels[:-1])
    ax.set_ylim(0, 1)
    
    # Add legend
    ax.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
    
    plt.title("Molecular Properties Radar Plot", pad=20)
    return fig

def lip2():
    st.title("Lipinski's Rule of Five Calculator")
    
    # Initialize session state variables
    if 'compounds' not in st.session_state:
        st.session_state.compounds = None
    if 'compound_names' not in st.session_state:
        st.session_state.compound_names = None
    if 'processed_data' not in st.session_state:
        st.session_state.processed_data = None
    if 'analyzed' not in st.session_state:
        st.session_state.analyzed = False

    # Input method selection
    input_method = st.radio("Choose input method:", ["Upload File", "Manual Input"])
    
    smiles_list = []
    
    if input_method == "Upload File":
        uploaded_file = st.file_uploader("Upload a text file containing SMILES notations", type=["txt"])
        if uploaded_file is not None:
            smiles_list = uploaded_file.read().decode("utf-8").splitlines()
            smiles_list = [s.strip() for s in smiles_list if s.strip()]
    else:
        smiles_input = st.text_area("Enter SMILES notations (one per line)")
        if smiles_input:
            smiles_list = [s.strip() for s in smiles_input.splitlines() if s.strip()]
    
    # Analyze button
    if st.button("Analyze Compounds"):
        if smiles_list:
            # Process SMILES and calculate properties
            compounds = [calculate_properties(smiles) for smiles in smiles_list]
            compounds = [c for c in compounds if c is not None]
            
            if not compounds:
                st.error("No valid SMILES notations found. Please check your input.")
                return
            
            compound_names = [f"Compound {i+1}" for i in range(len(compounds))]
            
            # Store in session state
            st.session_state.compounds = compounds
            st.session_state.compound_names = compound_names
            st.session_state.processed_data = pd.DataFrame([{
                "SMILES": smiles_list[i],
                "MW": c["MW"],
                "LogP": c["LogP"],
                "TPSA": c["TPSA"],
                "nRing": c["nRing"],
                "nHD": c["nHD"],
                "nHA": c["nHA"]
            } for i, c in enumerate(compounds)])
            st.session_state.analyzed = True

    # Show results if analysis has been performed
    if st.session_state.analyzed and st.session_state.compounds:
        # Only show distribution plots if there are multiple compounds
        if len(st.session_state.compounds) > 1:
            st.subheader("Property Distributions")
            fig = plot_distributions(st.session_state.processed_data)
            st.pyplot(fig)
            
            # Download button (only for multiple compounds)
            st.download_button(
                "Download Lipinski Data CSV",
                data=st.session_state.processed_data.to_csv(index=False),
                file_name="lipinski_data.csv",
                mime="text/csv"
            )

        # For multiple compounds, show compound selection
        if len(st.session_state.compounds) > 1:
            st.subheader("Select a compound for detailed analysis")
            selected_compound_name = st.selectbox(
                "Choose a compound",
                st.session_state.compound_names,
                key='compound_selector'
            )
            selected_idx = st.session_state.compound_names.index(selected_compound_name)
        else:
            # For single compound, just use the first one
            selected_idx = 0

        selected_compound = st.session_state.compounds[selected_idx]

        # Display compound details in columns
        col1, col2, col3 = st.columns([1, 1, 1])
        
        with col1:
            st.write("**Basic Properties:**")
            st.write(f"**Molecular Weight:** {selected_compound['MW']:.2f}")
            st.write(f"**LogP:** {selected_compound['LogP']:.2f}")
            st.write(f"**Polar Surface Area:** {selected_compound['TPSA']:.2f}")
            st.write(f"**Ring Count:** {selected_compound['nRing']}")
            st.write(f"**Follows Lipinski's Rule:** {selected_compound['FollowsLipinski']}")
            st.write(f"**Violations:** {selected_compound['Violations']}")
            st.write(f"**Atom Distribution:** {selected_compound['AtomDistribution']}")
        
        with col2:
            img = Draw.MolToImage(selected_compound['Mol'], size=(500, 500), dpi=1200)
            st.image(img, caption="2D Structure", use_column_width=True)
        
        with col3:
            radar_fig = plot_radar_normalized(selected_compound)
            st.pyplot(radar_fig)

if __name__ == "__main__":
    lip2()
