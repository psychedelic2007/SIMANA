import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdFMCS
from io import BytesIO

# Initialize session state variables
if 'uploaded_data' not in st.session_state:
    st.session_state.uploaded_data = None
if 'similarity_matrix' not in st.session_state:
    st.session_state.similarity_matrix = None
if 'smiles_list' not in st.session_state:
    st.session_state.smiles_list = None
if 'color_scheme' not in st.session_state:
    st.session_state.color_scheme = 'Blues'

def calculate_similarity_and_highlight(smiles1, smiles2):
    """Calculate Tanimoto similarity and highlight common substructure"""
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    if mol1 is None or mol2 is None:
        return None, None, None, None
    
    # Generate fingerprints
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    
    # Calculate Tanimoto similarity
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    
    # Find Maximum Common Substructure (MCS)
    mcs = rdFMCS.FindMCS([mol1, mol2])
    
    if mcs.numAtoms > 0:
        # Get the SMARTS pattern for the MCS
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        
        # Create match objects
        matches1 = mol1.GetSubstructMatch(mcs_mol)
        matches2 = mol2.GetSubstructMatch(mcs_mol)
        
        # Generate 2D depictions with highlighted substructures
        img1 = Draw.MolToImage(mol1, highlightAtoms=matches1)
        img2 = Draw.MolToImage(mol2, highlightAtoms=matches2)
        
        return similarity, img1, img2, True
    
    return similarity, None, None, False

def process_file(file_content):
    """Process uploaded file content and calculate similarity matrix"""
    smiles_list = [line.strip() for line in file_content.split('\n') if line.strip()]
    n = len(smiles_list)
    
    # Calculate similarity matrix
    similarity_matrix = np.zeros((n, n))
    mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) for mol in mols if mol is not None]
    
    for i in range(n):
        for j in range(n):
            if mols[i] is not None and mols[j] is not None:
                similarity = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                similarity_matrix[i, j] = similarity
    
    return similarity_matrix, smiles_list

def plot_heatmap(similarity_matrix, color_scheme):
    """Plot heatmap with given color scheme"""
    fig, ax = plt.subplots(figsize=(10, 8))
    mask = np.triu(np.ones_like(similarity_matrix), k=1)
    sns.heatmap(similarity_matrix, 
                mask=mask,
                cmap=color_scheme,
                vmin=0,
                vmax=1,
                annot=False,
                square=True,
                cbar_kws={'label': 'Tanimoto coefficient'})
    
    # Add labels
    num_compounds = similarity_matrix.shape[0]
    plt.xticks(range(num_compounds), [f'Compound {i+1}' for i in range(num_compounds)], rotation=45)
    plt.yticks(range(num_compounds), [f'Compound {i+1}' for i in range(num_compounds)], rotation=0)
    
    plt.title('Tanimoto Similarity Matrix')
    return fig

def handle_file_upload():
    """Handle file upload and process data"""
    uploaded_file = st.file_uploader("Upload a text file with SMILES (one per line)", type=["txt"])
    
    if uploaded_file is not None:
        # Only process the file if it hasn't been processed before
        if st.session_state.uploaded_data != uploaded_file.getvalue():
            st.session_state.uploaded_data = uploaded_file.getvalue()
            file_content = uploaded_file.getvalue().decode()
            st.session_state.similarity_matrix, st.session_state.smiles_list = process_file(file_content)
            return True
    return False

def tanimoto():
    """Main function for the Tanimoto Similarity Calculator page"""
    st.title("Molecular Tanimoto Similarity Calculator")
    
    # Color scheme options
    color_schemes = {
        'Blues': 'Blues',
        'Reds': 'Reds',
        'Greens': 'Greens',
        'Purples': 'Purples',
        'YlOrRd': 'YlOrRd',
        'YlGnBu': 'YlGnBu'
    }
    
    # Sidebar for input method selection
    input_method = st.sidebar.radio("Select Input Method", 
                                  ["Compare Two SMILES", "Upload File"])
    
    if input_method == "Compare Two SMILES":
        col1, col2 = st.columns(2)
        
        with col1:
            smiles1 = st.text_area("Enter SMILES for first molecule", height=100)
        with col2:
            smiles2 = st.text_area("Enter SMILES for second molecule", height=100)
        
        if st.button("Submit Comparison"):
            if smiles1 and smiles2:
                similarity, img1, img2, has_mcs = calculate_similarity_and_highlight(smiles1, smiles2)
                
                if similarity is not None:
                    st.write(f"Tanimoto Similarity: {similarity:.3f}")
                    
                    if has_mcs:
                        col1, col2 = st.columns(2)
                        with col1:
                            st.image(img1, caption="Molecule 1 (Common substructure highlighted)")
                        with col2:
                            st.image(img2, caption="Molecule 2 (Common substructure highlighted)")
                    else:
                        st.warning("No significant common substructure found")
                else:
                    st.error("Invalid SMILES strings provided")
    
    else:
        # File upload section
        file_processed = handle_file_upload()
        
        if st.session_state.similarity_matrix is not None:
            # Color scheme selector
            selected_color = st.selectbox(
                "Select color scheme",
                list(color_schemes.keys()),
                index=list(color_schemes.keys()).index(st.session_state.color_scheme)
            )
            st.session_state.color_scheme = selected_color
            
            # Plot heatmap
            fig = plot_heatmap(st.session_state.similarity_matrix, color_schemes[selected_color])
            st.pyplot(fig)
            
            # Add download button
            buf = BytesIO()
            plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
            st.download_button(
                label="Download Heatmap (300 DPI)",
                data=buf.getvalue(),
                file_name="tanimoto_similarity.png",
                mime="image/png"
            )
            
            # Compound comparison section
            st.subheader("Compare Individual Compounds")
            col1, col2 = st.columns(2)
            
            with col1:
                comp1 = st.selectbox("Select first compound", 
                                   [f"Compound {i+1}" for i in range(len(st.session_state.smiles_list))])
            with col2:
                comp2 = st.selectbox("Select second compound", 
                                   [f"Compound {i+1}" for i in range(len(st.session_state.smiles_list))])
            
            if st.button("Compare Selected Compounds"):
                idx1 = int(comp1.split()[-1]) - 1
                idx2 = int(comp2.split()[-1]) - 1
                
                similarity, img1, img2, has_mcs = calculate_similarity_and_highlight(
                    st.session_state.smiles_list[idx1], 
                    st.session_state.smiles_list[idx2]
                )
                
                if similarity is not None:
                    st.write(f"Tanimoto Similarity: {similarity:.3f}")
                    
                    if has_mcs:
                        col1, col2 = st.columns(2)
                        with col1:
                            st.image(img1, caption=f"{comp1} (Common substructure highlighted)")
                        with col2:
                            st.image(img2, caption=f"{comp2} (Common substructure highlighted)")
                    else:
                        st.warning("No significant common substructure found")

if __name__ == "__main__":
    tanimoto()
