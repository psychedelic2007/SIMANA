import streamlit as st
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.database.rcsb as rcsb
import biotite.structure.io.pdb as pdb
from biotite.structure.io import load_structure
from biotite.structure.io.pdb import PDBFile
from stmol import showmol
import py3Dmol
from io import StringIO
import requests

import requests

import requests
from datetime import datetime

def fetch_pdb_metadata(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        
        # Extract required information
        method = data.get('exptl', [{}])[0].get('method', 'N/A')
        resolution = data.get('rcsb_entry_info', {}).get('resolution_combined', [0])[0]
        
        # Fix R-Value Observed
        r_value = data.get('refine', [{}])[0].get('ls_R_factor_obs', 'N/A')
        if r_value == 'N/A':
            r_value = data.get('refine', [{}])[0].get('ls_R_factor_R_work', 'N/A')
        
        # Format Deposit Date
        deposit_date = data.get('rcsb_accession_info', {}).get('deposit_date', 'N/A')
        if deposit_date != 'N/A':
            deposit_date = datetime.fromisoformat(deposit_date.replace('Z', '+00:00')).strftime('%Y-%m-%d')
        
        # Correct Total Structure Weight
        total_weight = data.get('rcsb_entry_info', {}).get('molecular_weight', 'N/A')
        if total_weight != 'N/A':
            total_weight = f"{total_weight / 1000:.2f} kDa"
        
        atom_count = data.get('rcsb_entry_info', {}).get('deposited_atom_count', 'N/A')
        
        # Correct Modelled Residue Count
        modelled_residue_count = data.get('rcsb_entry_info', {}).get('deposited_model_count', 'N/A')
        
        deposited_residue_count = data.get('rcsb_entry_info', {}).get('deposited_polymer_monomer_count', 'N/A')
        unique_protein_chains = data.get('rcsb_entry_info', {}).get('polymer_entity_count_protein', 'N/A')
        unique_nucleic_acid_chains = data.get('rcsb_entry_info', {}).get('polymer_entity_count_nucleic_acid', 'N/A')
        
        return {
            'pdb_id': pdb_id,
            'description': data.get('struct', {}).get('title', 'N/A'),
            'method': method,
            'resolution': resolution,
            'r_value': r_value,
            'deposit_date': deposit_date,
            'total_weight': total_weight,
            'atom_count': atom_count,
            'modelled_residue_count': modelled_residue_count,
            'deposited_residue_count': deposited_residue_count,
            'unique_protein_chains': unique_protein_chains,
            'unique_nucleic_acid_chains': unique_nucleic_acid_chains
        }
    else:
        return None
        
def fetch_and_analyze(pdb_id, threshold_distance, dna_chain_ids, protein_chain_ids):
    # Fetch and load structure
    pdbx_file = pdbx.BinaryCIFFile.read(rcsb.fetch(pdb_id, "bcif"))
    structure = pdbx.get_structure(pdbx_file, model=1)
    
    # Separate structure into the DNA and protein chains
    dna_chain_ids = dna_chain_ids.split(',')
    protein_chain_ids = protein_chain_ids.split(',')
    dna = structure[np.isin(structure.chain_id, dna_chain_ids) & (structure.hetero == False)]
    
    proteins = []
    for chain_id in protein_chain_ids:
        proteins.append(structure[(structure.chain_id == chain_id) & (structure.hetero == False)])
    
    # Check if all protein chains are identical
    assert all(len(struc.get_residues(proteins[0])) == len(struc.get_residues(protein)) for protein in proteins)

    # Fast identification of contacts via a cell list:
    cell_list = struc.CellList(dna, cell_size=threshold_distance)

    # Sets to store the residue IDs of contact residues
    id_sets = [set() for _ in proteins]

    for protein, res_id_set in zip(proteins, id_sets):
        # For each atom in the protein chain, find all atoms in the DNA that are in contact with it
        contacts = cell_list.get_atoms(protein.coord, radius=threshold_distance)
        contact_indices = np.where((contacts != -1).any(axis=1))[0]
        contact_res_ids = protein.res_id[contact_indices]
        res_id_set.update(contact_res_ids)

    # Only consider residues that show contacts in all peptide chains
    common_ids = sorted(list(set.intersection(*id_sets)))

    # Generate output
    result = []
    for res_id in common_ids:
        res_name = proteins[0].res_name[proteins[0].res_id == res_id][0]
        result.append(res_name.capitalize() + str(res_id))

    return structure, result, common_ids

def render_mol(pdb_string, common_ids, chain_ids):
    pdbview = py3Dmol.view(width=800, height=500)
    pdbview.addModel(pdb_string, 'pdb')
    pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
    for chain_id in chain_ids:
        for res_id in common_ids:
            selection = {'chain': str(chain_id), 'resi': str(res_id)}
            pdbview.addStyle(selection, {'stick': {'color': 'red'}})
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.spin(True)
    return pdbview

# Initialize session state variables
if 'metadata' not in st.session_state:
    st.session_state.metadata = None
if 'structure' not in st.session_state:
    st.session_state.structure = None
if 'result' not in st.session_state:
    st.session_state.result = None
if 'common_ids' not in st.session_state:
    st.session_state.common_ids = None
if 'pdb_string' not in st.session_state:
    st.session_state.pdb_string = None
    
def pdi_try():
    st.title("Protein-DNA Contact Analysis")

    # Input fields
    col1, col2 = st.columns(2)
    with col1:
        pdb_id = st.text_input("Please Enter the PDB ID:", "2or1")
        threshold_distance = st.number_input("Please Give the Threshold Distance:", min_value=0.1, value=4.0)
        dna_chain_ids = st.text_input("Please Enter the DNA Chain ID(s) (comma separated):", "A,B")
    with col2:
        num_proteins = st.number_input("Please Enter the Number of Proteins:", min_value=1, max_value=10, value=2)
        protein_chain_ids = st.text_input("Please Enter the Protein Chain IDs (comma separated):", "L,R")

    # Submit button
    if st.button("Submit"):
        st.session_state.structure, st.session_state.result, st.session_state.common_ids = fetch_and_analyze(pdb_id, threshold_distance, dna_chain_ids, protein_chain_ids)
        st.session_state.metadata = fetch_pdb_metadata(pdb_id)
        
        # Convert structure to PDB string for visualization
        file = StringIO()
        pdb_file = PDBFile()
        pdb_file.set_structure(st.session_state.structure)
        pdb_file.write(file)
        st.session_state.pdb_string = file.getvalue()
        
    # Fetch and display metadata
    if st.session_state.metadata:
        st.write("***")
        col1, col2 = st.columns([2,1])
            
        with col1:
            st.subheader("General Information")
            st.write(f"**PDB ID:** {st.session_state.metadata['pdb_id']} ({st.session_state.metadata['description']})")
            st.write(f"**Method:** {st.session_state.metadata['method']}")
            st.write(f"**Resolution:** {st.session_state.metadata['resolution']} Å")
            st.write(f"**R-Value Observed:** {st.session_state.metadata['r_value']}")
            st.write(f"**Deposit Date:** {st.session_state.metadata['deposit_date']}")
            
        with col2:
            st.subheader("Macromolecule Content")
            st.write(f"**Total Structure Weight:** {st.session_state.metadata['total_weight']} Da")
            st.write(f"**Atom Count:** {st.session_state.metadata['atom_count']}")
            st.write(f"**Modelled Residue Count:** {st.session_state.metadata['modelled_residue_count']}")
            st.write(f"**Deposited Residue Count:** {st.session_state.metadata['deposited_residue_count']}")
            st.write(f"**Unique protein chains:** {st.session_state.metadata['unique_protein_chains']}")
            st.write(f"**Unique nucleic acid chains:** {st.session_state.metadata['unique_nucleic_acid_chains']}")

    st.write("***")
    
    # Display results if available
    if st.session_state.result is not None:
        col1, col2 = st.columns([1,2])

        with col1:
            st.subheader("Residues in contact with DNA:")
            table_data = [["Residue", "Residue"]]
            for i in range(0, len(st.session_state.result), 2):
                row = [st.session_state.result[i] if i < len(st.session_state.result) else ""]
                if i+1 < len(st.session_state.result):
                    row.append(st.session_state.result[i+1])
                else:
                    row.append("")
                table_data.append(row)
            
            st.table(table_data)
            
            # Viewer Settings
            st.subheader("Viewer Settings")
            displayed_residues = st.multiselect("Displayed Residues", st.session_state.result, default=st.session_state.result[:3])
            label_residues = st.checkbox("Label Residues", value=True)
            color_scheme = st.radio("Color Scheme", ["Standard", "Amino Acid"])
            cartoon_style = st.radio("Cartoon Style", ["Ribbon", "Trace"])
            cartoon_transparency = st.slider("Cartoon Transparency", 0.0, 1.0, 1.0)
            surface_transparency = st.slider("Surface Transparency", 0.0, 1.0, 0.0)
            rotate_structure = st.checkbox("Rotate Structure", value=False)

        with col2:
            st.write("PDB Structure with Highlighted Residues:")
            try:
                view = py3Dmol.view(width=800, height=500)
                view.addModel(st.session_state.pdb_string, 'pdb')
                
                # Set style for protein
                view.setStyle({'cartoon': {'style': cartoon_style.lower(), 'color': 'spectrum', 'opacity': cartoon_transparency}})
                
                # Set style for DNA
                view.setStyle({'hetflag': False, 'chains': dna_chain_ids.split(',')}, {'cartoon': {'style': 'nucleic', 'color': 'blue'}})
                
                # Highlight contact residues
                for chain_id in protein_chain_ids.split(','):
                    for res_id in st.session_state.common_ids:
                        selection = {'chain': str(chain_id), 'resi': str(res_id)}
                        view.addStyle(selection, {'stick': {'color': 'red', 'radius': 0.7}})
                        
                if surface_transparency > 0:
                    view.addSurface(py3Dmol.VDW, {'opacity': surface_transparency})
                if label_residues:
                    for res in displayed_residues:
                        view.addLabel(res, {'font': 'Arial', 'fontSize': 12, 'fontColor': 'black', 'backgroundOpacity': 0.2})
                view.setBackgroundColor('white')
                view.zoomTo()
                if rotate_structure:
                    view.spin(True)
                else:
                    view.spin(False)
        
                # Debug print statements
                print("Attempting to show molecule")
                showmol(view, height=500, width=800)
                print("Molecule display attempted")
            except Exception as e:
                st.error(f"Error rendering molecule: {str(e)}")

if __name__ == "__main__":
    pdi_try()
