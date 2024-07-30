import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import PDBParser
import csv
import io
import base64

def calculate_b_factors(pdb_content):
    parser = PDBParser(QUIET=True)
    pdb_io = io.StringIO(pdb_content)
    structure = parser.get_structure('structure', pdb_io)
    
    residues_b_factors = {}
    all_b_factors = []
    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = residue.get_id()
                b_factors = [atom.get_bfactor() for atom in residue]
                residues_b_factors[res_id] = b_factors
                all_b_factors.extend(b_factors)
    
    return residues_b_factors, all_b_factors

def plot_b_factors(all_residues_b_factors, labels, show_std_dev, customizations_curve):
    fig, ax = plt.subplots(figsize=(12, 6))
    
    for label, residues_b_factors in zip(labels, all_residues_b_factors):
        residue_ids = []
        b_factors_mean = []
        b_factors_std = []

        for res_id, b_factors in residues_b_factors.items():
            residue_ids.append(res_id[1])
            b_factors_mean.append(np.mean(b_factors))
            b_factors_std.append(np.std(b_factors))

        ax.plot(residue_ids, b_factors_mean, label=label, linewidth=customizations_curve['linewidth'], color=customizations_curve['colors'][0])
        
        if show_std_dev:
            ax.fill_between(residue_ids, 
                            [m - s for m, s in zip(b_factors_mean, b_factors_std)], 
                            [m + s for m, s in zip(b_factors_mean, b_factors_std)], 
                            alpha=0.2)

    ax.set_xlabel(customizations_curve['x_label'], fontsize=customizations_curve['x_label_size'])
    ax.set_ylabel(customizations_curve['y_label'], fontsize=customizations_curve['y_label_size'])
    ax.set_xticks(np.arange(min(residue_ids), max(residue_ids)+1, customizations_curve['x_tick_gap']))
    ax.set_yticks(np.arange(min(b_factors_mean), max(b_factors_mean)+1, customizations_curve['y_tick_gap']))
    ax.tick_params(axis='both', which='major', labelsize=customizations_curve['tick_size'])
    ax.tick_params(axis='x', rotation=customizations_curve['x_tick_rotation'])
    if customizations_curve['x_min'] is not None and customizations_curve['x_max'] is not None:
        ax.set_xlim([customizations_curve['x_min'], customizations_curve['x_max']])
    if customizations_curve['y_min'] is not None and customizations_curve['y_max'] is not None:
        ax.set_ylim([customizations_curve['y_min'], customizations_curve['y_max']])
    ax.set_title('B-factor per Residue')
    ax.legend()
    return fig

def plot_b_factor_distribution(all_b_factors_list, labels, customizations_dist):
    fig, ax = plt.subplots(figsize=(12, 6))
    
    for label, all_b_factors in zip(labels, all_b_factors_list):
        sns.kdeplot(all_b_factors, label=label, ax=ax, color=customizations_dist['colors'][0], alpha=customizations_dist['alpha'])

    ax.set_xlabel(customizations_dist['x_label'], fontsize=customizations_dist['x_label_size'])
    ax.set_ylabel(customizations_dist['y_label'], fontsize=customizations_dist['y_label_size'])
    ax.set_xticks(np.arange(min(all_b_factors_list[0]), max(all_b_factors_list[0])+1, customizations_dist['x_tick_gap']))
    ax.set_yticks(np.arange(0, len(all_b_factors_list[0]), customizations_dist['y_tick_gap']))
    ax.tick_params(axis='both', which='major', labelsize=customizations_dist['tick_size'])
    ax.tick_params(axis='x', rotation=customizations_dist['x_tick_rotation'])
    if customizations_dist['x_min'] is not None and customizations_dist['x_max'] is not None:
        ax.set_xlim([customizations_dist['x_min'], customizations_dist['x_max']])
    if customizations_dist['y_min'] is not None and customizations_dist['y_max'] is not None:
        ax.set_ylim([customizations_dist['y_min'], customizations_dist['y_max']])
    ax.set_title('Distribution of B-factors')
    ax.legend()
    return fig

def get_csv_download_link(all_residues_b_factors, labels):
    csv_buffer = io.StringIO()
    writer = csv.writer(csv_buffer)
    
    writer.writerow(['Label', 'Residue', 'Mean B-factor', 'Std B-factor'])
    for label, residues_b_factors in zip(labels, all_residues_b_factors):
        for res_id, b_factors in residues_b_factors.items():
            mean_b_factor = np.mean(b_factors)
            std_b_factor = np.std(b_factors)
            writer.writerow([label, res_id[1], mean_b_factor, std_b_factor])
    
    csv_string = csv_buffer.getvalue()
    b64 = base64.b64encode(csv_string.encode()).decode()
    return f'<a href="data:file/csv;base64,{b64}" download="b_factors.csv">Download CSV File</a>'

def bfactor():
    st.title("B-factor Calculator")

    if 'uploaded_files' not in st.session_state:
        st.session_state['uploaded_files'] = []
        st.session_state['labels'] = []
        st.session_state['uploaded_contents'] = []
        st.session_state['uploaded'] = False

    if 'num_files' not in st.session_state:
        st.session_state['num_files'] = 1

    if 'customizations_curve' not in st.session_state:
        st.session_state.customizations_curve = {
            'x_label': "Residue Number",
            'y_label': "B-factor Mean",
            'x_label_size': 12,
            'y_label_size': 12,
            'tick_size': 10,
            'x_tick_gap': 10,
            'y_tick_gap': 0.1,
            'linewidth': 1.5,
            'x_tick_rotation': 0,
            'colors': ["#000000"],
            'x_min': None,
            'x_max': None,
            'y_min': None,
            'y_max': None,
        }

    if 'customizations_dist' not in st.session_state:
        st.session_state.customizations_dist = {
            'x_label': "B-factor",
            'y_label': "Density",
            'x_label_size': 12,
            'y_label_size': 12,
            'tick_size': 10,
            'x_tick_gap': 0.1,
            'y_tick_gap': 0.01,
            'x_tick_rotation': 0,
            'colors': ["#000000"],
            'alpha': 0.5,
            'x_min': None,
            'x_max': None,
            'y_min': None,
            'y_max': None,
        }

    num_files = st.number_input("Number of files to upload", min_value=1, max_value=10, value=st.session_state['num_files'], key='num_files')

    if not st.session_state['uploaded']:
        for i in range(num_files):
            col1, col2 = st.columns(2)
            with col1:
                file = st.file_uploader(f"Choose PDB file {i+1}", type=['pdb'], key=f"file_{i}")
            with col2:
                label = st.text_input(f"Label for file {i+1}", value=f"File {i+1}", key=f"label_{i}")

            if file:
                pdb_content = file.getvalue().decode('utf-8')
                if file.name not in st.session_state['uploaded_files']:
                    st.session_state['uploaded_files'].append(file.name)
                    st.session_state['labels'].append(label)
                    st.session_state['uploaded_contents'].append(pdb_content)

    if st.button("Submit"):
        if len(st.session_state['uploaded_files']) == num_files:
            st.session_state['uploaded'] = True
        else:
            st.warning(f"Please upload all {num_files} PDB files before submitting.")

    if st.session_state['uploaded']:
        show_std_dev = st.checkbox("Show standard deviation", value=False, key='show_std_dev')

        all_residues_b_factors = []
        all_b_factors_list = []

        for pdb_content in st.session_state['uploaded_contents']:
            try:
                residues_b_factors, all_b_factors = calculate_b_factors(pdb_content)
                all_residues_b_factors.append(residues_b_factors)
                all_b_factors_list.append(all_b_factors)
            except Exception as e:
                st.error(f"An error occurred while processing a PDB file: {str(e)}")

        if all_residues_b_factors:
            st.write("B-factor per Residue Plot")
            fig1 = plot_b_factors(all_residues_b_factors, st.session_state['labels'], st.session_state['show_std_dev'], st.session_state.customizations_curve)
            st.pyplot(fig1)

            st.write("Distribution of B-factors Plot")
            fig2 = plot_b_factor_distribution(all_b_factors_list, st.session_state['labels'], st.session_state.customizations_dist)
            st.pyplot(fig2)

            st.markdown(get_csv_download_link(all_residues_b_factors, st.session_state['labels']), unsafe_allow_html=True)

            with st.expander("Curve Customizations"):
                col1, col2 = st.columns(2)
                with col1:
                    st.session_state.customizations_curve['x_label'] = st.text_input("X-axis label", value=st.session_state.customizations_curve['x_label'])
                    st.session_state.customizations_curve['y_label'] = st.text_input("Y-axis label", value=st.session_state.customizations_curve['y_label'])
                    st.session_state.customizations_curve['x_label_size'] = st.number_input("X-axis label size", value=st.session_state.customizations_curve['x_label_size'])
                    st.session_state.customizations_curve['y_label_size'] = st.number_input("Y-axis label size", value=st.session_state.customizations_curve['y_label_size'])
                    st.session_state.customizations_curve['tick_size'] = st.number_input("Tick size", value=st.session_state.customizations_curve['tick_size'])
                    st.session_state.customizations_curve['x_tick_gap'] = st.number_input("X-axis tick gap", value=st.session_state.customizations_curve['x_tick_gap'])
                with col2:
                    st.session_state.customizations_curve['y_tick_gap'] = st.number_input("Y-axis tick gap", value=st.session_state.customizations_curve['y_tick_gap'])
                    st.session_state.customizations_curve['linewidth'] = st.number_input("Line width", value=st.session_state.customizations_curve['linewidth'])
                    st.session_state.customizations_curve['x_tick_rotation'] = st.number_input("X-axis tick rotation", value=st.session_state.customizations_curve['x_tick_rotation'])
                    st.session_state.customizations_curve['x_min'] = st.number_input("X-axis min", value=st.session_state.customizations_curve['x_min'])
                    st.session_state.customizations_curve['x_max'] = st.number_input("X-axis max", value=st.session_state.customizations_curve['x_max'])
                    st.session_state.customizations_curve['y_min'] = st.number_input("Y-axis min", value=st.session_state.customizations_curve['y_min'])
                    st.session_state.customizations_curve['y_max'] = st.number_input("Y-axis max", value=st.session_state.customizations_curve['y_max'])

            with st.expander("Distribution Customizations"):
                col1, col2 = st.columns(2)
                with col1:
                    st.session_state.customizations_dist['x_label'] = st.text_input("X-axis label", value=st.session_state.customizations_dist['x_label'])
                    st.session_state.customizations_dist['y_label'] = st.text_input("Y-axis label", value=st.session_state.customizations_dist['y_label'])
                    st.session_state.customizations_dist['x_label_size'] = st.number_input("X-axis label size", value=st.session_state.customizations_dist['x_label_size'])
                    st.session_state.customizations_dist['y_label_size'] = st.number_input("Y-axis label size", value=st.session_state.customizations_dist['y_label_size'])
                    st.session_state.customizations_dist['tick_size'] = st.number_input("Tick size", value=st.session_state.customizations_dist['tick_size'])
                    st.session_state.customizations_dist['x_tick_gap'] = st.number_input("X-axis tick gap", value=st.session_state.customizations_dist['x_tick_gap'])
                with col2:
                    st.session_state.customizations_dist['y_tick_gap'] = st.number_input("Y-axis tick gap", value=st.session_state.customizations_dist['y_tick_gap'])
                    st.session_state.customizations_dist['x_tick_rotation'] = st.number_input("X-axis tick rotation", value=st.session_state.customizations_dist['x_tick_rotation'])
                    st.session_state.customizations_dist['x_min'] = st.number_input("X-axis min", value=st.session_state.customizations_dist['x_min'])
                    st.session_state.customizations_dist['x_max'] = st.number_input("X-axis max", value=st.session_state.customizations_dist['x_max'])
                    st.session_state.customizations_dist['y_min'] = st.number_input("Y-axis min", value=st.session_state.customizations_dist['y_min'])
                    st.session_state.customizations_dist['y_max'] = st.number_input("Y-axis max", value=st.session_state.customizations_dist['y_max'])
                    st.session_state.customizations_dist['alpha'] = st.number_input("Alpha (transparency)", value=st.session_state.customizations_dist['alpha'])

if __name__ == "__main__":
    bfactor()
