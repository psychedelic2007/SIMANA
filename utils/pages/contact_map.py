import streamlit as st
from matplotlib import pyplot as plt
import numpy as np
import tempfile
import os
import MDAnalysis as mda
from itertools import combinations
from io import BytesIO

# Function to calculate the contact map
def calculate_contact_map(pdb_file_path, cutoff=8.0):
    u = mda.Universe(pdb_file_path)
    protein = u.select_atoms('protein')
    ca_atoms = protein.select_atoms('name CA')
    num_residues = len(ca_atoms)
    contact_map = np.zeros((num_residues, num_residues))

    for i, j in combinations(range(num_residues), 2):
        distance = np.linalg.norm(ca_atoms.positions[i] - ca_atoms.positions[j])
        if distance < cutoff:
            contact_map[i, j] = 1
            contact_map[j, i] = 1
    
    return contact_map, ca_atoms

# Function to plot the contact map
def plot_contact_map(contact_map, ca_atoms, customizations):
    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(contact_map, cmap=customizations['cmap'], vmin=customizations['vmin'], vmax=customizations['vmax'])
    ax.set_xlabel(customizations['xlabel'], fontsize=customizations['label_fontsize'])
    ax.set_ylabel(customizations['ylabel'], fontsize=customizations['label_fontsize'])
    ax.set_xticks(np.arange(customizations['xlim_min'], customizations['xlim_max'] + 1, customizations['xticks_gap']))
    ax.set_yticks(np.arange(customizations['ylim_min'], customizations['ylim_max'] + 1, customizations['yticks_gap']))
    ax.tick_params(axis='both', which='major', labelsize=customizations['tick_labelsize'])
    ax.invert_yaxis()

    plt.colorbar(im, ax=ax, label='Contact')

    return fig

# Function to convert plot to a downloadable format
def save_plot_to_bytesio(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png")
    buf.seek(0)
    return buf

# Main Streamlit app
def contact_map():
    st.title('Contact Map Generator')

    if 'contact_map_data' not in st.session_state:
        st.session_state.contact_map_data = None

    if 'customizations' not in st.session_state:
        st.session_state.customizations = {
            'cmap': 'viridis',
            'vmin': 0.0,
            'vmax': 1.0,
            'xlabel': 'Residue Index',
            'ylabel': 'Residue Index',
            'label_fontsize': 15,
            'tick_labelsize': 12,
            'xticks_gap': 1,
            'yticks_gap': 1,
            'xlim_min': 0,
            'xlim_max': 100,
            'ylim_min': 0,
            'ylim_max': 100
        }

    # File upload section
    st.header('Upload Files')
    num_files = st.number_input('Enter number of files to upload', min_value=1, max_value=10, value=1)

    uploaded_files_pdb = []
    
    for i in range(num_files):
        pdb_file = st.file_uploader(f'Upload PDB file {i+1}', type=['pdb'])

        if pdb_file:
            pdb_temp = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
            pdb_file.seek(0)
            pdb_temp.write(pdb_file.read())
            pdb_temp.close()
            uploaded_files_pdb.append(pdb_temp.name)

    cutoff_distance = st.sidebar.number_input("Cutoff Distance (in Å)", min_value=0.0, value=8.0, step=0.1)

    # Customization options in sidebar
    st.sidebar.title('Customization Options')

    st.session_state.customizations['cmap'] = st.sidebar.selectbox('Select Colormap', ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn', 'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
                      'pink', 'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
                      'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic', 'twilight', 'twilight_shifted', 'hsv', 'Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2',
                      'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b', 'tab20c', 'flag', 'prism', 'ocean', 'gist_earth', 'terrain',
                      'gist_stern', 'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'gist_rainbow', 'rainbow', 'jet', 'turbo', 'nipy_spectral', 'gist_ncar'])
    st.session_state.customizations['vmin'] = st.sidebar.number_input('Minimum value for colormap', value=st.session_state.customizations['vmin'])
    st.session_state.customizations['vmax'] = st.sidebar.number_input('Maximum value for colormap', value=st.session_state.customizations['vmax'])
    st.session_state.customizations['xlabel'] = st.sidebar.text_input('X-axis label', value=st.session_state.customizations['xlabel'])
    st.session_state.customizations['ylabel'] = st.sidebar.text_input('Y-axis label', value=st.session_state.customizations['ylabel'])
    st.session_state.customizations['xticks_gap'] = st.sidebar.number_input('X-axis tick gap', value=st.session_state.customizations['xticks_gap'])
    st.session_state.customizations['yticks_gap'] = st.sidebar.number_input('Y-axis tick gap', value=st.session_state.customizations['yticks_gap'])
    st.session_state.customizations['xlim_min'] = st.sidebar.number_input('X-axis minimum range', value=st.session_state.customizations['xlim_min'])
    st.session_state.customizations['xlim_max'] = st.sidebar.number_input('X-axis maximum range', value=st.session_state.customizations['xlim_max'])
    st.session_state.customizations['ylim_min'] = st.sidebar.number_input('Y-axis minimum range', value=st.session_state.customizations['ylim_min'])
    st.session_state.customizations['ylim_max'] = st.sidebar.number_input('Y-axis maximum range', value=st.session_state.customizations['ylim_max'])

    # Calculate contact map when button is clicked
    if st.button('Generate Contact Maps') and uploaded_files_pdb:
        st.session_state.contact_map_data = []
        for pdb_path in uploaded_files_pdb:
            contact_map_data, ca_atoms = calculate_contact_map(pdb_path, cutoff=cutoff_distance)
            st.session_state.contact_map_data.append((contact_map_data, ca_atoms))

    # Plot contact maps if data is available and provide download option
    if st.session_state.contact_map_data:
        for i, (contact_map_data, ca_atoms) in enumerate(st.session_state.contact_map_data):
            fig = plot_contact_map(contact_map_data, ca_atoms, st.session_state.customizations)
            st.pyplot(fig)
            
            # Save plot to bytes buffer
            img_data = save_plot_to_bytesio(fig)
            
            # Provide download button
            st.download_button(
                label=f"Download Contact Map {i+1} as PNG",
                data=img_data,
                file_name=f"contact_map_{i+1}.png",
                mime="image/png"
            )

if __name__ == '__main__':
    contact_map()
