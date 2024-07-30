import streamlit as st
from matplotlib import pyplot as plt
import numpy as np
import tempfile
import os
import MDAnalysis as mda
from MDAnalysis.analysis.align import alignto

# Function to calculate DCCM
def calculate_dccm(uploaded_files_pdb, uploaded_files_xtc):
    dccms = []

    for pdb_path, xtc_path in zip(uploaded_files_pdb, uploaded_files_xtc):
        u = mda.Universe(pdb_path, xtc_path)
        ca_atoms = u.select_atoms('name CA')
        alignto(u, u, select='name CA')

        positions = np.zeros((len(u.trajectory), len(ca_atoms), 3))

        for i, ts in enumerate(u.trajectory):
            positions[i] = ca_atoms.positions

        mean_positions = positions.mean(axis=0)
        fluctuations = positions - mean_positions
        covariance_matrix = np.tensordot(fluctuations, fluctuations, axes=((0, 2), (0, 2)))
        dccm = np.corrcoef(covariance_matrix)
        dccms.append(dccm)

    average_dccm = np.mean(dccms, axis=0)

    # Clean up temporary files
    for pdb_path, xtc_path in zip(uploaded_files_pdb, uploaded_files_xtc):
        os.remove(pdb_path)
        os.remove(xtc_path)

    return average_dccm

# Function to plot DCCM
def plot_dccm(dccm, customizations):
    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(dccm, cmap=customizations['cmap'], vmin=customizations['vmin'], vmax=customizations['vmax'])
    ax.set_xlabel(customizations['xlabel'], fontsize=15)
    ax.set_ylabel(customizations['ylabel'], fontsize=15)
    ax.set_xticks(np.arange(customizations['xlim_min'], customizations['xlim_max'] + 1, customizations['xticks_gap']))
    ax.set_yticks(np.arange(customizations['ylim_min'], customizations['ylim_max'] + 1, customizations['yticks_gap']))
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.invert_yaxis()

    # Add color bar with label
    plt.colorbar(im, ax=ax, label='Correlation Coefficient')

    return fig

# Main Streamlit app
def dccm():
    st.title('DCCM Analysis')

    if 'dccm_data' not in st.session_state:
        st.session_state.dccm_data = None

    if 'customizations' not in st.session_state:
        st.session_state.customizations = {
            'cmap': 'viridis',
            'vmin': -1.0,
            'vmax': 1.0,
            'xlabel': 'Residue index',
            'ylabel': 'Residue index',
            'xticks_gap': 1,
            'yticks_gap': 1,
            'xlim_min': 0,
            'xlim_max': 100,
            'ylim_min': 0,
            'ylim_max': 100
        }

    uploaded_files_pdb = []
    uploaded_files_xtc = []

    # File upload section
    st.header('Upload Files')
    num_files = st.number_input('Enter number of files to upload', min_value=1, max_value=10, value=1)

    for i in range(num_files):
        pdb_file = st.file_uploader(f'Upload PDB file {i+1}', type=['pdb'])
        xtc_file = st.file_uploader(f'Upload XTC file {i+1}', type=['xtc'])

        if pdb_file and xtc_file:
            pdb_temp = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
            xtc_temp = tempfile.NamedTemporaryFile(delete=False, suffix='.xtc')

            pdb_file.seek(0)
            pdb_temp.write(pdb_file.read())
            xtc_file.seek(0)
            xtc_temp.write(xtc_file.read())

            pdb_temp.close()
            xtc_temp.close()

            uploaded_files_pdb.append(pdb_temp.name)
            uploaded_files_xtc.append(xtc_temp.name)

    # Calculate DCCM when button is clicked
    if st.button('Calculate DCCM') and uploaded_files_pdb and uploaded_files_xtc:
        st.session_state.dccm_data = calculate_dccm(uploaded_files_pdb, uploaded_files_xtc)

    # Customization options in sidebar
    st.sidebar.title('Customization Options')

    st.session_state.customizations['cmap'] = st.sidebar.selectbox('Select colormap', ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
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

    # Plot DCCM if data is available
    if st.session_state.dccm_data is not None:
        fig = plot_dccm(st.session_state.dccm_data, st.session_state.customizations)
        st.pyplot(fig)

if __name__ == '__main__':
    dccm()
