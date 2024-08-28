import streamlit as st
from matplotlib import pyplot as plt
import numpy as np
import tempfile
import os
from ramachandraw.utils import plot
import io

def generate_ramachandran_plot(pdb_data, cmap, alpha):
    # Create a temporary file to save the uploaded PDB data
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_file:
        temp_file.write(pdb_data)
        temp_file_path = temp_file.name

    # Create a temporary file for the plot
    with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as plot_file:
        plot(
            temp_file_path,
            cmap=cmap,
            alpha=alpha,
            dpi=300,
            save=True,
            show=False,
            filename=plot_file.name
        )
        plot_file_path = plot_file.name

    # Read the plot image
    with open(plot_file_path, "rb") as f:
        image_data = f.read()

    # Clean up the temporary files
    os.remove(temp_file_path)
    os.remove(plot_file_path)

    return image_data

# Main Streamlit app
def rama():
    st.title('Ramachandran Plot Generator')

    if 'ramachandran_data' not in st.session_state:
        st.session_state.ramachandran_data = None

    if 'customizations' not in st.session_state:
        st.session_state.customizations = {
            'cmap': 'viridis',
            'alpha': 0.75
        }

    uploaded_file = st.file_uploader('Upload a PDB file', type=['pdb'])
    if uploaded_file is not None:
        st.session_state.ramachandran_data = uploaded_file.getvalue()

    # Customization options in sidebar
    st.sidebar.title('Customization Options')

    st.session_state.customizations['cmap'] = st.sidebar.selectbox(
        'Select colormap',
        ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
         'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn', 'binary',
         'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink', 'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia', 'hot',
         'afmhot', 'gist_heat', 'copper', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 'Spectral',
         'coolwarm', 'bwr', 'seismic', 'twilight', 'twilight_shifted', 'hsv', 'Pastel1', 'Pastel2', 'Paired', 'Accent',
         'Dark2', 'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b', 'tab20c', 'flag', 'prism', 'ocean', 'gist_earth',
         'terrain', 'gist_stern', 'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'gist_rainbow', 'rainbow', 'jet',
         'turbo', 'nipy_spectral', 'gist_ncar']
    )
    st.session_state.customizations['alpha'] = st.sidebar.slider(
        'Set transparency (alpha)', min_value=0.0, max_value=1.0, value=st.session_state.customizations['alpha'], step=0.01
    )

    # Calculate and plot Ramachandran plot
    if st.session_state.ramachandran_data is not None:
        image_data = generate_ramachandran_plot(
            st.session_state.ramachandran_data,
            st.session_state.customizations['cmap'],
            st.session_state.customizations['alpha']
        )
        st.image(image_data, caption='Ramachandran Plot', use_column_width=True)
        
        # Provide download button for the plot
        st.download_button(
            label="Download Ramachandran Plot",
            data=image_data,
            file_name="ramachandran_plot.png",
            mime="image/png"
        )

if __name__ == '__main__':
    rama()
