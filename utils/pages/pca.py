import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
import base64
import io
import seaborn as sns
from sklearn.decomposition import PCA

def parse_xvg(file_path):
    x_data = []
    y_data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('@') or line.startswith('#'):
                continue
            columns = line.strip().split()
            x_data.append(float(columns[0]))
            y_data.append(float(columns[1]))
    return x_data, y_data

def calculate_pca(uploaded_files, customizations):
    all_data = []
    for uploaded_file in uploaded_files:
        with open("temp.xvg", "wb") as f:
            f.write(uploaded_file.getbuffer())
        
        x, y = parse_xvg("temp.xvg")
        all_data.extend(list(zip(x, y)))
    
    all_data = np.array(all_data)
    
    pca = PCA(n_components=2)
    pca_data = pca.fit_transform(all_data)
    
    return pca_data

def plot_pca(pca_data, customizations):
    fig, ax = plt.subplots(figsize=(15, 12))
    scatter = ax.scatter(pca_data[:, 0], pca_data[:, 1], c=pca_data[:, 0], cmap=customizations['cmap'], vmin=customizations['vmin'], vmax=customizations['vmax'])
    fig.colorbar(scatter, ax=ax)
    
    ax.set_xlabel(customizations['x_label'], fontsize=customizations['x_label_size'])
    ax.set_ylabel(customizations['y_label'], fontsize=customizations['y_label_size'])
    ax.set_title('Principal Component Analysis')
    
    ax.tick_params(axis='x', labelsize=customizations['tick_size'], rotation=customizations['x_tick_rotation'])
    ax.tick_params(axis='y', labelsize=customizations['tick_size'])
    ax.xaxis.set_major_locator(plt.MultipleLocator(customizations['x_tick_gap']))
    ax.yaxis.set_major_locator(plt.MultipleLocator(customizations['y_tick_gap']))
    ax.set_xlim(customizations['x_min'], customizations['x_max'])
    ax.set_ylim(customizations['y_min'], customizations['y_max'])
    
    return fig

def get_download_link(fig, plot_name):
    format_options = ['PDF', 'PNG', 'SVG', 'EPS', 'JPG']
    format_selected = st.selectbox(f'Select download format for {plot_name}:', format_options, key=f"{plot_name}_format")
    
    if format_selected:
        buf = io.BytesIO()
        fig.savefig(buf, format=format_selected.lower())
        buf.seek(0)
        b64 = base64.b64encode(buf.read()).decode()
        href = f'<a href="data:application/{format_selected.lower()};base64,{b64}" download="{plot_name}.{format_selected.lower()}">Download</a>'
        return href
    return None

def pca():
    st.title("Principal Component Analysis")

    if 'customizations' not in st.session_state:
        st.session_state.customizations = {
            'x_label': "PC1",
            'y_label': "PC2",
            'x_label_size': 12,
            'y_label_size': 12,
            'tick_size': 10,
            'x_tick_gap': 1,
            'y_tick_gap': 1,
            'x_tick_rotation': 0,
            'x_min': None,
            'x_max': None,
            'y_min': None,
            'y_max': None,
            'cmap': 'viridis',
            'vmin': None,
            'vmax': None,
        }

    if 'data' not in st.session_state:
        st.session_state.data = None

    num_files = st.number_input('Enter number of files to upload', min_value=1, max_value=10, value=1, key='num_files')
    
    uploaded_files = []
    for i in range(num_files):
        uploaded_file = st.file_uploader(f"Upload .xvg file {i+1}", type=['xvg'], key=f'file_{i}')
        if uploaded_file:
            uploaded_files.append(uploaded_file)

    if st.button('Submit') and uploaded_files:
        pca_data = calculate_pca(uploaded_files, st.session_state.customizations)
        st.session_state.data = pca_data

    st.sidebar.header("PCA Plot Customization")
    st.session_state.customizations['x_label'] = st.sidebar.text_input("X-axis Label", value=st.session_state.customizations['x_label'])
    st.session_state.customizations['y_label'] = st.sidebar.text_input("Y-axis Label", value=st.session_state.customizations['y_label'])
    st.session_state.customizations['x_label_size'] = st.sidebar.number_input("X-axis Label Size", value=st.session_state.customizations['x_label_size'])
    st.session_state.customizations['y_label_size'] = st.sidebar.number_input("Y-axis Label Size", value=st.session_state.customizations['y_label_size'])
    st.session_state.customizations['tick_size'] = st.sidebar.number_input("Tick Size", value=st.session_state.customizations['tick_size'])
    st.session_state.customizations['x_tick_gap'] = st.sidebar.number_input("X-axis Tick Gap", value=st.session_state.customizations['x_tick_gap'])
    st.session_state.customizations['y_tick_gap'] = st.sidebar.number_input("Y-axis Tick Gap", value=st.session_state.customizations['y_tick_gap'])
    st.session_state.customizations['x_tick_rotation'] = st.sidebar.slider("X-axis Tick Rotation", min_value=0, max_value=90, value=st.session_state.customizations['x_tick_rotation'])
    st.session_state.customizations['x_min'] = st.sidebar.number_input("X-axis Min", value=st.session_state.customizations['x_min'])
    st.session_state.customizations['x_max'] = st.sidebar.number_input("X-axis Max", value=st.session_state.customizations['x_max'])
    st.session_state.customizations['y_min'] = st.sidebar.number_input("Y-axis Min", value=st.session_state.customizations['y_min'])
    st.session_state.customizations['y_max'] = st.sidebar.number_input("Y-axis Max", value=st.session_state.customizations['y_max'])
    cmap_options = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn', 'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
                      'pink', 'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
                      'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic', 'twilight', 'twilight_shifted', 'hsv', 'Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2',
                      'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b', 'tab20c', 'flag', 'prism', 'ocean', 'gist_earth', 'terrain',
                      'gist_stern', 'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'gist_rainbow', 'rainbow', 'jet', 'turbo', 'nipy_spectral', 'gist_ncar']
    st.session_state.customizations['cmap'] = st.sidebar.selectbox("Color Map", cmap_options, index=cmap_options.index(st.session_state.customizations['cmap']))
    st.session_state.customizations['vmin'] = st.sidebar.number_input("Color Map Min", value=st.session_state.customizations['vmin'])
    st.session_state.customizations['vmax'] = st.sidebar.number_input("Color Map Max", value=st.session_state.customizations['vmax'])

    if st.session_state.data is not None and len(st.session_state.data) > 0:
        pca_data = st.session_state.data
        
        fig_pca = plot_pca(pca_data, st.session_state.customizations)
        st.pyplot(fig_pca)
        
        download_link_pca = get_download_link(fig_pca, "pca_plot")
        if download_link_pca:
            st.markdown(download_link_pca, unsafe_allow_html=True)

if __name__ == '__main__':
    pca()
