import streamlit as st
import base64

st.set_page_config(layout="wide", page_title="SIMANA", page_icon="utils/images/logo.svg")  # Set the page layout to wide

def load_image(image_path):
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode()

# Load images and convert to base64
satyam_image = load_image("utils/images/sat.jpg")
arshad_image = load_image("utils/images/arshad.jpg")

def home_page():
    # Create two columns
    col1, col2 = st.columns([1, 2])  # Adjust the proportions as needed

    # Left column content (GIF)
    with col1:
        st.image("utils/images/protein.gif", use_column_width=True)  # Set the width as needed
        
    # Right column content
    with col2:
        st.markdown("""
        <h1>Welcome to <span style="color:cyan;">SIM</span>ulation <span style="color:cyan;">ANA</span>lysis</h1>
        <p>Navigate the Depths of Data with Ease</p>
        <p>SIMANA: Where Complex Simulations Find Simplicity in Analysis</p>
        """, unsafe_allow_html=True)
    
    st.write('''***''')
    
    st.title("What SIMANA Does?")
    
    col1, col2 = st.columns([2, 2])
    
    with col1:
    	st.markdown("""
    	<h1>Data Analytics & Plotting</h1>
    	<p>✓ Enabling Seamless Plotting and Customized Analysis of MD Simulation Data, SIMANA empowers users to effortlessly visualize their data, tailoring plots to their exact specifications.</p>
	<p>✓ Leveraging the power of Python, alongside essential libraries like Matplotlib, Seaborn, and Numpy, we simplify the graphical analysis process, ensuring every insight is within reach.</p>
    	""", unsafe_allow_html=True)
    
    with col2:
    	st.image("utils/images/Data_Analysis.svg", use_column_width=True)
    	
    st.write('''***''')
    
    col1, col2 = st.columns([2, 2])
    
    with col1:
    	st.image("utils/images/working.gif", use_column_width=True)
    
    with col2:
    	st.markdown("""
    	<h1>Key Features</h1>
    	<p>✓ Seamless visualization of MD Simulation Data</p>
	<p>✓ Customizable plotting options</p>
	<p>✓ Dynamic plotting capabilities tailored to user specifications</p>
	<p>✓ Interactive plotting tools for in-depth exploration</p>
	<p>✓ Intuitive user interface for effortless navigation</p>
	<p>✓ Export options for saving plots in various formats</p>
	<p>✓ Built-in statistical analysis tools for quantifying simulation outcomes</p>
    	""", unsafe_allow_html=True)
    	
    st.write('''***''')
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
    	st.markdown("""
   	<h1>Tutorial</h1>
    	<p>✓ Introduction to SIMANA: Getting Started</p>
    	<p>✓ Uploading and Managing Simulation Data</p>
    	<p>✓ Plotting Basics: Creating Your First Plot</p>
    	<p>✓ Customizing Plot Styles and Options</p>
    	<p>✓ Understanding Simulation Data Formats and Structures</p>
    	<p>✓ Interpreting and Analyzing Simulation Results</p>
    	""", unsafe_allow_html=True)
    
    with col2:
    	st.image("utils/images/tutorial.gif", use_column_width=True)
    
    st.write('''***''')

    st.title("Authors")
	
    # Load images and convert to base64
    satyam_image = load_image("utils/images/sat.jpg")
    arshad_image = load_image("utils/images/arshad.jpg")
	
    # Define custom CSS for circular images
    circular_image_css = """
    <style>
        .circular-image {
            border-radius: 50%;
            width: 150px;
            height: 150px;
            object-fit: cover;
        }
    </style>
    """
    
    st.markdown(circular_image_css, unsafe_allow_html=True)
    
    col1, col2 = st.columns([2, 2])
    
    with col1:
        st.markdown(
            f"""
            <div style="text-align: center;">
                <img src="data:image/jpeg;base64,{satyam_image}" class="circular-image">
                <h3>Satyam Sangeet</h3>
                <p>University of Sydney (Ph.D)</p>
                <p>CompObelisk (Founder)</p>
                <p><a href="https://scholar.google.com.au/citations?user=GgF3yTYAAAAJ&hl=en" target="_blank" style="text-decoration: none;">Google Scholar</a> |
                <a href="https://www.researchgate.net/profile/Satyam-Sangeet" target="_blank" style="text-decoration: none;">ResearchGate</a> |
                <a href="https://github.com/psychedelic2007" target="_blank" style="text-decoration: none;">Github</a> |
                <a herf="https://psychedelic2007.github.io" target="_blank" style="text-decoration: none;">Personel Webpage</a></p>
            </div>
            """,
            unsafe_allow_html=True,
        )

    with col2:
        st.markdown(
            f"""
            <div style="text-align: center;">
                <img src="data:image/jpeg;base64,{arshad_image}" class="circular-image">
                <h3>Arshad Khan</h3>
                <p>National Tsing Hua University (PhD)</p>
                <p> CompObelisk (Co-Founder)</p>
                <p><a href="https://scholar.google.com.au/citations?user=mQl293UAAAAJ&hl=en" target="_blank" style="text-decoration: none;">Google Scholar</a> |
                <a href="https://www.researchgate.net/profile/Arshad-Khan-74" target="_blank" style="text-decoration: none;">ResearchGate</a></p>
            </div>
            """,
            unsafe_allow_html=True,
        )
    
    st.write('''***''')	
    footer = """
    <div style="position: relative; bottom: 0; width: 100%; padding: 10px 0; text-align: center; margin: auto;">
        <h3>| SIMANA © 2024 | Contact </h3>
    </div>
    """
    st.markdown(footer, unsafe_allow_html=True)

if __name__ == "__main__":
    home_page()
