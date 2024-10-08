import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import streamlit as st
st.set_page_config(layout="wide", page_title="SIMANA", page_icon="utils/images/logo.svg")  # Set the page layout to wide
from streamlit.components.v1 import html
from multiapp import MultiApp
from utils.pages.home_page import home_page
from utils.pages.rmsd import rmsd
from utils.pages.rmsf import rmsf
from utils.pages.rg import rg
from utils.pages.sasa import sasa
from utils.pages.hbond import hbond
from utils.pages.dccm import dccm
from utils.pages.pca import pca
from utils.pages.bfactor import bfactor
from utils.pages.boiled_egg import boiled_egg
from utils.pages.lip import lip
from utils.pages.rama import rama
from utils.pages.contact_map import contact_map
from utils.pages.tutorial import tutorial

app = MultiApp()
app.add_app("Home Page", home_page)
app.add_app("Root Mean Square Deviation", rmsd)
app.add_app("Root Mean Square Fluctuation", rmsf)
app.add_app("Radius of Gyration", rg)
app.add_app("Solvent Accessible Surface Area", sasa)
app.add_app("Hydrogen Bond", hbond)
app.add_app("Dynamic Cross Correlation Analysis", dccm)
app.add_app("Prinicpal Component Analysis", pca)
app.add_app("Ramachandran Map", rama)
app.add_app("Contact Map", contact_map)
app.add_app("BFactor Analysis", bfactor)
app.add_app("Boiled Egg Analysis", boiled_egg)
app.add_app("Lipinski Calculation", lip)
app.add_app("Tutorial", tutorial)
app.run()
