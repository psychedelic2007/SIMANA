import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import streamlit as st
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

app = MultiApp()
app.add_app("Home Page", home_page)
app.add_app("Root Mean Square Deviation", rmsd)
app.add_app("Root Mean Square Fluctuation", rmsf)
app.add_app("Radius of Gyration", rg)
app.add_app("Solvent Accessible Surface Area", sasa)
app.add_app("Hydrogen Bond", hbond)
app.add_app("Dynamic Cross Correlation Analysis", dccm)
app.add_app("Prinicpal Component Analysis", pca)
app.add_app("BFactor Analysis", bfactor)
app.add_app("Boiled Egg Analysis", boiled_egg)
app.run()
