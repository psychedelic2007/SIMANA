import streamlit as st

# Function to display a tutorial section
def display_tutorial(title, description, formula, formula_description, input_requirements, output_description, interpretation):
    with st.expander(title):
        st.markdown(description)
        st.latex(formula)
        st.markdown(formula_description)
        st.markdown("**Input Requirements:**")
        st.markdown(input_requirements)
        st.markdown("**Output Files:**")
        st.markdown(output_description)
        st.markdown("**Interpretation:**")
        st.markdown(interpretation)

def tutorial():
    # Title of the page
    col1, col2 = st.columns([1, 3])
    with col1:
    	st.image("utils/images/video-lesson.png", use_column_width=True)
    with col2:
        st.markdown("""
        <h1>Welcome to <span style="color:cyan;">SIM</span>ulation <span style="color:cyan;">ANA</span>lysis Tutorial Page</h1>
        <p>Here you will find detailed explanations and instructions on various topics.</p>
        """, unsafe_allow_html=True)

    # RMSD Section
    rmsd_description = """
    RMSD (Root Mean Square Deviation) is a measure of the average distance between atoms 
    (usually the backbone atoms) of superimposed proteins. It is commonly used to quantify 
    the structural similarity between two proteins or the deviation of a structure over time 
    in a molecular dynamics simulation.
    """
    rmsd_formula = r"""
    \text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^{N} \left( \mathbf{r}_i(t) - \mathbf{r}_i^{ref} \right)^2}
    """
    rmsd_formula_description = r"""
    **Where**,
    - $\mathbf{r}_i$ is the position of the 𝑖th atom in the current conformation.
    - $\mathbf{r}_i^{ref}$ is the position of the 𝑖th atom in the reference conformation.
    - $N$ is the total number of atoms.
    """
    rmsd_input_requirements = """
    - A data file (.xvg) is required.
    - The input files should be properly formatted and preprocessed as needed.
    """
    rmsd_output_description = """
    - **RMSD Plot**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - **RMSD Distribution Plot**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - All plots are provided at 300 DPI.
    """
    rmsd_interpretation = """
    - **RMSD Plot**: This plot displays the RMSD value of the protein or the system versus the simulation time. 
      If the fluctuations increase, it suggests that the system is becoming less stable; a decrease indicates greater stability.
    - **RMSD Distribution Plot**: This plot shows the RMSD values on the x-axis and their frequency on the y-axis, 
      providing insight into where the majority of RMSD values are clustered.
    """

    display_tutorial("Root Mean Square Deviation (RMSD)", rmsd_description, rmsd_formula, rmsd_formula_description, rmsd_input_requirements, rmsd_output_description, rmsd_interpretation)

    # RMSF Section
    rmsf_description = """
    RMSF (Root Mean Square Fluctuation) measures the flexibility of individual residues 
    or atoms over time in a molecular dynamics simulation. It provides insight into the 
    dynamic behavior of specific parts of a protein.
    """
    rmsf_formula = r"""
    \text{RMSF} = \sqrt{\frac{1}{T} \sum_{t=1}^{T} \left( \mathbf{r}_i(t) - \langle \mathbf{r}_i \rangle \right)^2}
    """
    rmsf_formula_description = r"""
    **Where**,
    - $\mathbf{r}_i(t)$ is the position of atom $i$ at time $t$.
    - $\langle \mathbf{r}_i \rangle$ is the average position of atom $i$ over the trajectory.
    - $T$ is the total number of time frames in the trajectory
    """
    rmsf_input_requirements = """
    - A data file (.xvg) is required.
    - The input files should be properly formatted and preprocessed as needed.
    """
    rmsf_output_description = """
    - **RMSF Plot**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - **RMSD Distribution**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - All plots are provided at 300 DPI.
    """
    rmsf_interpretation = """
    - **RMSF Plot**: This plot illustrates the flexibility of each residue in the protein over time. 
      Higher RMSF values indicate greater flexibility, while lower values suggest rigidity.
    """

    display_tutorial("Root Mean Square Fluctuation (RMSF)", rmsf_description, rmsf_formula, rmsf_formula_description, rmsf_input_requirements, rmsf_output_description, rmsf_interpretation)

    # Radius of Gyration (RG) Section
    rg_description = """
    Radius of Gyration (Rg) is a measure of the compactness of a molecular structure. 
    It is the root mean square distance of the collection of atoms from their common center of mass.
    """
    rg_formula = r"""
    R_g = \sqrt{\frac{\sum_{i=1}^{N} m_i \left( \mathbf{r}_i - \mathbf{R}_{cm} \right)^2}{\sum_{i=1}^{N} m_i}}
    """
    rg_formula_description = r"""
    **Where**,
    - $\mathbf{r}_i$ is the position of $ith$ atom.
    - $\mathbf{R}_{cm}$ is the center of mass of the system
    - $N$ is the total number of atoms
    """
    rg_input_requirements = """
    - A data file (.xvg) is required.
    - The input files should be properly formatted and preprocessed as needed.
    """
    rg_output_description = """
    - **Radius of Gyration Plot**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - **Radius of Gyration Distribution**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - All plots are provided at 300 DPI.
    """
    rg_interpretation = """
    - **Radius of Gyration Plot**: This plot indicates the compactness of the protein structure over time. 
      A lower Rg value suggests a more compact structure, while a higher value indicates a more extended conformation.
    """

    display_tutorial("Radius of Gyration (RG)", rg_description, rg_formula, rg_formula_description, rg_input_requirements, rg_output_description, rg_interpretation)

    # Solvent Accessible Surface Area (SASA) Section
    sasa_description = """
    SASA (Solvent Accessible Surface Area) is the surface area of a biomolecule that is accessible to a solvent. 
    It is an important measure in understanding protein folding and stability.
    """
    sasa_formula = r"""
    \text{SASA} = \sum_{i=1}^{N} A_i
    """
    sasa_formula_description = r"""
    **Where**,
    - $A_i$ is the Solvent accessible area of the $ith$ atom.
    """
    sasa_input_requirements = """
    - A data file (.xvg) is required.
    - The input files should be properly formatted and preprocessed as needed.
    """
    sasa_output_description = """
    - **SASA Plot**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - **SASA Distribution**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - All plots are provided at 300 DPI.
    """
    sasa_interpretation = """
    - **SASA Plot**: This plot shows the solvent-accessible surface area of the protein over time. 
      Changes in SASA can indicate folding or unfolding events within the protein structure.
    """

    display_tutorial("Solvent Accessible Surface Area (SASA)", sasa_description, sasa_formula, sasa_formula_description, sasa_input_requirements, sasa_output_description, sasa_interpretation)

    # Hydrogen Bond Section
    hbond_description = """
    Hydrogen bonds are weak interactions that play a crucial role in stabilizing the structure of biomolecules. 
    This section helps in analyzing the number and pattern of hydrogen bonds during a molecular dynamics simulation.
    """
    hbond_formula = r"""
    \text{Criteria for H-Bond: } \text{Distance} < 3.5\ \text{Å and Angle} > 120^\circ
    """
    hbond_formula_description = r"""
    **Where**,
    - **$Distance < 3.5 Å$** criterion specifies that for a hydrogen bond to be considered, the distance between the hydrogen atom (donor) and the acceptor atom (typically a lone pair of electrons on another atom) must be less than 3.5 angstroms (Å). This is the straight-line distance between the hydrogen atom and the acceptor atom.
    - **$Angle > 120^\circ$** criterion states that the angle formed between the donor atom, the hydrogen atom, and the acceptor atom must be greater than 120 degrees. Specifically, this is the angle between the line connecting the donor to the hydrogen and the line connecting the hydrogen to the acceptor. This angle is often referred to as the donor-hydrogen-acceptor (D-H-A) angle. The angle ensures that the hydrogen bond is in a favorable geometric arrangement, which is typically a linear or nearly linear configuration.
    """
    hbond_input_requirements = """
    - A data file (.xvg) is required.
    - The input files should be properly formatted and preprocessed as needed.
    """
    hbond_output_description = """
    - **Hydrogen Bond Plot**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - **Hydrogen Bond Distribution**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - All plots are provided at 300 DPI.
    """
    hbond_interpretation = """
    - **Hydrogen Bond Plot**: This plot displays the number of hydrogen bonds over time. 
      A stable number of hydrogen bonds suggests a stable structure, while fluctuations might indicate changes in stability.
    """

    display_tutorial("Hydrogen Bond", hbond_description, hbond_formula, hbond_formula_description, hbond_input_requirements, hbond_output_description, hbond_interpretation)
    
    #DCCM Section
    dccm_description = """
    The Dynamic Cross-Correlation Matrix (DCCM) is used to analyze the correlated motions between different atoms or residues in a protein during a molecular dynamics simulation. 
    It helps to understand how the movements of one part of the protein are related to those of another part.
    """
    dccm_formula = r"""
    \text{DCCM}_{ij} = \frac{\langle (x_i(t) - \langle x_i \rangle)(x_j(t) - \langle x_j \rangle) \rangle}{\sqrt{\langle (x_i(t) - \langle x_i \rangle)^2 \rangle \langle (x_j(t) - \langle x_j \rangle)^2\rangle}}
    """
    dccm_formula_description = r"""
    **Where**,
    - $x_i(t)$ represents the position (or other coordinate) of atoms $i$ at time $t$ in the simulation.
    - $\langle x_i \rangle$ is the average position (or coordinate) of atom $i$ over the entire simulation time. It’s calculated as:
                                                                       $\lt x_{i}\gt = \frac{1}{T}\sum_{t=1}^{T}x_{i}(t)$
    where $T$ is the total number of time frames in teh simulation
    - $\langle (x_i(t) - \langle x_i \rangle)(x_j(t) - \langle x_j \rangle) \rangle$ term represents the covariance between the positions of atoms $i$ and $j$ over time. It measures how changes in the position of atom $i$ are related to changes in the position of atom $j$. Covariance indicates whether the movements of the two atoms are correlated. A positive value suggests that when one atom moves in one direction, the other tends to move in the same direction. A negative value suggests they move in opposite directions.
    - $\langle (x_i(t) - \langle x_i \rangle)^2 \rangle$ term is the variance of the position of atom $i$ over time. It measures the extent to which the position of atom $i$ fluctuates around its average position.
    - $\sqrt{{\langle (x_{i}(t)-\langle x_{i}\rangle )^{2}\rangle}{\langle (x_{j}(t)-\langle x_{j}\rangle )^{2}\rangle}} $. This denominator normalizes the covariance, so the DCCM values range between -1 and 1. It ensures that the DCCM reflects the strength of correlation independent of the magnitude of fluctuations.
    """
    dccm_input_requirements = """
    - A trajectory file (.xtc) and a reference structure (.pdb) are required.
    - The input files should be properly formatted and preprocessed as needed.
    """
    dccm_output_description = """
    - **DCCM Matrix**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - All plots are provided at 300 DPI.
    """
    dccm_interpretation = """
    - **DCCM Matrix**: The matrix shows the correlation of movements between different atoms or residues. Positive values indicate correlated motions, while negative values indicate anti-correlated motions. The intensity of the color can help identify regions with significant correlated or anti-correlated motions.
    """
    display_tutorial("Dynamic Cross Correlation Matrix (DCCM)", dccm_description, dccm_formula, dccm_formula_description, dccm_input_requirements, dccm_output_description, dccm_interpretation)

    #PCA Section
    pca_description = """
    Principal Component Analysis (PCA) is a dimensionality reduction technique that transforms the data into a set of orthogonal components (principal components) 
    that capture the most variance in the data. PCA is used to simplify the analysis of complex data sets by reducing the number of variables while preserving as much 
    information as possible.
    """
    pca_formula = r"""
    \text{PC}_k = \sum_{i=1}^{N} w_{ik} x_i
    """
    pca_formula_description = r"""
    **Where**,
    - $PC_{k}$ represents the $kth$ principal component. Principal components are new variables constructed as linear combinations of the original variables, and they capture the directions of maximum variance in the data.
    - $w_{ik}$ are the weights (or loadings) for the $ith$ variable on the $kth$ principal component. They are elements of the eigenvector associated with the $kth$ principal component. Each weight determines the contribution of the $ith$ original variable to the $kth$ principal component.
    - $x_{i}$ represents the $ith$ original variable or feature in the dataset. It is the value of the $ith$ variable for a given observation.
    - $sum_{i=1}^{N}$ this summation adds up the contributions of all original variables $x_{i}$, weighted by their corresponding weights $w_{ik}$, to compute the value of the $kth$ principal component.
    """
    pca_input_requirements = """
    - A data file (.xvg) is required. Typically the .xvg file that you obtaine after performing **_gmx anaeig_**.
    - The data should be preprocessed and normalized if necessary.
    """
    pca_output_description = """
    - **PCA Plot**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - All plots are provided at 300 DPI.
    """
    pca_interpretation = """
    - **PCA Plot**: The plot displays the data in terms of principal components. The first few principal components usually capture the most significant variance in the data. 
    Examining the scatter plots can reveal clusters, trends, or outliers.
    """
    display_tutorial("Principal Component Analysis (PCA)", pca_description, pca_formula, pca_formula_description, pca_input_requirements, pca_output_description, pca_interpretation)

    #Ramachandran section
    ramachandran_description = """
    The Ramachandran Plot is a graphical representation of the φ (phi) and ψ (psi) dihedral angles of the amino acid residues in a protein structure. 
    It is used to evaluate the conformational angles and validate the geometry of protein structures.
    """
    ramachandran_formula = r"""
    \text{Ramachandran Plot} = \text{Plot of } \phi \text{ vs. } \psi \text{ dihedral angles}
    """
    ramachandran_formula_description = r"""
    **Where**,
    - **$\phi$ Angle** is the dihedral angle around the bond between the nitrogen (N) and alpha carbon ($C_{alpha}$) atoms of a protein backbone. It describes the rotation about this bond.
    - **$\psi$ Angle** is the dihedral angle around the bond between the alpha carbon ($C_{alpha}$) and carbonyl carbon (C=O) atoms. It describes the rotation about this bond.
    - In a Ramachandran plot, the x-axis represents the $\phi$ dihedral angle, while the y-axis represents the $\psi$ dihedral angle. Each point on the plot corresponds to a specific $\phi$ and $\psi$ angle pair for a given residue in the protein structure.
    - The plot shows various regions where certain combinations of $\phi$ and $\psi$ angles are allowed based on steric hindrance and favorable interactions. These regions typically correspond to common protein secondary structures, such as alpha-helices and beta-sheets. For example, alpha-helices tend to cluster in one region of the plot, while beta-sheets cluster in another. Unfavorable or disallowed regions are those where steric clashes would occur, making those angles less favorable for stable protein structures.
    """
    ramachandran_input_requirements = """
    - A PDB file is required.
    - The file should be properly formatted and preprocessed as needed.
    """
    ramachandran_output_description = """
    - **Ramachandran Plot**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - All plots are provided at 300 DPI.
    """
    ramachandran_interpretation = """
    - **Ramachandran Plot**: The plot shows the distribution of φ and ψ angles. Regions in the plot correspond to allowed and disallowed conformations. 
    The plot helps in assessing the quality of the protein structure; a higher number of residues in the favored regions indicates a good-quality structure.
    """
    display_tutorial("Ramachandran Analysis", ramachandran_description, ramachandran_formula, ramachandran_formula_description, ramachandran_input_requirements, ramachandran_output_description, ramachandran_interpretation)

    # Contact Map Section
    contact_map_description = """
    A contact map is a representation of the distances between pairs of residues in a protein structure. 
    It is used to understand the spatial proximity and interactions within the protein.
    """
    contact_map_formula = r"""
    \text{Contact} = \begin{cases} 
    1 & \text{if } d_{ij} \leq \text{cutoff} \\ 
    0 & \text{otherwise} 
    \end{cases}
    """
    contact_map_formula_description = r"""
    **Where**,
    - **Contact** is a binary value indicating whether two atoms (or residues) are considered to be in contact or not.
    - $d_{ij}$ is the distance between atom $i$ and atom $j$. It is typically measured in units such as angstroms (Å).
    - **Cutoff** is a predefined distance threshold (cutoff value) used to determine whether the atoms are in contact. It is a critical value that specifies the maximum distance within which two atoms are considered to be interacting or in contact.
    - If $d_{ij} \leq cutoff$ i.e. the distance between atoms $i$ and $j$ is less than or equal to the cutoff distance, then the Contact value is set to 1. This means that the atoms are considered to be in contact. **Otherwise**, If the distance between atoms $i$ and $j$ is greater than the cutoff distance, then the Contact value is set to 0. This means that the atoms are not considered to be in contact.
    - This classification helps in analyzing interactions, such as hydrogen bonds, van der Waals interactions, or other contact-based phenomena.
    - **Interaction Analysis**: By applying this formula to all pairs of atoms or residues in a structure, researchers can generate a contact map or matrix, which helps in understanding the spatial organization and interactions within the molecule
    """
    contact_map_input_requirements = """
    - A PDB file is required.
    - The input files should be properly formatted and preprocessed as needed.
    """
    contact_map_output_description = """
    - **Contact Map Plot**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - All plots are provided at 300 DPI.
    """
    contact_map_interpretation = """
    - **Contact Map Plot**: This plot visualizes the contacts between residues within the protein structure. 
      Persistent contacts indicate stable interactions, while fluctuating contacts may suggest dynamic regions.
    """

    display_tutorial("Contact Map", contact_map_description, contact_map_formula, contact_map_formula_description, contact_map_input_requirements, contact_map_output_description, contact_map_interpretation)

    # B-Factor Section
    bfactor_description = """
    B-Factor (also known as Debye-Waller factor) is a measure of the atomic displacement or flexibility 
    in a crystal structure. It is often used to assess the quality of protein structures.
    """
    bfactor_formula = r"""
    B_i = \frac{8 \pi^2}{3} \langle u_i^2 \rangle
    """
    bfactor_formula_description = r"""
    **Where**,
    - $B_{i}$ represents the B-factor (or temperature factor) of atom $i$. It quantifies the extent of thermal vibration or disorder of the atom within the crystal structure. Higher B-factors indicate greater atomic displacement or flexibility.
    - $\frac{8\pi^{2}}{3}$ is a constant factor that converts the mean square displacement $(\langle u^{2}_{i} \rangle)$ into the B-factor. The constant derives from the relation between the atomic mean square displacement and the B-factor in crystallography.
    - $(\langle u^{2}_{i} \rangle)$ denotes the mean square displacement of atom $i$. It represents the average of the squared displacements of the atom from its equilibrium position due to thermal vibrations or other sources of disorder.
    - **Thermal Vibrations**: The B-factor provides insight into how much an atom is vibrating around its average position. It is important in understanding the dynamic aspects of the molecular structure. Atoms with high B-factors are typically more disordered or flexible, whereas atoms with low B-factors are more rigid and well-ordered.
    - **Crystallography**: In X-ray crystallography, the B-factor is used to account for the thermal motion of atoms and the accuracy of the electron density map. It helps refine the atomic positions and improve the quality of the crystallographic model.
    """
    bfactor_input_requirements = """
    - A PDB file is required.
    - The input files should be properly formatted and preprocessed as needed.
    """
    bfactor_output_description = """
    - **B-Factor Plot**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - All plots are provided at 300 DPI.
    """
    bfactor_interpretation = """
    - **B-Factor Plot**: This plot indicates the atomic displacement within the protein structure. 
      Higher B-factor values correspond to greater atomic flexibility.
    """

    display_tutorial("B-Factor", bfactor_description, bfactor_formula, bfactor_formula_description, bfactor_input_requirements, bfactor_output_description, bfactor_interpretation)
    
    #Boiled egg section
    boiled_egg_description = """
    The Boiolod Egg Plot is used to visualize the flexibility of residues in a protein by plotting their average B-factors. It provides insight into the regions of the protein 
    that exhibit high or low flexibility.
    """
    boiled_egg_formula = r"""
    \begin{align}
    \text{WLogP} &= \log_{10}\left( \frac{C_{\text{oct}}}{C_{\text{water}}} \right) \\
    \text{TPSA} &= \sum_{i=1}^{n} \text{SA}_{i}
    \end{align}
    """
    boiled_egg_formula_description = r"""
    **where**,
    - $WLogP = log_{10}\left( \frac{C_{oct}}{C_{water}} \right)$ corresponds to the Lipophilicity of the molecule.
    	- **LogP (Partition Coefficient)**: measures the distribution of a compound between an octanol phase and a water phase. It is an indicator of the compound's lipophilicity.
    	- **$C_{oct}$**: Concentration of the compound in octanol.
    	- **$C_{water}$**: Concentration of the compound in water.
    	- **Purpose**: LogP is used to assess the hydrophobic (lipophilic) nature of molecules, which is critical for understanding drug absorption, distribution, and permeability.
    - $TPSA = \sum_{i=1}^{N}SA_{i}$ represents the polarity of the molecule.
    	- **TPSA** is the sum of the polar surface areas of all polar atoms (usually oxygen, nitrogen, and their hydrogens) in a molecule.
    	- **Surface Area**: This is typically calculated by considering the contributions of each polar atom and its associated hydrogen atoms.
    	- **Purpose**: TPSA is a measure of the molecular polarity and hydrogen bonding potential, which affects the molecule’s solubility, permeability, and overall drug-likeness.
    - **Lipophilicity (LogP)**: A measure of how well a molecule partitions between an organic solvent (usually octanol) and water. Higher LogP values indicate higher lipophilicity.
    - **Polarity (TPSA)**: The polar surface area of a molecule, which provides an estimate of how polar a molecule is.
    - **Plot Characteristics**: Molecules are plotted in a two-dimensional space where the x-axis represents LogP and the y-axis represents TPSA. This plot helps in visualizing how the molecule's properties compare with ideal values for drug-like behavior.
    """
    boiled_egg_input_requirements = """
    - SMILE notation of the ligand molecule can be eneterd in the text box. **_Make sure you enter one SMILES per line_**.
    """
    boiled_egg_output_description = """
    - **Boiolod Egg Plot**: Can be downloaded in the following formats: JPG, PNG, SVG, EPS, PDF.
    - All plots are provided at 300 DPI.
    """
    boiled_egg_interpretation = """
    - **Boiolod Egg Plot**: The plot shows the average B-factors of residues. Higher B-factors indicate higher flexibility, while lower B-factors suggest rigidity. 
    This helps in identifying flexible and rigid regions within the protein structure.
    """
    display_tutorial("Boiled Egg Analysis", boiled_egg_description, boiled_egg_formula, boiled_egg_formula_description, boiled_egg_input_requirements, boiled_egg_output_description, boiled_egg_interpretation)
    
    #Lipinski section
    lipinski_description = """
    Lipinski's Rule of Five is used to evaluate the drug-likeness of a compound. It assesses factors like molecular weight, lipophilicity, and the number of hydrogen bond donors 
    and acceptors to predict the likelihood of oral bioavailability.
    """
    lipinski_formula = r"""
    \text{Rule of Five:} \\
    \text{1. Molecular weight } \leq 500 \\
    \text{2. LogP } \leq 5 \\
    \text{3. Hydrogen bond donors } \leq 5 \\
    \text{4. Hydrogen bond acceptors } \leq 10
    """
    lipinski_formula_description = r"""
    - **Molecular Weight $\le 500$**
    	- **Explanation**: This criterion states that for a drug to have good oral bioavailability, its molecular weight should ideally be less than or equal to 500 Daltons. Larger molecules may have difficulties in crossing cellular membranes and may have poor absorption in the gastrointestinal tract.
    - **LogP $\le 5$**
    	- **Explanation**: LogP (or partition coefficient) measures the lipophilicity of the compound, indicating how well it partitions between octanol and water. A LogP value of 5 or less suggests that the compound is balanced in terms of hydrophobicity and hydrophilicity, which is generally favorable for absorption and distribution. Compounds with higher LogP values may be too lipophilic, which could lead to poor solubility in aqueous environments.
    - **Hydrogen Bond Donors $\le 5$**
    	- **Explanation**: This criterion limits the number of hydrogen bond donors (e.g., -OH and -NH groups) in a molecule to 5 or fewer. Excessive hydrogen bonding can increase the molecule’s polarity, making it less likely to permeate cell membranes and potentially affecting oral bioavailability.
    - **Hydrogen Bond Acceptors $\le 10$**
    	- **Explanation**: This criterion limits the number of hydrogen bond acceptors (e.g., oxygen and nitrogen atoms) in a molecule to 10 or fewer. Like hydrogen bond donors, excessive hydrogen bond acceptors can also affect a compound’s ability to cross biological membranes and impact its absorption and permeability.
    """
    lipinski_input_requirements = """
    - The input should be a chemical structure file in SMILES format.
    - If you want to calculate Lipinski for a single molecule then paste the corresponding SMILES formula in the provided box.
    - If you want to perform Lipinski's calculation for bulk amount of ligands then upload a **_.txt_** file containing SMILES notation of each ligand. Make sure that the text file has only one SMILES per line.
    """
    lipinski_output_description = """
    - **Lipinski Calculation Report**: Can be downloaded in the following formats: CSV.
    - The CSV file will contain seven columns: [SMILES | Molecular Weight | Hydrogen Donors | Hydrogen Acceptors | LogP | Follows Lipinski or not | Violations]
    """
    lipinski_interpretation = """
    - **Lipinski Calculation Report**: The report evaluates whether the compound adheres to Lipinski’s Rule of Five. 
    Compounds that meet these criteria are more likely to be orally bioavailable. Deviations from the rule suggest further testing or modification.
    """
    display_tutorial("Lipinski Analysis", lipinski_description, lipinski_formula, lipinski_formula_description, lipinski_input_requirements, lipinski_output_description, lipinski_interpretation)

if __name__ == "__main__":
    tutorial()
