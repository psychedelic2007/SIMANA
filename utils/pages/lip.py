import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
from io import StringIO, BytesIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_lipinski(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Calculate Lipinski's parameters
    mol_weight = Descriptors.MolWt(mol)
    num_h_donors = Descriptors.NumHDonors(mol)
    num_h_acceptors = Descriptors.NumHAcceptors(mol)
    logp = Descriptors.MolLogP(mol)

    # Check if the molecule follows Lipinski's Rule of Five
    violations = []
    if mol_weight > 500:
        violations.append("MolWt > 500")
    if num_h_donors > 5:
        violations.append("HDonors > 5")
    if num_h_acceptors > 10:
        violations.append("HAcceptors > 10")
    if logp > 5:
        violations.append("LogP > 5")

    follows_rule = "No" if violations else "Yes"
    violation_details = ", ".join(violations) if violations else "--"

    return mol_weight, num_h_donors, num_h_acceptors, logp, follows_rule, violation_details

def plot_distributions(data):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))

    sns.histplot(data['MolWt'], ax=axs[0, 0], kde=True, linewidth=2.0, edgecolor="black")
    kde = axs[0, 0].lines[0]
    kde.set_linewidth(3.0)
    axs[0, 0].set_title('Molecular Weight')

    sns.histplot(data['HDonors'], ax=axs[0, 1], kde=True, linewidth=2.0, edgecolor="black")
    kde = axs[0, 1].lines[0]
    kde.set_linewidth(3.0)
    axs[0, 1].set_title('Number of Hydrogen Donors')

    sns.histplot(data['HAcceptors'], ax=axs[1, 0], kde=True, linewidth=2.0, edgecolor="black")
    kde = axs[1, 0].lines[0]
    kde.set_linewidth(3.0)
    axs[1, 0].set_title('Number of Hydrogen Acceptors')

    sns.histplot(data['LogP'], ax=axs[1, 1], kde=True, linewidth=2.0, edgecolor="black")
    kde = axs[1, 1].lines[0]
    kde.set_linewidth(3.0)
    axs[1, 1].set_title('LogP')

    plt.tight_layout()

    return fig

def lip():
    st.title("Lipinski's Rule of Five Calculator")

    # File upload
    uploaded_file = st.file_uploader("Upload a text file containing SMILES notations", type=["txt"])

    if uploaded_file is not None:
        smiles_list = uploaded_file.read().decode("utf-8").splitlines()

        if st.button("Submit"):
            data = {
                "SMILES": [],
                "MolWt": [],
                "HDonors": [],
                "HAcceptors": [],
                "LogP": [],
                "Follows Lipinski": [],
                "Violations": []
            }

            for smiles in smiles_list:
                if smiles.strip():  # Ensure it's not an empty line
                    result = calculate_lipinski(smiles)
                    if result:
                        mol_weight, num_h_donors, num_h_acceptors, logp, follows_rule, violation_details = result
                        data["SMILES"].append(smiles)
                        data["MolWt"].append(mol_weight)
                        data["HDonors"].append(num_h_donors)
                        data["HAcceptors"].append(num_h_acceptors)
                        data["LogP"].append(logp)
                        data["Follows Lipinski"].append(follows_rule)
                        data["Violations"].append(violation_details)

            df = pd.DataFrame(data)

            # Store the DataFrame and plot in session state
            st.session_state['processed_data'] = df
            st.session_state['plot_figure'] = plot_distributions(df)

            # Display the plot
            st.pyplot(st.session_state['plot_figure'])

            # Provide download buttons
            st.download_button(
                label="Download Lipinski's Parameters CSV",
                data=st.session_state['processed_data'].to_csv(index=False),
                file_name="lipinski_parameters.csv",
                mime="text/csv"
            )

            buf = BytesIO()
            st.session_state['plot_figure'].savefig(buf, format="png", dpi=300)  # Save with 300 DPI
            buf.seek(0)

            st.download_button(
                label="Download Distribution Plot",
                data=buf,
                file_name="lipinski_distribution.png",
                mime="image/png"
            )

if __name__ == "__main__":
    lip()
