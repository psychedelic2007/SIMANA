import streamlit as st
from utils.operations.rmsd import rmsd

def analysis():
    if 'operation' not in st.session_state:
        st.session_state.operation = None

    if st.session_state.operation is None:
        st.title("Analysis Page")
        st.write("Choose an operation to perform:")
        
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("RMSD"):
                st.session_state.operation = "rmsd"
                st.experimental_rerun()
            if st.button("Multiplication"):
                st.session_state.operation = "multiplication"
                st.experimental_rerun()
        
        with col2:
            if st.button("Subtraction"):
                st.session_state.operation = "subtraction"
                st.experimental_rerun()
            if st.button("Division"):
                st.session_state.operation = "division"
                st.experimental_rerun()
    else:
        if st.session_state.operation == "rmsd":
            rmsd()
        #elif st.session_state.operation == "subtraction":
        #    subtraction()
        #elif st.session_state.operation == "multiplication":
        #    multiplication()
        #elif st.session_state.operation == "division":
        #    division()
        
        if st.button("Back to Operations"):
            st.session_state.operation = None
            st.experimental_rerun()
