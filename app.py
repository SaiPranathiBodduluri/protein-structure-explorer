
import streamlit as st
import requests
import streamlit.components.v1 as components
import json
import base64
import pandas as pd
from io import StringIO

st.set_page_config(
    page_title="Protein Structure Explorer",
    page_icon="🧬",
    layout="wide"
)

st.markdown("""
<style>
.block-container {
    padding-top: 2rem;
    max-width: 1200px;
}
.hero {
    background: linear-gradient(120deg, #0f766e, #164e63);
    padding: 30px;
    border-radius: 16px;
    color: white;
    margin-bottom: 25px;
}
.hero h1 {
    color: white;
    margin-bottom: 8px;
}
.hero p {
    color: white;
    font-size: 16px;
}
div[data-testid="stMetric"] {
    background: #f0fdfa;
    padding: 15px;
    border-radius: 12px;
    border: 1px solid #ccfbf1;
}
</style>
""", unsafe_allow_html=True)

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/"

st.markdown("""
<div class="hero">
    <h1>🧬 Protein Structure Explorer</h1>
    <p>Search proteins, explore amino acid sequences, and visualize
    predicted 3D structures in one place.</p>
</div>
""", unsafe_allow_html=True)

st.subheader("🔎 Search Protein")

search_type = st.selectbox(
    "Choose search method",
    [
        "UniProt Protein ID",
        "Protein Name",
        "Organism",
        "Amino Acid Sequence"
    ]
)

if search_type == "UniProt Protein ID":
    query = st.text_input("Enter UniProt ID", "P69905")
elif search_type == "Protein Name":
    query = st.text_input("Enter protein name", "Insulin")
elif search_type == "Organism":
    query = st.text_input("Enter organism", "Homo sapiens")
else:
    query = st.text_area(
        "Enter amino acid sequence",
        placeholder="Paste a protein sequence here"
    )


def search_proteins(search_type, query):
    if search_type == "UniProt Protein ID":
        response = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{query.strip().upper()}.json",
            timeout=30
        )
        response.raise_for_status()
        item = response.json()

        names = item.get("proteinDescription", {})
        recommended = names.get("recommendedName", {})
        protein_name = recommended.get("fullName", {}).get("value")

        if not protein_name:
            subtypes = names.get("submissionNames", [])
            if subtypes:
                protein_name = subtypes[0].get("fullName", {}).get("value")

        sequence = item.get("sequence", {}).get("value", "")

        return [{
            "Entry": item.get("primaryAccession", query.strip().upper()),
            "Protein names": protein_name or "Unknown",
            "Organism": item.get("organism", {}).get(
                "scientificName", "Unknown"
            ),
            "Length": len(sequence),
            "Sequence": sequence
        }]

    if search_type == "Protein Name":
        search_query = f'"{query.strip()}"'

    elif search_type == "Organism":
        search_query = f'organism_name:"{query.strip()}"'

    else:
        sequence = "".join(query.split()).upper()
        valid_letters = set("ACDEFGHIKLMNPQRSTVWY")

        if not sequence or not set(sequence) <= valid_letters:
            st.error("Enter a valid amino acid sequence.")
            return []

        search_query = f"sequence:{sequence}"

    response = requests.get(
        UNIPROT_API,
        params={
            "query": search_query,
            "format": "tsv",
            "size": 10,
            "fields": "accession,protein_name,organism_name,length,sequence"
        },
        timeout=40
    )
    response.raise_for_status()

    df = pd.read_csv(StringIO(response.text), sep="\t")
    return df.to_dict("records")


if st.button("Search Protein", type="primary", use_container_width=True):
    if not query.strip():
        st.warning("Please enter a search value.")
    else:
        try:
            with st.spinner("Searching protein database..."):
                results = search_proteins(search_type, query)

            st.session_state["protein_results"] = results
            st.session_state["loaded_protein"] = None

            if results:
                st.success(f"Found {len(results)} protein result(s).")
            else:
                st.warning("No matching proteins found.")

        except requests.RequestException as e:
            st.error(f"Could not connect to UniProt: {e}")
        except Exception as e:
            st.error(f"Search failed: {e}")


results = st.session_state.get("protein_results", [])

if results:
    st.subheader("📋 Search Results")

    labels = [
        f"{p.get('Entry', '')} — {p.get('Protein names', 'Unknown')}"
        for p in results
    ]

    selected = st.selectbox(
        "Select a protein",
        range(len(results)),
        format_func=lambda i: labels[i]
    )

    protein = results[selected]
    protein_id = protein.get("Entry", "")

    st.dataframe(
        pd.DataFrame([{
            "UniProt ID": protein_id,
            "Protein Name": protein.get("Protein names", "Unknown"),
            "Organism": protein.get("Organism", "Unknown"),
            "Sequence Length": protein.get("Length", "")
        }]),
        hide_index=True,
        use_container_width=True
    )

    if st.button(
        "View 3D Structure",
        type="primary",
        use_container_width=True
    ):
        try:
            with st.spinner("Retrieving AlphaFold structure..."):
                response = requests.get(
                    ALPHAFOLD_API + protein_id,
                    timeout=30
                )
                response.raise_for_status()
                structures = response.json()

                if not structures:
                    st.warning(
                        "No AlphaFold structure is available for this protein."
                    )
                else:
                    structure = structures[0]
                    pdb_response = requests.get(
                        structure["pdbUrl"], timeout=30
                    )
                    pdb_response.raise_for_status()

                    sequence = (
                        structure.get("uniprotSequence", "")
                        or protein.get("Sequence", "")
                    )

                    st.session_state["loaded_protein"] = {
                        "id": protein_id,
                        "name": structure.get(
                            "uniprotDescription",
                            protein.get("Protein names", "Unknown")
                        ),
                        "organism": structure.get(
                            "organismScientificName",
                            protein.get("Organism", "Unknown")
                        ),
                        "sequence": sequence,
                        "pdb": pdb_response.text
                    }

        except Exception as e:
            st.error(f"Could not retrieve the structure: {e}")


loaded = st.session_state.get("loaded_protein")

if loaded:
    st.divider()
    st.subheader("🧬 Protein Analysis Results")
    st.markdown(f"### {loaded['name']}")

    c1, c2, c3 = st.columns(3)
    c1.metric("UniProt ID", loaded["id"])
    c2.metric("Sequence Length", f"{len(loaded['sequence'])} amino acids")
    c3.metric("Structure Status", "Available")

    st.write("**Organism:**", loaded["organism"])

    with st.expander("View Amino Acid Sequence"):
        st.code(loaded["sequence"])

    st.subheader("🔬 Interactive 3D Protein Structure")

    pdb_json = json.dumps(loaded["pdb"])

    html = """
    <div id="viewer" style="width:100%;height:480px;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
    const pdb = PDB_PLACEHOLDER;
    const viewer = $3Dmol.createViewer(
        document.getElementById("viewer"),
        {backgroundColor:"white"}
    );
    viewer.addModel(pdb, "pdb");
    viewer.setStyle({}, {cartoon:{color:"spectrum"}});
    viewer.zoomTo();
    viewer.render();
    </script>
    """.replace("PDB_PLACEHOLDER", pdb_json)

    components.html(html, height=500)

    encoded = base64.b64encode(loaded["pdb"].encode()).decode()

    st.markdown(
        f'<a download="{loaded["id"]}_structure.pdb" '
        f'href="data:chemical/x-pdb;base64,{encoded}">'
        '📥 Download Protein Structure (PDB)</a>',
        unsafe_allow_html=True
    )

st.divider()
st.caption(
    "Protein Structure Explorer | UniProt protein search | "
    "AlphaFold predicted structures"
)
