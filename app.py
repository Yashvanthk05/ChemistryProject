from flask import Flask, request, jsonify, render_template
from rdkit import Chem
from rdkit.Chem import AllChem

app = Flask(__name__)

# Mapping for compounds with SMILES and chemical formulas
compound_smiles_mapping = {
    "paracetamol": {
        "smiles": "CC(=O)NC1=CC=C(O)C=C1",
        "name": "Paracetamol",
        "formula": "C8H9NO2",
        "description": "Paracetamol, also known as acetaminophen, is an analgesic and antipyretic drug used to relieve mild to moderate pain and reduce fever.",
    },
    "methane": {
        "smiles": "C",
        "name": "Methane",
        "formula": "CH4",
        "description": "Methane is a simple hydrocarbon and the primary component of natural gas. It is used as a fuel and is also a potent greenhouse gas."
    },
    "ethane": {
        "smiles": "CC",
        "name": "Ethane",
        "formula": "C2H6",
        "description": "Ethane is a hydrocarbon and an alkane used in the petrochemical industry. It is also a component of natural gas and can be used as a fuel."
    },
    "ferrocene": {
        "smiles": "c1cc[cH-]c1.[Fe]",
        "name": "Ferrocene",
        "formula": "C10H10Fe",
        "description": "Ferrocene is an organometallic compound with an iron atom sandwiched between two cyclopentadienyl rings. It is used in various industrial applications and organometallic chemistry."
    },
    "grignard reagent": {
        "smiles": "CC[Mg]Br",
        "name": "Ethylmagnesium bromide",
        "formula": "C2H5MgBr",
        "description": "Grignard reagents are organometallic compounds commonly used in organic synthesis to form carbon-carbon bonds."
    },
    "teflon": {
        "smiles": "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
        "name": "Teflon",
        "formula": "C2F4",
        "description": "Teflon, also known as polytetrafluoroethylene, is a polymer known for its non-stick properties and chemical resistance. It is used in cookware and industrial applications."
    },
    "bakelite": {
        "smiles": "c1cc(c(c(c1)OC(C)C)OC(C)C)OC(C)C",
        "name": "Bakelite",
        "formula": "C19H20O4",
        "description": "Bakelite is one of the first synthetic plastics, known for its heat resistance and durability. It is used in electronics and other industrial applications."
    },
    "haemoglobin": {
        "smiles": "CCC1=C2C(=CC=C1)C(=C2)C=O",
        "name": "Heme Group",
        "formula": "C34H32FeN4O4",
        "description": "Haemoglobin is a complex protein with an iron-containing heme group. It is responsible for transporting oxygen in the blood."
    },
    "benzene": {
        "name": "Benzene",
        "smiles": "c1ccccc1",
        "formula": "C6H6",
        "description": "Benzene is an aromatic hydrocarbon with a six-membered carbon ring. It is a common precursor in the chemical industry.",
    },
    "caffeine": {
        "name": "Caffeine",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "formula": "C8H10N4O2",
        "description": "Caffeine is a stimulant found in coffee, tea, and other beverages. It has a complex structure with multiple nitrogen atoms.",
    },
    "aspirin": {
        "name": "Aspirin",
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "formula": "C9H8O4",
        "description": "Aspirin, or acetylsalicylic acid, is a widely used pain reliever and anti-inflammatory drug.",
    },
    "ferrocene": {
        "name": "Ferrocene",
        "smiles": "c1cc[cH-]c1.[Fe]",
        "formula": "C10H10Fe",
        "description": "Ferrocene is an organometallic compound with iron at its core, commonly used in organometallic chemistry.",
    },
    "chlorophyll": {
        "name": "Chlorophyll",
        "smiles": "CC1=CC2=CC=C3C4=C5C6=C(C=C(C5)C4=NC3=NC2=C1)C=CC(C6)=C",
        "formula": "C55H72MgN4O5",
        "description": "Chlorophyll is a green pigment found in plants, essential for photosynthesis. It has a complex structure with magnesium at its center.",
    },
}

@app.route('/')
def home():
    return render_template("index.html")

@app.route('/generate_3d', methods=['POST'])
def generate_3d():
    data = request.get_json()
    compound_name = data.get("compound", "").strip().lower()

    # Get the compound details from the predefined mapping
    compound_details = compound_smiles_mapping.get(compound_name)

    if not compound_details:
        return jsonify({"success": False, "message": "Compound not found."})

    smiles = compound_details["smiles"]

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({"success": False, "message": "Invalid SMILES."})

        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)

        # Convert to SDF or another suitable format
        model_sdf = Chem.MolToMolBlock(mol)

        return jsonify({"success": True, "model": model_sdf, "name": compound_details["name"], "formula": compound_details["formula"],
        "description": compound_details["description"]})

    except Exception as e:
        return jsonify({"success": False, "message": str(e)})

if __name__ == "__main__":
    app.run(debug=True)
