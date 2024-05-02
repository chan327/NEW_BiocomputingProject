from flask import Flask, request, jsonify, render_template, url_for
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from io import BytesIO
import os

app = Flask(__name__)

# Function to extract data from XML line
def extract_data_from_line(line, tag):
    start_tag = f"<{tag}>"
    end_tag = f"</{tag}>"
    start_index = line.find(start_tag) + len(start_tag)
    end_index = line.find(end_tag)
    if start_index > -1 and end_index > -1:
        return line[start_index:end_index].strip()
    return None

# Initialize lists to hold extracted information
names = []
cas_numbers = []
descriptions = []
smiles = []

# Parse XML file
with open('full database.xml', encoding='UTF-8') as file:
    is_within_drug = False
    for line in file:
        if '<drug ' in line:
            is_within_drug = True
            current_drug = {'name': None, 'cas_number': None, 'description': None, 'smiles': None}
        elif '</drug>' in line and is_within_drug:
            if all(value is not None for value in current_drug.values()):
                names.append(current_drug['name'])
                cas_numbers.append(current_drug['cas_number'])
                descriptions.append(current_drug['description'])
                smiles.append(current_drug['smiles'])
            is_within_drug = False
        elif is_within_drug:
            if '<name>' in line and not current_drug['name']:
                current_drug['name'] = extract_data_from_line(line, 'name')
            elif '<cas-number>' in line and not current_drug['cas_number']:
                current_drug['cas_number'] = extract_data_from_line(line, 'cas-number')
            elif '<description>' in line and not current_drug['description']:
                current_drug['description'] = extract_data_from_line(line, 'description')
            elif '<kind>SMILES</kind>' in line and not current_drug['smiles']:
                value = next(file)
                current_drug['smiles'] = extract_data_from_line(value, 'value')

# Create dictionary of drugs
drugs_dict = {}
for i in range(len(names)):
    drugs_dict[names[i]] = {
        'cas_number': cas_numbers[i],
        'description': descriptions[i],
        'smiles': smiles[i]
    }

# Function to generate molecule image
def generate_molecule_image(smiles_code):
    mol = Chem.MolFromSmiles(smiles_code)
    if mol is None:
        return None
    img = Draw.MolToImage(mol)
    img_bytes = BytesIO()
    img.save(img_bytes, format='JPEG')
    img_bytes.seek(0)
    return img_bytes

# Function to compare two smiles strings
def two_smiles(smiles1, smiles2):
    """-----------------------------------------------------------------------------------------
    Compare the SMILES strings for two different drugs. Determine the Tanimoto similarity scores for two types of
    fingerprints: Morgan fingerprints and substructure fingerprints.

    :param smiles1: string, for the inputted drug of interest
    :param smiles2: string, for the drug in the comprehensive set that will be compared with the inputted drug
    :return: sims, dictionary containing two key-value pairs (tanimoto_morgan and tanimoto_sub)
    -----------------------------------------------------------------------------------------"""

    sims = {}
    # Load molecules from SMILES strings
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    # Check if both molecules were successfully created
    if mol1 is None or mol2 is None:
        # One or both of the SMILES strings can't be parsed, return None
        return None

    try:
        # Aromatize molecules
        mol1 = Chem.Mol(mol1)
        mol2 = Chem.Mol(mol2)
        Chem.SanitizeMol(mol1)
        Chem.SanitizeMol(mol2)
        Chem.Kekulize(mol1)
        Chem.Kekulize(mol2)

        # Calculate similarity between Morgan fingerprints, radius of 2
        morgan1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
        morgan2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)

        tanimoto_morgan = DataStructs.TanimotoSimilarity(morgan1, morgan2)

        # Calculate similarity between Substructure fingerprints
        substrcture1 = Chem.PatternFingerprint(mol1)
        substructure2 = Chem.PatternFingerprint(mol2)

        tanimoto_sub = DataStructs.TanimotoSimilarity(substrcture1, substructure2)
        sims['tanimoto_morgan'] = tanimoto_morgan
        sims['tanimoto_sub'] = tanimoto_sub

    # Exception handling
    except Exception as e:
        print(f"Error occurred: {e}")
        return None

    return sims


# Function to determine SMILES similarity among all drugs in the comprehensive dataset
def smiles_similarity(input, drug_dict):
    """-----------------------------------------------------------------------------------------
    Compile the similarity scores of the inputted drug's SMILES string with that of all other drugs in the dataset

    :param input: string, containing the SMILES string for the drug of interest
    :param drug_dict: dict, containing information for all drugs
    :return: similarity_list, list of dictionaries containing the similarity scores for each drug
    -----------------------------------------------------------------------------------------"""
    similarity_list = []
    # Iterate over each drug in the drug data
    for key, value in drug_dict.items():
        # Check if the drug contains a SMILES string
        if isinstance(value, dict) and 'smiles' in value:
            similarity_dict = {}
            # Calculate similarity score between input SMILES and drug's SMILES
            try:
                similarity_score = two_smiles(input, value['smiles'])
                # Store similarity score for each drug
                if similarity_score is not None:
                    similarity_dict[key] = similarity_score
                    similarity_list.append(similarity_dict)
            except Exception as e:
                # print(f"Error occurred for drug '{key}': {e}")
                pass
    return similarity_list

# Function to find the top 10 most similar drugs by SMILES
def most_similar_drugs(similarity_scores):
    """-----------------------------------------------------------------------------------------
    Sort the smiles_smiliarity() function output to return the top 10 most similar drugs to the inputted drug

    :param similarity_scores: list of dictionaries containing the similarity scores for each drug
    :return: ten_drugs_list: list of 10 dictionaries with information about the most similar drugs
    -----------------------------------------------------------------------------------------"""

    # Combine the list of dictionaries into a dictionary
    scores = {key: value for dict in similarity_scores for key, value in dict.items()}
    # Sort drugs by similarity score in descending order
    # Remove the first item in the list because the most similar drug in the list to the inputted drug will be the
    # inputted drug itself
    sorted_drugs = sorted(scores.items(), key=lambda x: max(x[1].values()), reverse=True)[1:]
    # Store the top 10 most similar drugs and their similarity scores
    top_ten = sorted_drugs[:10]
    ten_drugs_list = []
    for i, (drug, score) in enumerate(top_ten[:10], start=1):
        drug_dict = {
            "rank": i,
            "drug_name": drug,
            "morgan_similarity": score['tanimoto_morgan'],
            "substructure_similarity": score['tanimoto_sub']
        }
        ten_drugs_list.append(drug_dict)

    return ten_drugs_list


@app.route('/')
def index():
    return render_template('index.html', name='PyCharm')

@app.route('/item/<name>')
def item_page(name):
    drug = drugs_dict.get(name)
    if drug:
        # Define variables relevant to SMILES comparison
        smiles_code = drug['smiles']
        smiles_similar_drugs = smiles_similarity(smiles_code, drugs_dict)
        smiles_ten = most_similar_drugs(smiles_similar_drugs)
        # Creating molecule image
        img_bytes = generate_molecule_image(smiles_code)
        if img_bytes:
            img_bytes.seek(0)
            temp_dir = os.path.join(app.root_path, 'static', 'temp')
            if not os.path.exists(temp_dir):
                os.makedirs(temp_dir)
            image_path = os.path.join(temp_dir, f'{name}.jpg')
            with open(image_path, 'wb') as img_file:
                img_file.write(img_bytes.read())
            image_url = url_for('static', filename=f'temp/{name}.jpg')
            return render_template('results.html', item_name=name, drug=drug, molecule_image=image_url,
                                   similar_drugs=smiles_ten)
    return "Item not found", 404

@app.route('/search', methods=['GET'])
def search():
    query = request.args.get('query', '').lower()
    results = [name for name in names if query in name.lower()]
    return jsonify(results[:10])

if __name__ == '__main__':
    app.run(debug=True, port=8001)

