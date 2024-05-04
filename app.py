from flask import Flask, request, jsonify, render_template, url_for
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import MolToImage
from io import BytesIO
import os
import pubchempy as pcp
import html


app = Flask(__name__)

# Function to extract data from XML line
def extract_data_from_line(line, tag):
    """-----------------------------------------------------------------------------------------
   Extract relevant data from the XML file by recognizing information stored within tags

    :param line: line of XML file
    :param tag: string, specifying what will be extracted
    :return: string containing relevant information
    -----------------------------------------------------------------------------------------"""
    start_tag = f"<{tag}>"
    end_tag = f"</{tag}>"
    start_index = line.find(start_tag) + len(start_tag)
    end_index = line.find(end_tag)
    # Start and end tag on the same line
    if start_index > -1 and end_index > -1:
        return line[start_index:end_index].strip()

    # Start tag on one line, no end tag
    elif start_index > -1 and end_index == -1:
        data = line[start_index:].strip()
        for next_line in file:
            end_index = next_line.find(end_tag)
            if end_index > -1:
                data += " " + next_line[:end_index].strip()
                break
            else:
                data += " " + next_line.strip()
        return data
    else:
        return None


def process_text(input_text):
    """-----------------------------------------------------------------------------------------
    Process the strings returned after passing through html.unescape() so that newlines are preserved, and text
    in square brackets is removed (originally these were hyperlinks in the DrugBank website)

    :param input_text: string, containing unescaped text
    :return: processed_text, string containing html tags
    -----------------------------------------------------------------------------------------"""
    # Remove the parts of the string that fall within square brackets
    new_text = input_text.split('[')
    new_text_list = []
    for part in new_text:
        end_bracket_index = part.find(']')
        if end_bracket_index != -1:
            new_text_list.append(part[end_bracket_index + 1:])
        else:
            new_text_list.append(part)
    # Join the processed parts back into a string
    new_text_list = ''.join(new_text_list)

    # Replace '**' with a newline indicated by <br>
    parts = new_text_list.split('**')
    new_text_list = []
    for i, part in enumerate(parts):
        if i != 0:
            # Add newline to everything but the first part
            new_text_list.append('\n')
        new_text_list.append(part.strip())
    # Join the parts back into a string, with <br> as an HTML tag
    processed_text = '<br>'.join(new_text_list)

    return processed_text

# Initialize lists to hold extracted drug information
names = []
cas_numbers = []
descriptions = []
smiles = []
pharmacodynamics = []
mechanism = []
toxicity = []

# Parse XML file
with open('full database.xml', encoding='UTF-8') as file:
    is_within_drug = False
    for line in file:
        if '<drug ' in line:
            is_within_drug = True
            # Create initial dictionary to store information
            current_drug = {'name': None, 'cas_number': None, 'description': None, 'smiles': None, 'pharmacodynamics':
                            None, 'mechanism': None, 'toxicity': None}
        elif '</drug>' in line and is_within_drug:
            if all(value is not None for value in current_drug.values()):
                names.append(current_drug['name'])
                cas_numbers.append(current_drug['cas_number'])
                descriptions.append(current_drug['description'])
                smiles.append(current_drug['smiles'])
                pharmacodynamics.append(current_drug['pharmacodynamics'])
                mechanism.append(current_drug['mechanism'])
                toxicity.append(current_drug['toxicity'])
            is_within_drug = False
        elif is_within_drug:
            # Extract information based on inputted tags and store as dictionary values
            if '<name>' in line and not current_drug['name']:
                current_drug['name'] = extract_data_from_line(line, 'name')
            elif '<cas-number>' in line and not current_drug['cas_number']:
                current_drug['cas_number'] = extract_data_from_line(line, 'cas-number')
            elif '<description>' in line and not current_drug['description']:
                original_desc = extract_data_from_line(line, 'description')
                original_desc = html.unescape(original_desc)
                current_drug['description'] = process_text(original_desc)
            elif '<kind>SMILES</kind>' in line and not current_drug['smiles']:
                value = next(file)
                current_drug['smiles'] = extract_data_from_line(value, 'value')
            elif '<pharmacodynamics>' in line and not current_drug['pharmacodynamics']:
                original_pharm = extract_data_from_line(line, 'pharmacodynamics')
                original_pharm = html.unescape(original_pharm)
                current_drug['pharmacodynamics'] = process_text(original_pharm)
            elif '<mechanism-of-action>' in line and not current_drug['mechanism']:
                original_mech = extract_data_from_line(line, 'mechanism-of-action')
                original_mech= html.unescape(original_mech)
                current_drug['mechanism'] = process_text(original_mech)
            elif '<toxicity>' in line and not current_drug['toxicity']:
                original_tox = extract_data_from_line(line, 'toxicity')
                original_tox = html.unescape(original_tox)
                current_drug['toxicity'] = process_text(original_tox)

# Create dictionary of drugs
drugs_dict = {}
for i in range(len(names)):
    drugs_dict[names[i]] = {
        'cas_number': cas_numbers[i],
        'description': descriptions[i],
        'smiles': smiles[i],
        'pharmacodynamics': pharmacodynamics[i],
        'mechanism': mechanism[i],
        'toxicity': toxicity[i]
    }

# Function to generate molecule image
def generate_molecule_image(smiles_code):
    """-----------------------------------------------------------------------------------------
   Generate a molecule image from a drug's SMILES string

    :param smiles_code: string, smiles string stored in drug dictionary
    :return: two JPEG images, one of 2D construction of molecule and another with stereochemistry
    -----------------------------------------------------------------------------------------"""
    mol = Chem.MolFromSmiles(smiles_code)
    if mol is None:
        return None
    # Generate 2D structure
    img = Draw.MolToImage(mol)
    img_bytes = BytesIO()
    # Save image as JPEG
    img.save(img_bytes, format='JPEG')
    img_bytes.seek(0)

    # Create another structure, using 3D coordinates to determine stereochemistry
    mol_plus = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_plus)
    AllChem.MMFFOptimizeMolecule(mol_plus)

    # Convert these coordinates back into an image
    img_plus = MolToImage(mol_plus)
    img_bytes_plus = BytesIO()
    # Save image as JPEG
    img_plus.save(img_bytes_plus, format='JPEG')
    img_bytes_plus.seek(0)

    return img_bytes, img_bytes_plus
    # For some drugs, there was a ValueError: Bad Conformer ID when using AllChem.MMFFOptimizeMolecule(mol_plus). To
    # generate only the 2D structure, uncomment out the line below.
    # return img_bytes



""" This is the work we completed trying to use PyMol to generate a 3D representation of the molecule. We weren't able
to get this part completed successfully, but here was our progress. 

    def generate_3d_molecule(smiles_code):
    
    Generate a 3D representation of a molecule from its SMILES code using Py3Dmol.

    :param smiles_code: SMILES code of the molecule.
    :return: HTML string containing the 3D visualization of the molecule.
    
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(smiles_code, 'sdf')
    viewer.setStyle({'stick': {}})
    viewer.setBackgroundColor('white')
    viewer.zoomTo()
    return viewer.render()
Code leads to SMILES parse Error and is showing up as a blank box on the HTML file, unsure of how to fix
"""


# Function to compare two smiles strings
def two_smiles(smiles1, smiles2):
    """-----------------------------------------------------------------------------------------
    Compare the SMILES strings for two different drugs. Determine the Tanimoto similarity scores for two types of
    fingerprints: Morgan fingerprints and substructure fingerprints.
    Websites referenced for help:
    https://ctr.fandom.com/wiki/Report_the_similarity_between_two_structures#Indigo/Python
    https://www.researchgate.net/post/How_to_compare_a_set_of_SMILES_structures

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
            # Exception handling, program continues execution
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

def gather_pubchem_info(input):
    """-----------------------------------------------------------------------------------------
    Use the pubchempy library to parse through the PubChem of the inputted drug
    https://pubchempy.readthedocs.io/en/latest/guide/properties.html

    :param input: string, name of drug of interest
    :return: pubchem_dict: dictionary containing 13 key-value pairs with chemical information
    -----------------------------------------------------------------------------------------"""

    pubchem_results = pcp.get_properties(('MolecularFormula', 'MolecularWeight',
                               'InChI', 'InChIKey', 'IUPACName', 'TPSA', 'Complexity', 'Charge', 'XLogP',
                                'ExactMass', 'MonoisotopicMass', 'HBondDonorCount', 'HBondAcceptorCount'),
                              input, 'name')

    pubchem_dict = pubchem_results[0]
    return pubchem_dict


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

        # Define variable for PubChem information
        pubchem = gather_pubchem_info(name)

        # Creating molecule images
        # Remove img_bytes_plus if the search query returns a ValueError, to access the page with only one image
        # img_bytes = generate_molecule_image(smiles_code)
        img_bytes, img_bytes_plus = generate_molecule_image(smiles_code)
        # if img_bytes:
        if img_bytes and img_bytes_plus:
            img_bytes.seek(0)
            temp_dir = os.path.join(app.root_path, 'static', 'temp')
            if not os.path.exists(temp_dir):
                os.makedirs(temp_dir)
            image_path = os.path.join(temp_dir, f'{name}.jpg')
            with open(image_path, 'wb') as img_file:
                img_file.write(img_bytes.read())

            # Save 3D image
            # Comment out this section to access the web page with only one image
            image_path_plus = os.path.join(temp_dir, f'{name}_plus.jpg')
            with open(image_path_plus, 'wb') as img_file:
                img_file.write(img_bytes_plus.read())

            image_url = url_for('static', filename=f'temp/{name}.jpg')
            # Comment out the line below to access the web page with only one image
            image_url_plus = url_for('static', filename=f'temp/{name}_3d.jpg')

            return render_template('results.html', item_name=name, drug=drug, molecule_image=image_url,
                                   molecule_image_plus=image_url_plus, similar_drugs=smiles_ten, pubchem=pubchem,
                                   pubmed_search = f"https://pubmed.ncbi.nlm.nih.gov/?term={name}",
                                   google_scholar = f"https://scholar.google.com/scholar?hl=en&as_sdt=0%2C15&q={name}&oq=)",
                                   science_direct = f"https://www.sciencedirect.com/search?qs={name}")
        # return render_template('results.html', item_name=name, drug=drug, molecule_image=image_url,
        #                        similar_drugs=smiles_ten, pubchem=pubchem,
        #                        pubmed_search=f"https://pubmed.ncbi.nlm.nih.gov/?term={name}",
        #                        google_scholar=f"https://scholar.google.com/scholar?hl=en&as_sdt=0%2C15&q={name}&oq=)",
        #                        science_direct=f"https://www.sciencedirect.com/search?qs={name}")
    return "Item not found", 404

@app.route('/search', methods=['GET'])
def search():
    query = request.args.get('query', '').lower()
    results = [name for name in names if query in name.lower()]
    return jsonify(results[:10])

if __name__ == '__main__':
    app.run(debug=True, port=8001)


