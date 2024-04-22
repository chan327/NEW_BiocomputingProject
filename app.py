from flask import Flask, request, jsonify, render_template, url_for
from rdkit import Chem
from rdkit.Chem import Draw
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

@app.route('/')
def index():
    return render_template('index.html', name='PyCharm')

@app.route('/item/<name>')
def item_page(name):
    drug = drugs_dict.get(name)
    if drug:
        smiles_code = drug['smiles']
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
            return render_template('results.html', item_name=name, drug=drug, molecule_image=image_url)
    return "Item not found", 404

@app.route('/search', methods=['GET'])
def search():
    query = request.args.get('query', '').lower()
    results = [name for name in names if query in name.lower()]
    return jsonify(results[:10])

if __name__ == '__main__':
    app.run(debug=True, port=8001)