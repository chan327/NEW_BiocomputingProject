from flask import Flask, request, jsonify, render_template

app = Flask(__name__)


def extract_data_from_line(line, tag):
    # This function extracts data enclosed within the specified XML tags in the line
    start_tag = f"<{tag}>"
    end_tag = f"</{tag}>"
    start_index = line.find(start_tag) + len(start_tag)
    end_index = line.find(end_tag)
    if start_index > -1 and end_index > -1:
        return line[start_index:end_index].strip()
    return None


# Initialize variables to hold the extracted information in separate lists
names = []
cas_numbers = []
descriptions = []
smiles = []

with open('full database.xml') as file:
    is_within_drug = False
    for line in file:
        # Check if we're entering a new drug entry
        if '<drug ' in line:
            is_within_drug = True
            current_drug = {'name': None, 'cas_number': None, 'description': None, 'smiles': None}
        elif '</drug>' in line and is_within_drug:
            # End of the current drug entry
            if all(value is not None for value in current_drug.values()):  # Ensure all fields are populated
                names.append(current_drug['name'])
                cas_numbers.append(current_drug['cas_number'])
                descriptions.append(current_drug['description'])
                smiles.append(current_drug['smiles'])
            is_within_drug = False
        elif is_within_drug:
            # Extract data if we're within a drug entry
            if '<name>' in line and not current_drug['name']:
                current_drug['name'] = extract_data_from_line(line, 'name')
            elif '<cas-number>' in line and not current_drug['cas_number']:
                current_drug['cas_number'] = extract_data_from_line(line, 'cas-number')
            elif '<description>' in line and not current_drug['description']:
                current_drug['description'] = extract_data_from_line(line, 'description')
            elif '<kind>SMILES</kind>' in line and not current_drug['smiles']:
                value = next(file)
                current_drug['smiles'] = extract_data_from_line(value, 'value')


@app.route('/')
def index():
    return render_template('index.html', name='PyCharm')

@app.route('/search', methods=['GET'])
def search():
    query = request.args.get('query', '').lower()
    results = [names for names in names if query in names.lower()]
    return jsonify(results[:10])

@app.route('/item/<name>')
def item_page(name):
    if name in names:
        return render_template('results.html', item_name=name)
    else:
        return "Item not found", 404

if __name__ == '__main__':
    app.run(debug=True, port=8001)


