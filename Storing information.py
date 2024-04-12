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
list_of_drugs = []

with open('full database.xml') as file:
    is_within_drug = False
    for line in file:
        # Check if we're entering a new drug entry
        if '<drug ' in line:
            is_within_drug = True
            current_drug = {'name': None, 'cas_number': None, 'description': None, 'smiles': None}
            list_of_drugs.append(current_drug)
        elif '</drug>' in line and is_within_drug:
            # End of the current drug entry
            if current_drug['name'] is not None and current_drug['cas_number'] is not None and current_drug['description'] is not None:
                # Ensure the fields for name, cas, and description are populated
                names.append(current_drug['name'])
                cas_numbers.append(current_drug['cas_number'])
                descriptions.append(current_drug['description'])
                if current_drug['smiles'] is not None:
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
# print(smiles)
# print(len(names), len(cas_numbers), len(descriptions), len(smiles))
print(list_of_drugs[20])
