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

with open('full database.xml') as file:
    is_within_drug = False
    for line in file:
        # Check if we're entering a new drug entry
        if '<drug ' in line:
            is_within_drug = True
            current_drug = {'name': None, 'cas_number': None, 'description': None}
        elif '</drug>' in line and is_within_drug:
            # End of the current drug entry
            if all(value is not None for value in current_drug.values()):  # Ensure all fields are populated
                names.append(current_drug['name'])
                cas_numbers.append(current_drug['cas_number'])
                descriptions.append(current_drug['description'])
            is_within_drug = False
        elif is_within_drug:
            # Extract data if we're within a drug entry
            if '<name>' in line and not current_drug['name']:
                current_drug['name'] = extract_data_from_line(line, 'name')
            elif '<cas-number>' in line and not current_drug['cas_number']:
                current_drug['cas_number'] = extract_data_from_line(line, 'cas-number')
            elif '<description>' in line and not current_drug['description']:
                current_drug['description'] = extract_data_from_line(line, 'description')

output_file_path = 'drugs_data.txt'

# Write the data to the file
with open(output_file_path, 'w') as file:
    for name, cas_number, description in zip(names, cas_numbers, descriptions):
        # Writing each drug's data followed by a delimiter line
        file.write(f"Name: {name}\n")
        file.write(f"CAS Number: {cas_number}\n")
        file.write(f"Description: {description}\n")
        file.write('-' * 50 + '\n')  # Delimiter line for readability
