import csv
import json

product_smiles = []

with open('/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/mols_iter2/single_pt_mols_iter2.csv', newline='') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    next(csvreader)  # Skip the header row if there is one
    for row in csvreader:
        product_smiles.append(row[1])


pathways_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/gen_mols/synthesis_pathways.json'
with open(pathways_filepath, 'r') as file:
    all_pathways = json.load(file)

reaction_pathways = {}
not_found = set()
for product in product_smiles:
    try:
        reaction_pathways[product] = all_pathways[product]
    except:
        not_found.add(product)
        continue

iter2_pathways_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/mols_iter2/iter2_synthesis_pathways.json'
with open(iter2_pathways_filepath, 'w') as outfile:
    json.dump(reaction_pathways, outfile, indent=4)
