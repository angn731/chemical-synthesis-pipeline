# imports
import json
from GroupByTemp import transform_data2

filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/best_six_plate_seq/MFBO_selected_mols_six_plate_seq_85.json'
with open(filepath, 'r') as jsonfile:
    wellplates = json.load(jsonfile)

filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/best_six_plate_seq/six_plate_reactants.json'
with open(filepath, 'r') as jsonfile:
    reagents = json.load(jsonfile)

current_reagents = set(reagents)

plate_four_rxns = []
plate_three_rxns = []
plate_two_rxns = []
plate_one_rxns = []
plate_zero_rxns = []

for k,v in wellplates.items():
    if k.startswith('4_'):
        plate_four_rxns = v
    if k.startswith('3_'):
        plate_three_rxns = v
    if k.startswith('2_'):
        plate_two_rxns = v
    if k.startswith('1_'):
        plate_one_rxns = v
    if k.startswith('0_'):
        plate_zero_rxns = v


# print(f'reactants: {reactants}')
# print(f'plate 1: {plate_one_rxns}')

def set_products(prev_plate):
    """
    Takes the previous plate, which is a list of dicts of rxns,
    and returns a set of the products produced.
    """
    prods = set()
    for rxn in prev_plate:
        for rxn_smiles in rxn.keys():
            prod = rxn_smiles.split('>>')[1]
            prods.add(prod)
    return prods

def get_reactants(current_plate):
    """
    Takes the current plate, which is a list of dicts of rxns,
    and returns a set of the reactants needed.
    """
    reacts = set()
    for rxn in current_plate:
        for rxn_smiles in rxn.keys():
            rxn_reactants = rxn_smiles.split('>>')[0].split('.')
            reacts.update(rxn_reactants)
    return reacts

def get_all_reactants(plates):
    """
    Takes a list of plates and returns all the reactants needed.
    """
    reacts = set()
    for plate in plates:
        plate_reactants = get_reactants(plate)
        reacts.update(plate_reactants)
    return reacts


if __name__ == "__main__":

    # all the intermediates/final prods produced in the first 4 plates
    prev_products = set_products(plate_three_rxns) | set_products(plate_two_rxns) | set_products(plate_one_rxns) | set_products(plate_zero_rxns)
    # set of all the current available chemicals
    avail_reactants = prev_products | current_reagents
    # set of all the reactants needed
    needed_reactants = get_all_reactants([plate_four_rxns, plate_three_rxns, plate_two_rxns, plate_one_rxns, plate_zero_rxns])

    results = needed_reactants - avail_reactants

    # results are below -- these are the missing chemicals
    # {'Nc1c(Cl)cccc1Cl', 'N#Cc1ccccc1B(O)O', 'B(C1=CSC=C1)(O)O',
    #  'CC1=CC(=C(C(=C1)C)N)C', 'O=C(O)c1cccc(F)c1C(=O)O', 'Nc1cccnc1N',
    #  'OB(O)c1ccc(F)cc1', 'C=CO', 'O=[N+]([O-])c1ccc(F)cc1', 'Nc1ncc([N+](=O)[O-])cn1',
    #  'NS(=O)(=O)O', 'C1=CN=C(N=C1)N', 'Nc1ncn[nH]1', 'N#CN', 'Brc1ccccc1', 'CC(C)(C)O', 'Clc1ccccc1Cl'}

    print(results)

    # inp_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/new_rxns_to_pathways.json'
    # with open(inp_filepath, 'r') as jsonfile:
    #     rxns_to_pathways = json.load(jsonfile)
