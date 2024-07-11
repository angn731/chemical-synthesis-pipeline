# imports
import json
import os
from GroupByTemp import transform_data

def iterate_reactions(rxns, products):
    """
    Helper function that returns a set of the products
    produced in a list of rxn dictionaries.
    """
    update_products = set()
    # eaxh rxn is a dict with one k,v pair
    for rxn in rxns:
        for k in rxn.keys():
            potential_product = k.split('>>')[1]
            if potential_product in products:
                update_products.add(potential_product)
    return update_products


def set_products(wellplates, num, products):
    """
    Returns a set of the products produced in the
    first "num" wellplates of a sequence.
    """
    products_produced = set()
    for i in range(0, num):
        for k in wellplates.keys():
            if k.startswith(str(i)):
                rxns = wellplates[k]
                update_products = iterate_reactions(rxns, products)
                products_produced.update(update_products)
                break
    return products_produced


def indicate_products(wellplates, num, set_prods, rxns_to_pathways):
    """
    Takes as input a wellplate sequence, the "num" of wellplates
    to check, a set of the products produced in the specified
    wellplates, and a dict mapping rxn SMILES to a list of the rxn
    synthesis pathway it's in.
    Mutates the wellplates dictionary by adding a key to each reaction
    that specifies the final product and another key that specifies
    whether the reaction is the final reaction.
    """
    plate_ids = [n for n in range(num)]
    updates = 0
    final_prods = set()
    for plate, rxns in wellplates.items():
        if int(plate[0]) in plate_ids:
            # rxns is a list of dicts, each with 1 k,v pair
            for rxn in rxns:
                for rxn_smiles, conditions in rxn.items():
                    pathway = rxns_to_pathways[rxn_smiles]
                    conditions["pathway product"] = pathway[-1].split('>>')[1]
                    conditions["final product"] = False
                    # updates += 1
                    if rxn_smiles.split('>>')[1] in set_prods:
                        conditions["final product"] = True
                        final_prods.add(rxn_smiles.split('>>')[1])
                        updates += 1
    return wellplates, final_prods, updates


def identify_reactants(set_products, products_data):
    """
    Given a set of products and data about them, returns a set of necessary
    reactants.
    """
    reactants = set()
    print(f'products data: {products_data}')
    for product in set_products:
        data = products_data[product]
        templates = data['templates']
        # print(f'templates: {templates}')
        for template in templates:
            reactants.update(template['reactants'])
    return reactants


if __name__ == "__main__":
    pathways_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/MFBO_selected_mols_filtered_pathways.json'
    with open(pathways_filepath, 'r') as jsonfile:
        filtered_pathways = json.load(jsonfile)
    products = set(filtered_pathways.keys())

    rxns_top_conditions, all_rxns, rxns_to_pathways = transform_data(filtered_pathways)
    # print(rxns_to_pathways)

    wellplate_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/best_six_plate_seq/MFBO_selected_mols_six_plate_seq_85.json'
    with open(wellplate_filepath, 'r') as jsonfile:
        wellplate_seq = json.load(jsonfile)

    # check for number of products produced
    set_prods = set_products(wellplate_seq, 6, products)
    # print(f'set_prods: {set_prods}')

    # updated_wellplates, final_prods, updates = indicate_products(wellplate_seq, 6, set_prods, rxns_to_pathways)

    # save mutated dict
    # dir = "/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols"
    # filename = f"MFBO_selected_mols_six_plate_seq_product_key.json"
    # filepath = os.path.join(dir, filename)
    # with open(filepath, 'w') as outfile:
    #     json.dump(updated_wellplates, outfile, indent=4)
    # print(wellplate_seq)

    # save a set of the necessary reactants
    products_data_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/MFBO_selected_mols_for_synthesis.json'
    with open(products_data_filepath, 'r') as jsonfile:
        products_data = json.load(jsonfile)
    # print(f'products_data: {products_data}')

    dir = "/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/best_six_plate_seq"
    filename = f"six_plate_reactants.json"
    filepath = os.path.join(dir, filename)
    reactants = identify_reactants(set_prods, products_data)
    # print(products_data['CS(=O)(=O)N(S(=O)(=O)c1ccc(C=Cc2ccccc2)cc1)S(=O)(=O)c1cccs1'])
    # with open(filepath, 'w') as outfile:
    #     json.dump(list(reactants), outfile, indent=4)

    # small test case
    # products_data = {
    # "CS(=O)(=O)N(S(=O)(=O)c1ccc(C=Cc2ccccc2)cc1)S(=O)(=O)c1cccs1": {
    #     "HDAC_Docking": ["0.003787613349295299", "0.15"],
    #     "LogP": ["4.156231893789057", "0.33392553915015455"],
    #     "Toxicity": ["0.6682438342203894", "0.04982561078733472"],
    #     "similarity": ["0.49042339295468546", 0],
    #     "rank": 0,
    #     "initial_mol": "CS(N)(=O)=O",
    #     "templates": [
    #     {"_id": "Sulfonamide", "reactants": ["O=S(=O)(Cl)c1ccc(Br)cc1"]},
    #     {"_id": "Sulfonamide", "reactants": ["O=S(=O)(Cl)c1cccs1"]},
    #     {"_id": "Suzuki_vinyl_aryl", "reactants": ["OB(O)/C=C/c1ccccc1"]}
    #     ]
    # }}

    # set_prods = {"CS(=O)(=O)N(S(=O)(=O)c1ccc(C=Cc2ccccc2)cc1)S(=O)(=O)c1cccs1"}

    # print(identify_reactants(set_prods, products_data))
