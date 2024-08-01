# imports
import json
import os
import GroupByTemp

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
            if k.startswith(f'{str(i)}_'):
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
        if int(plate[0]) in plate_ids and plate[1] == '_':
            # rxns is a list of dicts, each with 1 k,v pair
            for rxn in rxns:
                # there is only 1 rxn_smiles and conditions pair
                for rxn_smiles, conditions in rxn.items():
                    # print(rxn_smiles)
                    pathway = rxns_to_pathways[rxn_smiles]
                    conditions["pathway product"] = pathway[-1].split('>>')[1]
                    conditions["final product"] = False
                    # updates += 1
                    if rxn_smiles.split('>>')[1] == conditions["pathway product"]:
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
    # print(f'products data: {products_data}')
    for product in set_products:
        data = products_data[product]
        templates = data['templates']
        # print(f'templates: {templates}')
        for template in templates:
            reactants.update(template['reactants'])
    return reactants


if __name__ == "__main__":
    pathways_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/mols_iter2/iter2_filtered_pathways.json'
    with open(pathways_filepath, 'r') as jsonfile:
        filtered_pathways = json.load(jsonfile)
    products = set(filtered_pathways.keys())

    wellplate_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/mols_iter2/iter2_best_wellplate_seq.json'
    with open(wellplate_filepath, 'r') as jsonfile:
        wellplate_seq = json.load(jsonfile)

    rxns_to_pathways_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/mols_iter2/iter2_tagged_rxns_to_pathways.json'
    with open(rxns_to_pathways_filepath, 'r') as jsonfile:
        rxns_to_pathways = json.load(jsonfile)

    # check for number of products produced
    set_prods = set_products(wellplate_seq, 4, products)
    # print(len(set_prods))

    updated_wellplates, final_prods, updates = indicate_products(wellplate_seq, 4, set_prods, rxns_to_pathways)
    # print(len(final_prods))

    # save mutated dict
    dir = "/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/mols_iter2"
    filename = f"iter2_indicate_products.json"
    filepath = os.path.join(dir, filename)
    with open(filepath, 'w') as outfile:
        json.dump(updated_wellplates, outfile, indent=4)

    # save a set of the necessary reactants
    products_data_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/MFBO_selected_mols_for_synthesis.json'
    with open(products_data_filepath, 'r') as jsonfile:
        products_data = json.load(jsonfile)
    # print(f'products_data: {products_data}')
