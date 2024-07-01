from rdkit import Chem
from rdkit.Chem import AllChem
import json
import ForwardPred
import os

# output: dictionary, where keys are SMILES strings representing reaction products and the values are tuples
# each element of the tuple is a dictionary produced by the forward prediction code.
# so first element would be the first reaction, second element would be the second reaction, etc.

templates_dict_path = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/templates_dict.json'

with open(templates_dict_path, 'r') as file:
    rxn_templates = json.load(file)


def react(reactant1_smiles, reactant2_smiles, template_id):
    """
    Helper function that returns a SMILES of the product,
    given two reactant SMILES and a reaction template ID.
    """
    reactant1 = Chem.MolFromSmiles(reactant1_smiles)
    reactant2 = Chem.MolFromSmiles(reactant2_smiles)
    reaction_comps = rxn_templates[template_id].split(">>")
    reaction_smarts = reaction_comps[1] + ">>" + reaction_comps[0]
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    products = rxn.RunReactants((reactant1, reactant2))
    product_smiles = ""
    if products:
        product = products[0][0]  # Assuming one product
        product_smiles = Chem.MolToSmiles(product)
    else:
        products = rxn.RunReactants((reactant2, reactant1))
        try:
            product = products[0][0]
            product_smiles = Chem.MolToSmiles(product)
        # if neither reactant order works
        except IndexError:
            return None
    return product_smiles

def react_one_reactant(reactant_smiles, template_id):
    """
    Helper function that returns a SMILES of the product,
    given one reactant SMILES and a reaction template ID.
    """
    reactant = Chem.MolFromSmiles(reactant_smiles)
    reaction_comps = rxn_templates[template_id].split(">>")
    reaction_smarts = reaction_comps[1] + ">>" + reaction_comps[0]
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    products = rxn.RunReactants((reactant, ))
    product_smiles = ""
    try:
        product = products[0][0]  # Assuming one product
        product_smiles = Chem.MolToSmiles(product)
    except IndexError:
        return None
    return product_smiles

def make_rxn_smiles(rxn_molecules):
    """
    Helper function that takes in a tuple of SMILES and returns a reaction SMILES.
    Compatible with one or two reactants.
    """
    length = len(rxn_molecules)
    if length == 3:
        return rxn_molecules[0] + '.' + rxn_molecules[1] + ">>" + rxn_molecules[2]
    if length == 2:
        return rxn_molecules[0] + ">>" + rxn_molecules[1]

# templates is a list of dictionaries, each with two keys
# _id is a string
# reactants is a list, where we can assume that there's only one element
# next step is given a reaction product (a key in the top-level dict), generate a tuple of reaction SMILES strings
# can feed each of these strings in forward pred in a different function

def synthesis_pathway(product, data):
    """
    Helper function that given a SMILES string of the product and a dictionary
    of associated data, returns a tuple of reactions needed to synthesize the product.
    Reactions are represented as SMILES strings.
    """
    print(f'desired product: {product}')
    print(f'data: {data}')
    templates = data["templates"]
    reactant1 = data["initial_mol"]
    rxn_pathway = ()
    # loop through all the reaction templates, where each template is a dict
    for idx, temp in enumerate(templates):
        template_id = temp["_id"]
        # assume there is only one reactant
        rxn_molecules = tuple()
        multiple_prods = False
        one_reactant = False
        try:
            reactant2 = temp["reactants"][0]
            intermediate = react(reactant1, reactant2, template_id)
            rxn_molecules = (reactant1, reactant2, intermediate)
        # the intermediate from the previous step is the only reactant
        except IndexError:
            one_reactant = True
            intermediate = react_one_reactant(reactant1, template_id)
            rxn_molecules = (reactant1, intermediate)
        # if one of the reactions in the sequence of reactions doesn't work, return immediately
        if intermediate is None:
            return None
        rxn_smiles = make_rxn_smiles(rxn_molecules)
        rxn_pathway += (rxn_smiles,)
        # update reactant 1 to react the intermediate by the previous step
        # with the reactant in the next step
        reactant1 = intermediate
        # if isinstance(intermediate, list):
        #     matches_product = False
        #     for possible_prod in intermediate:
        #         if Chem.MolToSmiles(Chem.MolFromSmiles(product)) == possible_prod:
        #             matches_product = True
        #             if one_reactant:
        #                 rxn_molecules = (reactant1, possible_prod)
        #                 rxn_smiles = make_rxn_smiles(rxn_molecules)
        #                 rxn_pathway += (rxn_smiles,)
        #             else:
        #                 rxn_molecules = (reactant1, reactant2, possible_prod)
        #                 rxn_smiles = make_rxn_smiles(rxn_molecules)
        #                 rxn_pathway += (rxn_smiles,)
        #     if not matches_product:
        #         raise ValueError("The final product does not match the desired product.")
    # the final product should be the desired product
    if Chem.MolToSmiles(Chem.MolFromSmiles(reactant1)) != product:
        print(f'pathway: {rxn_pathway}')
        print(f'produced product: {reactant1}')
        raise ValueError("The final product does not match the desired product.")
    return rxn_pathway

# forward pred took a json file that was a list of reaction SMILES strings, then output a dict
# where the keys were reaction SMILES strings and the values were a list of dicts with conditions and
# the probability of generating the desired product


# ////////////////////////////////////////////
# GENERATING SYNTHESIS PATHWAYS AND CONDITIONS
# ////////////////////////////////////////////

def read_file(file_path):
    """
    Helper function that reads in a json file mapping a molecule to data like
    drug activity, similarity score, LogP, toxicity, and templates (to determine
    synthesis pathway).
    """
    with open(file_path, 'r') as jsonfile:
        products_data = json.load(jsonfile)
    return products_data


def save_file(file_path, data):
    """
    Helper function that saves data into a json file.
    """
    with open(file_path, 'w') as outfile:
        json.dump(data, outfile, indent=4)
    print(f"File saved to {file_path}")


def make_synthesis_pathways(products_data):
    """
    Helper function that returns a 2-elem tuple, where the
    1st elem is a list of all the reactions and the
    2nd elem is a dictionary mapping each product SMILES to a tuple
    of ordered reactions SMILES necessary to synthesize the product.
    """
    # store synthesis pathways in a dict
    reaction_pathways = {}
    # simultaneously store a list of all the reactions
    all_reactions = []
    for prod, data in products_data.items():
        # rxn_pathway is a tuple of reaction strings
        pathway = synthesis_pathway(prod, data)
        # if the returned pathway is None, one of the reactions didn't
        # produce a valid result, so skip the product entirely
        if pathway is None:
            continue
        reaction_pathways[prod] = synthesis_pathway(prod, data)
        rxns = list(reaction_pathways[prod])
        all_reactions.extend(rxns)
    print(f'all_reactions: {all_reactions}')
    print(f'reaction_pathways: {reaction_pathways}')
    return all_reactions, reaction_pathways


def pathways_with_conditions(dir, file_name, all_reactions, reaction_pathways):
    """
    Helper function that takes in a list of all reactions (all_reactions) and a dict with
    desired product mapped to a tuple of SMILES strings representing reactions (reaction_pathways).
    Generates conditions and top products for each set of conditions for each reaction.
    Then filters to keep sets of conditions that generate the desired product.
    Finally, loops through all of the tuples in reaction_pathways and grabs a list of dictionaries
    of conditions for each of the reactions in the tuples.
    Returns a dictionary containing this information.
    """
    rxn_contexts_and_preds = ForwardPred.contexts_and_preds(all_reactions)
    print('finished predicting contexts and top products')
    # save file
    file_path1 = file_name + "_contexts_and_preds.json"
    file_path1 = os.path.join(dir, file_path1)
    save_file(file_path1, rxn_contexts_and_preds)
    print(f'saved contexts and top products to {file_path1}')

    valid_conditions = ForwardPred.compare_products(rxn_contexts_and_preds)
    print('finished comparing top products to desired product')
    # save file
    file_path2 = file_name + "_matching_prods.json"
    file_path2 = os.path.join(dir, file_path2)
    save_file(file_path2, valid_conditions)
    print(f'saved conditions with matching products to {file_path2}')

    pathways_with_conditions = {}
    for product, pathway in reaction_pathways.items():
        pathway_with_conditions = tuple()
        for reaction in pathway:
            reaction_dict = {}
            reaction_dict[reaction] = valid_conditions[reaction]
            pathway_with_conditions += (reaction_dict,)
        pathways_with_conditions[product] = pathway_with_conditions
    print('finished adding conditions to pathways')
    pathways_path = file_name + "_pathways_w_conditions.json"
    save_file(pathways_path, pathways_with_conditions)
    print(f"File saved to {pathways_path}")

    return pathways_with_conditions


if __name__ == "__main__":
    # test case for react()
    # reactant1 = "C1=CC(=C(N=C1)N)Br"
    # reactant2 = "O=C(O)c1ccc(B(O)O)cc1"
    # template_id = "ChanLam"
    # print(react(reactant1, reactant2, template_id))

    # test case for synthesis_pathway()
    # product = "O=S(=O)(N(c1ccccc1)c1ncccn1)C(F)(F)F"
    # data = {"QED": ["0.8735248151923747", 0],
    #     "similarity": ["0.579330638082168", 0],
    #     "rank": 3,
    #     "initial_mol": "B(C1=CC=CC=C1)(O)O",
    #     "templates": [{"_id": "ChanLam", "reactants": ["Nc1ncccn1"]},
    #                   {"_id": "Sulfonamide", "reactants": ["O=S(=O)(Cl)C(F)(F)F"]}]}
    # print(synthesis_pathway(product, data))

    # test case template for testing two reactants and a template ID
    # reactant1 = Chem.MolFromSmiles('CS(=O)(=O)NS(=O)(=O)c1ccc(Br)cc1')
    # reactant2 = Chem.MolFromSmiles('NS(=O)(=O)O')
    # reaction_comps = rxn_templates['Mitsunobu_amine'].split(">>")
    # reaction_smarts = reaction_comps[1] + ">>" + reaction_comps[0]
    # rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    # products = rxn.RunReactants((reactant1, reactant2))
    # print(products)
    # product = products[0][0]
    # product_smiles = Chem.MolToSmiles(product)
    # print(product_smiles)

    # load json file with products and their associated data
    products_path = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/MFBO_selected_mols_for_synthesis.json'
    products_data = read_file(products_path)
    all_reactions, reaction_pathways = make_synthesis_pathways(products_data)
    dir = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols'
    file_name = 'MFBO_selected_mols'
    pathways_with_conditions(dir, file_name, all_reactions, reaction_pathways)
