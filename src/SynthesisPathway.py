from rdkit import Chem
from rdkit.Chem import AllChem
import json
import ForwardPred
import os

# Templates to run reactions
templates_dict_path = './templates_dict.json'

with open(templates_dict_path, 'r') as file:
    rxn_templates = json.load(file)


def react(reactants, template_id):
    """
    Helper function that returns a SMILES of the product,
    given a tuple of reactant SMILES and a reaction template ID.
    """
    reactants = [Chem.MolFromSmiles(reactant) for reactant in reactants]
    reactants = tuple(reactants)
    reaction_comps = rxn_templates[template_id].split(">>")
    reaction_smarts = reaction_comps[1] + ">>" + reaction_comps[0]
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    products = rxn.RunReactants(reactants)
    product_smiles = ""
    if products:
        product = products[0][0]  # Assuming one product
        product_smiles = Chem.MolToSmiles(product)
    else:
        reversed_reactants = reactants[::-1]
        products = rxn.RunReactants(reversed_reactants)
        try:
            product = products[0][0]
            product_smiles = Chem.MolToSmiles(product)
        # if neither reactant order works
        except IndexError:
            return None
    return product_smiles


def matching_product(products, desired_product):
    """
    Helper function that loops through products and returns an RDKit
    object of the product that matches the desired product.
    """
    desired_product = Chem.MolToSmiles(Chem.MolFromSmiles(desired_product))
    product = None
    for prod in products:
        if Chem.MolToSmiles(prod[0]) == desired_product:
            product = prod[0]
    if product is not None:
        return product
    else:
        return None


def react(reactants, template_id, check_prods, desired_product=None):
    """
    Helper function that returns a SMILES of the product,
    given a tuple of reactant SMILES and a reaction template ID.
    """
    reactants = [Chem.MolFromSmiles(reactant) for reactant in reactants]
    reactants = tuple(reactants)
    reaction_comps = rxn_templates[template_id].split(">>")
    reaction_smarts = reaction_comps[1] + ">>" + reaction_comps[0]
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    products = rxn.RunReactants(reactants)
    product_smiles = ""

    # Products are found in the first ordering of the reactants
    if products:
        if not check_prods:
            product = products[0][0]  # Assuming one product
        else:
            product = matching_product(products, desired_product)

    # Otherwise reverse the products and try running reaction again
    else:
        reversed_reactants = reactants[::-1]
        products = rxn.RunReactants(reversed_reactants)
        try:
            if not check_prods:
                product = products[0][0]
            else:
                product = matching_product(products, desired_product)

        # If neither reactant order works
        except IndexError:
            return None
    if product is None:
        return None
    product_smiles = Chem.MolToSmiles(product)
    return product_smiles


def make_rxn_smiles(rxn_molecules):
    """
    Helper function that takes in a 2-element tuple where the first elt
    is a tuple of reactant SMILES and the second tuple is the product.
    Returns a reaction SMILES.
    """
    reactants_tup = rxn_molecules[0]
    product = rxn_molecules[1]
    rxn_smiles = ""
    for idx, reactant in enumerate(reactants_tup):
        if idx < len(reactants_tup) - 1:
            rxn_smiles += reactant
            rxn_smiles += "."
        else:
            rxn_smiles += reactant
    tail = ">>" + product
    rxn_smiles += tail
    return rxn_smiles


def synthesis_pathway(product, data):
    """
    Helper function that given a SMILES string of the product and a dictionary
    of associated data, returns a tuple of reactions needed to synthesize the product.
    Reactions are represented as SMILES strings.
    """
    templates = data["templates"]
    reactant1 = data["initial_mol"]
    rxn_pathway = ()

    # Loop through all the reaction templates, where each template is a dict
    for idx, temp in enumerate(templates):
        template_id = temp["_id"]
        rxn_molecules = tuple()
        check_prods = False
        if idx == len(templates) - 1:
            check_prods = True

        # Assume there is only one reactant
        try:
            # If there is no reactant, move to the except block
            reactant2 = temp["reactants"][0]
            if not check_prods:
                intermediate = react((reactant1, reactant2), template_id, check_prods)
            else:
                intermediate = react((reactant1, reactant2), template_id, check_prods, product)
            rxn_molecules = ((reactant1, reactant2), intermediate)

        # The intermediate from the previous step is the only reactant
        except IndexError:
            if not check_prods:
                intermediate = react((reactant1,), template_id, check_prods)
            else:
                intermediate = react((reactant1,), template_id, check_prods, product)
            rxn_molecules = ((reactant1,), intermediate)

        # If one of the reactions in the sequence of reactions doesn't work, return immediately
        if intermediate is None:
            return None
        rxn_smiles = make_rxn_smiles(rxn_molecules)
        rxn_pathway += (rxn_smiles,)

        # Update reactant 1 to react the intermediate by the previous step
        # with the reactant in the next step
        reactant1 = intermediate

    # Check that the final product is the desired product
    if Chem.MolToSmiles(Chem.MolFromSmiles(product)) != reactant1 and product != reactant1:
        raise ValueError("The final product does not match the desired product.")
    return rxn_pathway


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
    Helper function that returns a 2-elt tuple, where the
    1st elt is a list of all the reactions and the
    2nd elt is a dictionary mapping each product SMILES to a tuple
    of ordered reactions SMILES necessary to synthesize the product.
    """
    # Store synthesis pathways in a dict
    reaction_pathways = {}

    # Simultaneously store a list of all the reactions
    all_reactions = []

    for prod, data in products_data.items():
        # rxn_pathway is a tuple of reaction strings
        pathway = synthesis_pathway(prod, data)

        # If the returned pathway is None, one of the reactions didn't
        # produce a valid result, so skip the product entirely
        if pathway is None:
            continue
        reaction_pathways[prod] = pathway
        rxns = list(reaction_pathways[prod])
        all_reactions.extend(rxns)
    return all_reactions, reaction_pathways


def pathways_with_conditions(dir, file_name, all_reactions, reaction_pathways):
    """
    Helper function that takes a list of reactions (`all_reactions`) and a dictionary (`reaction_pathways`)
    mapping desired products to tuples of reaction SMILES.

    Generates conditions and top products for each reaction, filters conditions to retain those that
    produce the desired product, and collects valid conditions for each reaction in `reaction_pathways`.

    Returns a dictionary mapping each reaction to its valid conditions.
    """
    rxn_contexts_and_preds = ForwardPred.contexts_and_preds(all_reactions)
    print('finished predicting contexts and top products')

    # Save file
    file_path1 = file_name + "_contexts_and_preds.json"
    file_path1 = os.path.join(dir, file_path1)
    save_file(file_path1, rxn_contexts_and_preds)
    print(f'saved contexts and top products to {file_path1}')

    valid_conditions = ForwardPred.compare_products(rxn_contexts_and_preds)
    print('finished comparing top products to desired product')

    # Save file
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

    # Load products data from JSON file
    products_path = './gen_mols/gen_mols_Jun2024.json'
    products_data = read_file(products_path)

    # Generate synthesis pathways and reactions
    all_reactions, reaction_pathways = make_synthesis_pathways(products_data)

    # Save synthesis pathways to a JSON file
    synthesis_pathways_path = './gen_mols/synthesis_pathways.json'
    save_file(synthesis_pathways_path, reaction_pathways)

    # Save reactions to request to a separate JSON file
    reactions_path = './gen_mols/rxns_to_request.json'
    save_file(reactions_path, all_reactions)

    # Define directory and file name for pathways with conditions
    output_dir = './MFBO_selected_mols'
    output_file_name = 'MFBO_selected_mols'

    # Generate pathways with conditions and save to files
    pathways_with_conditions(output_dir, output_file_name, all_reactions, reaction_pathways)

    print("Process completed successfully.")
