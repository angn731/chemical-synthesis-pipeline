from src.api_client import APIClient
import glob
import json
import yaml
import os
import csv

# sets ip API client to communicate with server
hostname = r'https://askcos.mit.edu:7000/'
client = APIClient(hostname, verify=False)

# ////////
# CONTEXT
# ////////

def read_file(original_filename, directory):
    """
    Helper function that reads in a json file storing a list of reaction strings
    and returns the list.
    """
    target_directory = os.path.join(directory, original_filename)
    with open(target_directory, 'r') as jsonfile:
        incoming_rxns = json.load(jsonfile)
        reactions_to_request = incoming_rxns
    return reactions_to_request


def save_file(file_path, data):
    """
    Helper function that saves data into a json file.
    """
    with open(file_path, 'w') as outfile:
        json.dump(data, outfile, indent=4)
    print(f"File saved to {file_path}")


def forward_request(reactants, rxn_conditions, size=4):
    """
    Helper function that takes a list of reactants and a list of possible
    reaction conditions (`rxn_conditions`), where each element in `rxn_conditions`
    is a dictionary representing a specific set of conditions.

    Returns a list in which each dictionary from `rxn_conditions` is updated to
    include a new key-value pair: the top 5 products predicted under the corresponding
    set of conditions.
    """
    keys = ["prob", "smiles"]
    request_groups = [rxn_conditions[i:i + size]
                    for i in range(0, len(rxn_conditions), size)]

    # Initialize the index for the element in rxn_conditions
    current_idx = 0

    # Loop through dictionaries
    for request_group in request_groups:
        print(f'forward pred current id {current_idx}')
        current_ids = []
        for conditions in request_group:

            # Make a call to the forward predictor using reactants, reagents, and solvent
            params = {'reactants':reactants,
                            'reagents': conditions["reagent"],
                            'solvent': conditions["solvent"],
                            'num_results': 5}
            result = client.post('forward', data = params)
            current_ids.append(result['task_id'])

            # Attempt to retrieve the results for each task
            try:
                task_results = [client.get_result(tid, timeout=600, interval=5) for tid in current_ids]
                # print('retrieving forward preds')
            except KeyError as e:
                print(f'KeyError: {e}')
                print(f'params: {params}')
                continue
            except json.JSONDecodeError:
                print('\tJSON DECODE ERROR!')
                continue

            # Parse task_results, a list in which each element outputs the top 5 products
            for idx, task_result in enumerate(task_results):
                top_products = []
                for product in task_result["output"]:
                    top_products.append({key: product[key] for key in keys})
                rxn_conditions[current_idx + idx]["top_products"] = top_products

        # Before moving on to the next request_group, update the index pointer in rxn_conditions
        current_idx += size

    return rxn_conditions


def contexts_and_preds(reactions_to_request):
    """
    Helper function that takes a list of SMILES strings representing reactions and returns a dictionary.

    Each key is a SMILES string of a reaction, and the corresponding value is a list of dictionaries,
    where each dictionary contains the top reaction conditions (from the ASCKOS context endpoint)
    and the associated top products (from the forward predictor ASCKOS endpoint).
    """
    # Dictionary to store context data for each reaction
    contexts = {}
    current_contexts = {}
    conditions = ["temperature", "solvent", "reagent", "catalyst"]
    size = 4
    request_groups = [reactions_to_request[i:i + size]
                    for i in range(0, len(reactions_to_request), size)]
    mol_id = 0
    req_group_id = 0
    for request_group in request_groups:
        req_group_id += 1
        for rxn_smiles in request_group:
            mol_id += 1
            print(f'on reaction {mol_id}')
            current_ids = []
            params = {'reactants':rxn_smiles.split('>>')[0],
                            'products': rxn_smiles.split('>>')[1],
                            'return_scores':True,
                            'num_results': 5}
            result = client.post('context', data = params)

            # Host server assigns each task an ID, which is stored in the current_ids list,
            # where each list element is a tuple
            current_ids.append((rxn_smiles, result['task_id']))

            # attempt to retrieve task results using task ID
            try:
                task_results = [client.get_result(tid[1], timeout=600, interval=5) for tid in current_ids]
                # print('retrieving contexts')
            except KeyError as e:
                print(f'KeyError: {e}')
                print(f'params: {params}')
                continue
            except json.JSONDecodeError:
                print('\tJSON DECODE ERROR!')
                continue

            # Update contexts dictionary with task results
            for current_index, task_result in enumerate(task_results):
                # parse task_result, which is a dictionary
                rxn_result = []
                # Each possibility is a dictionary
                for possibility in task_result["output"]:
                    rxn_result.append({key: possibility[key] for key in conditions})
                # Each element in rxn_result is a dictionary with one set of possible conditions
                smiles_str = current_ids[current_index][0]
                rxn_result = forward_request(smiles_str.split('>>')[0], rxn_result)
                # current_ids[current_index][0] is a SMILES string
                contexts[smiles_str]= rxn_result
                current_contexts[smiles_str]= rxn_result

        # Save data to file for every three groups
        if req_group_id % 3 == 0:
            file_path = f'./data/gen_mols_contexts_and_preds/group_{req_group_id // 3}.json'
            save_file(file_path, current_contexts)
            current_contexts = {}
    return contexts


def all_contexts_and_preds(folder_path, output_path):
    """
    In case of an error when running contexts_and_preds(), combine data from all the intermediate
    contexts_and_preds files into one file.
    """
    combined_data = {}

    # Iterate through all files in the folder
    for filename in os.listdir(folder_path):
        if filename.startswith('group'):
            file_path = os.path.join(folder_path, filename)
            with open(file_path, 'r') as file:
                data = json.load(file)
                # Assume each JSON file contains a single dictionary
                combined_data.update(data)

    # Write combined data to a new JSON file
    with open(output_path, 'w') as file:
        json.dump(combined_data, file, indent=4)

    print(f"Combined data has been written to {output_path}")


# ////////////////////////////////////////
# COMPARE TOP PRODUCTS TO TARGET PRODUCT
# ////////////////////////////////////////

def compare_products(contexts):
    """
    Helper function that returns a dictionary mapping each reaction SMILES to a list of
    predicted reaction conditions that generate the desired product.

    Each list element is a dictionary containing the reaction conditions and the associated
    probability of producing the desired product, as predicted by the forward prediction endpoint.
    """

    # Initialize dict
    matching_product = {}
    condition_params = ["temperature", "solvent", "reagent", "catalyst"]
    for rxn_smiles in contexts:
        print('next reaction')
        rxn_product = rxn_smiles.split('>>')[1]
        valid_conditions = []
        # Each rxn_condition is a dict
        for rxn_condition in contexts[rxn_smiles]:
            top_products = rxn_condition["top_products"]
            # Each prod is a dict
            for prod in top_products:
                if prod["smiles"] == rxn_product:
                    valid = {key: rxn_condition[key] for key in condition_params}
                    valid["prob"] = prod["prob"]
                    valid_conditions.append(valid)
        matching_product[rxn_smiles] = valid_conditions
    return matching_product


def match_pathways_to_conditions(reaction_pathways, valid_conditions):
    """
    Helper function that maps each reaction in a pathway to a list of top conditions
    for producing the specified product.

    Parameters:
    - `reaction_pathways`: A dictionary mapping product SMILES to tuples representing their synthesis pathways.
    - `valid_conditions`: A dictionary mapping reaction SMILES to the top conditions for producing the specified product.
    """

    pathways_with_conditions = {}
    for product, pathway in reaction_pathways.items():
        pathway_with_conditions = tuple()
        for reaction in pathway:
            reaction_dict = {}
            reaction_dict[reaction] = valid_conditions[reaction]
            pathway_with_conditions += (reaction_dict,)
        pathways_with_conditions[product] = pathway_with_conditions
    return pathways_with_conditions


def pathways_with_conditions(dir, file_name, all_reactions, reaction_pathways):
    """
    Takes a list of reactions (`all_reactions`) and a dictionary (`reaction_pathways`)
    mapping desired products to tuples of reaction SMILES.

    Generates conditions and top products for each set of conditions, filters conditions to retain those
    producing the desired product, and collects valid conditions for each reaction in `reaction_pathways`.

    Returns a dictionary of filtered conditions for all reactions.
    """
    rxn_contexts_and_preds = contexts_and_preds(all_reactions)
    print('finished predicting contexts and top products')
    # save file
    file_path1 = file_name + "_all_contexts_and_preds.json"
    file_path1 = os.path.join(dir, file_path1)
    save_file(file_path1, rxn_contexts_and_preds)
    print(f'saved contexts and top products to {file_path1}')

    valid_conditions = compare_products(rxn_contexts_and_preds)
    print('finished comparing top products to desired product')
    # save file
    file_path2 = file_name + "_matching_prods.json"
    file_path2 = os.path.join(dir, file_path2)
    save_file(file_path2, valid_conditions)
    print(f'saved conditions with matching products to {file_path2}')

    pathways_with_conditions = match_pathways_to_conditions(reaction_pathways, valid_conditions)
    print('finished adding conditions to pathways')
    pathways_path = file_name + "_pathways_w_conditions.json"
    save_file(pathways_path, pathways_with_conditions)
    print(f"File saved to {pathways_path}")

    return pathways_with_conditions


if __name__ == "__main__":

    # Step 1: Load the list of reaction SMILES strings
    original_filename = 'rxns_to_request.json'
    directory = "./gen_mols"
    reactions_to_request = read_file(original_filename, directory)

    # Step 2: Load synthesis pathways
    pathways_filepath = os.path.join(directory, 'synthesis_pathways.json')
    with open(pathways_filepath, 'r') as file:
        reaction_pathways = json.load(file)

    # Step 3: Load all contexts and predictions
    folder_path = os.path.join(directory, 'gen_mols_contexts_and_preds')
    output_path = os.path.join(folder_path, 'gen_mols_all_contexts_and_preds.json')
    with open(output_path, 'r') as file:
        all_contexts_and_predictions = json.load(file)

    # Step 4: Filter reactions for iteration 2 pathways
    iter2_pathways_filepath = './mols_iter2/iter2_synthesis_pathways.json'
    with open(iter2_pathways_filepath, 'r') as file:
        iter2_pathways = json.load(file)

    iter2_reactions = set()
    for product, pathway in iter2_pathways.items():
        iter2_reactions.update(pathway)

    iter2_contexts_and_preds = {rxn: data for rxn, data in all_contexts_and_predictions.items() if rxn in iter2_reactions}

    # Step 5: Save filtered contexts and predictions for iteration 2
    iter2_contexts_and_preds_filepath = './mols_iter2/iter2_contexts_and_preds.json'
    save_file(iter2_contexts_and_preds_filepath, iter2_contexts_and_preds)

    # Step 6: Compare products and save valid conditions
    valid_conditions = compare_products(iter2_contexts_and_preds)
    matching_prods_filepath = './mols_iter2/iter2_matching_prods.json'
    save_file(matching_prods_filepath, valid_conditions)

    # Step 7: Map synthesis pathways to valid conditions and save
    pathways_w_conditions = match_pathways_to_conditions(iter2_pathways, valid_conditions)
    pathways_w_conditions_filepath = './mols_iter2/iter2_pathways_w_conditions.json'
    save_file(pathways_w_conditions_filepath, pathways_w_conditions)
