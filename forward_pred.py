from api_client import APIClient
import glob
import json
import yaml
import os
import csv

# sets ip API client to communicate with server
hostname = r'https://askcos.mit.edu:7000/'
client = APIClient(hostname, verify=False)

original_filename = 'test_smiles_askcos_forward_pred.json'  # json of smiles to try and synthesize
directory = r"/Users/angelinaning/Downloads/jensen_lab_urop"
target_directory = os.path.join(directory, original_filename)

# //////////////////
# CONTEXT
# //////////////////

# first create list of SMILE strings of rxns to search by reading in json files
# reactions_to_request = []
# files = glob.glob(os.path.join(directory, "*_rxns_*.json"))
# for file in files:
#     with open(file, 'r') as jsonfile:
#         incoming_rxns = json.load(jsonfile)
#         reactions_to_request.extend(incoming_rxns)
with open(target_directory, 'r') as jsonfile:
    incoming_rxns = json.load(jsonfile)
    reactions_to_request = incoming_rxns
    print(f'reactions_to_request: {reactions_to_request}')

# helper function for forward request groups
def forward_request(reactants, rxn_conditions, size=4):
    """
    rxn_conditions is a list where each element is a dictionary with one set of possible conditions.
    Returns a list where each dictionary element is updated with a key-value pair specifying
    the top 5 products, given the specified set of conditions.
    """
    print('running forward pred')
    keys = ["prob", "smiles"]
    request_groups = [rxn_conditions[i:i + size]
                    for i in range(0, len(rxn_conditions), size)]
    # initialize the index for the element in rxn_conditions
    current_idx = 0
    # loop through dictionaries
    for request_group in request_groups:
        current_ids = []
        for conditions in request_group:
            # make a call to the forward predictor using reactants, reagents, and solvent
            params = {'reactants':reactants,
                            'reagents': conditions["reagent"],
                            'solvent': conditions["solvent"],
                            'num_results': 5}
            result = client.post('forward', data = params)
            current_ids.append(result['task_id'])
            # attempt to retrieve the results for each task
            try:
                task_results = [client.get_result(tid, timeout=600, interval=5) for tid in current_ids]
                print('retrieving forward preds')
            except KeyError as e:
                print(f'KeyError: {e}')
                print(f'params: {params}')
                continue
            except json.JSONDecodeError:
                print('\tJSON DECODE ERROR!')
                continue
            # parse task_results, a list in which each element outputs the top 5 products
            for idx, task_result in enumerate(task_results):
                top_products = []
                for product in task_result["output"]:
                    top_products.append({key: product[key] for key in keys})
                rxn_conditions[current_idx + idx]["top_products"] = top_products
        # before moving on to the next request_group, update the index pointer in rxn_conditions
        current_idx += size

    return rxn_conditions

# then make a request to host
# dict to store context data for each rxn
contexts = {}
conditions = ["temperature", "solvent", "reagent", "catalyst"]
size = 4
request_groups = [reactions_to_request[i:i + size]
                  for i in range(0, len(reactions_to_request), size)]
for request_group in request_groups:
    for rxn_smiles in request_group:
        current_ids = []
        params = {'reactants':rxn_smiles.split('>>')[0],
                          'products': rxn_smiles.split('>>')[1],
                          'return_scores':True,
                          'num_results': 5}
        result = client.post('context', data = params)
        # host server assigns each task an ID, which is stored in the current_ids list,
        # where each list element is a tuple
        current_ids.append((rxn_smiles, result['task_id']))
        # attempt to retrieve task results using task ID
        try:
            task_results = [client.get_result(tid[1], timeout=600, interval=5) for tid in current_ids]
            print('retrieving contexts')
        except KeyError as e:
            print(f'KeyError: {e}')
            print(f'params: {params}')
            continue
        except json.JSONDecodeError:
            print('\tJSON DECODE ERROR!')
            continue
        # update contexts dictionary with task results
        for current_index, task_result in enumerate(task_results):
            # parse task_result, which is a dictionary
            rxn_result = []
            # each possibility is a dictionary
            for possibility in task_result["output"]:
                rxn_result.append({key: possibility[key] for key in conditions})
            # each element in rxn_result is a dictionary with one set of possible conditions
            smiles_str = current_ids[current_index][0]
            rxn_result = forward_request(smiles_str.split('>>')[0], rxn_result)
            # current_ids[current_index][0] is a smiles string
            contexts[smiles_str]= rxn_result

# Combine directory and filename
file_path = os.path.join(directory, f'{original_filename.split('.')[0]}_contexts_and_pred_prods.json')

# Save the contexts dictionary to a JSON file
with open(file_path, 'w') as outfile:
    json.dump(contexts, outfile, indent=4)

print(f"Contexts and predicted products saved to {file_path}")


# ////////////////////////////////////////
# COMPARE TOP PRODUCTS TO TARGET PRODUCT
# ////////////////////////////////////////

# compare top products to target product and output a file that maps each rxn SMILES string
# to the condition(s) that generates the desired product (including probability); if none of
# top product match, output an empty list

# original_filename = 'test_smiles_askcos_forward_pred_contexts_and_pred_prods.json'  # json of smiles to try and synthesize
# target_directory = os.path.join(directory, original_filename)

with open(file_path, 'r') as jsonfile:
    contexts = json.load(jsonfile)

# initialize dict
matching_product = {}
condition_params = ["temperature", "solvent", "reagent", "catalyst"]
for rxn_smiles in contexts:
    print('next reaction')
    rxn_product = rxn_smiles.split('>> ')[1]
    # print(f'rxn_product: {rxn_product}')
    valid_conditions = []
    # each rxn_condition is a dict
    for rxn_condition in contexts[rxn_smiles]:
        # print(f'condition: {rxn_condition}')
        top_products = rxn_condition["top_products"]
        # each prod is a dict
        for idx, prod in enumerate(top_products):
            if prod["smiles"] == rxn_product:
                valid = {key: rxn_condition[key] for key in condition_params}
                valid["prob"] = prod["prob"]
                valid_conditions.append(valid)
    matching_product[rxn_smiles] = valid_conditions

# Combine directory and filename
file_path = os.path.join(directory, f'{original_filename.split('_contexts_and_pred_prods')[0]}_contexts_matching_prods.json')

# Save the contexts dictionary to a JSON file
with open(file_path, 'w') as outfile:
    json.dump(matching_product, outfile, indent=4)

print(f"Matching product contexts saved to {file_path}")


# if __name__ == "__main__":
    # test case for forward prediction
    # reactants = "CC1CNCC(C)O1.O=[N+]([O-])c1ccc(F)nc1"
    # rxn_conditions = [
    #     {
    #         "temperature": 82.41651153564453,
    #         "solvent": "CN(C)C=O",
    #         "reagents": "O=C([O-])[O-].[K+]",
    #         "catalyst": ""
    #     },
    #     {
    #         "temperature": 88.05548095703125,
    #         "solvent": "CS(C)=O",
    #         "reagents": "O=C([O-])[O-].[K+]",
    #         "catalyst": ""
    #     },
    #     {
    #         "temperature": 81.92193603515625,
    #         "solvent": "CC#N",
    #         "reagents": "O=C([O-])[O-].[K+]",
    #         "catalyst": ""
    #     },
    #     {
    #         "temperature": 87.54207611083984,
    #         "solvent": "CS(C)=O",
    #         "reagents": "",
    #         "catalyst": ""
    #     },
    #     {
    #         "temperature": 81.73251342773438,
    #         "solvent": "CN(C)C=O",
    #         "reagents": "",
    #         "catalyst": ""
    #     }
    # ]
    # print(forward_request(reactants, rxn_conditions))
