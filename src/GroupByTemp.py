# imports
import json
import os
import numpy as np


# /////////////////////////////
# FILTERING SYNTHESIS PATHWAYS
# /////////////////////////////

def save_file(file_path, data):
    """
    Helper function that saves data into a json file.
    """
    with open(file_path, 'w') as outfile:
        json.dump(data, outfile, indent=4)
    print(f"File saved to {file_path}")


def check_empty_or_small_prob(rxn_dict):
    """
    Helper function that checks if the list of conditions for a reaction is empty.

    Input:
    - A dictionary with a single key-value pair, where:
    - Key: Reaction SMILES (string)
    - Value: A list of dictionaries, each representing a set of conditions that produce the desired product.

    Returns:
    - A Boolean indicating whether the list of conditions is empty.
    """
    for k,rxns_list in rxn_dict.items():
        if rxns_list == []:
            return True
        probs = [condition["prob"] for condition in rxns_list]
        if max(probs) < 0.05:
            return True
        else:
            return False


def filter_data(pathways_with_conditions):
    """
    Helper function that filters out products with invalid reaction steps.

    Input:
    - A dictionary where:
    - Key: Product SMILES (string)
    - Value: A list of reactions, each represented as a dictionary with a single key-value pair:
        - Key: Reaction SMILES (string)
        - Value: A list of dictionaries, each specifying a set of conditions that produce the desired product.

    Returns:
    - A new dictionary in the same format, excluding products whose synthesis pathways contain reactions without valid conditions.
    """
    filtered_pathways = {}
    for product, rxns_list in pathways_with_conditions.items():
        valid = True
        for rxn in rxns_list:
            if check_empty_or_small_prob(rxn):
                valid = False # stop checking remaining rxns in pathway
                break
        # continue to next product if one the rxns doesn't have valid conditions
        if not valid:
            continue
        else:
            filtered_pathways[product] = rxns_list
    return filtered_pathways


# ////////////////////////////////////////////////////
# TRANSFORMING DATA TO STREAMLINE TEMPERATURE BINNING
# ///////////////////////////////////////////////////


def increment_count(rxns_count, reaction_smiles):
    """
    Helper function that directly mutates rxns_count to update the
    count associated with a particular reaction.
    """
    if reaction_smiles in rxns_count.keys():
        rxns_count[reaction_smiles] += 1
    else:
        rxns_count[reaction_smiles] = 0
    count = rxns_count[reaction_smiles]
    return count


def transform_data(filtered_pathways):
    """
    Helper function that takes in the final output from SynthesisPathway and
    returns a tuple of 3 elements.
    1st element: Dictionary where each key is a reaction SMILES and each value
    is a dict storing the set of conditions with the highest probability of success.
    2nd element: Set of all the reaction strings.
    3rd element: Dictionary where each key is a reaction SMILES and each value is
    a list of lists containing the synthesis pathways it's in.
    """
    # Initialize variables
    rxns_count = {}
    rxns_top_conditions = {}
    all_rxns = set()
    rxns_to_pathways = {}

    # pathway is a list
    for product, pathway in filtered_pathways.items():
        current_pathway = []
        # each pathway is a list where elements are dicts
        for reaction in pathway:
            # conditions is a list of dictionaries
            for reaction_smiles, conditions in reaction.items():
                max_prob = 0
                top_condition = None
                for condition in conditions:
                    if condition["prob"] > max_prob:
                        max_prob = condition["prob"]
                        top_condition = condition
                rxns_top_conditions[reaction_smiles] = top_condition
                # modify reaction smiles
                count = increment_count(rxns_count, reaction_smiles)
                # tag each reaction with an ID to distinguish same rxns in
                # different pathways from each other
                modified_reaction_smiles = f'{str(count)}_{reaction_smiles}'
                all_rxns.add(modified_reaction_smiles)
                current_pathway.append(modified_reaction_smiles)
        # reaction is already tagged
        for reaction in current_pathway:
            rxns_to_pathways[reaction] = current_pathway
    return rxns_top_conditions, all_rxns, rxns_to_pathways


# //////////////////////////////
# ASSIGNING REACTIONS TO PLATES
# //////////////////////////////


def make_next_group(completed, uncompleted, rxns_to_pathways):
    """
    Helper function that takes in a set of completed reactions, a set of uncompleted reactions,
    and a dictionary of reactions mapped to the synthesis pathway they're part of.

    Returns the next group of reactions that can be completed.
    """
    next_group = []
    for current_reaction in uncompleted:
        pathway = rxns_to_pathways[current_reaction]
        idx = pathway.index(current_reaction)

        if idx == 0 or pathway[idx - 1] in completed:
            next_group.append(current_reaction)

    return next_group


def sort_rxns_by_temp(next_group, rxns_top_conditions):
    """
    Helper function that takes in the next group of possible reactions and returns
    a list of tuples, in which the reactions are sorted in order of temperature.
    """
    reactions_list = []

    # Each reaction is tagged with a unique ID
    for reaction in next_group:
        rxn_without_id = reaction[2:]
        rxn_and_temp_tuple = (reaction, rxns_top_conditions[rxn_without_id]["temperature"])
        reactions_list.append(rxn_and_temp_tuple)

    # Sort the list of tuples by temperature
    sorted_reactions = sorted(reactions_list, key=lambda x: x[1])
    return sorted_reactions


def log_bins(sorted_next_group, num_bins=8):
    """
    Helper function that processes a list of unsorted, uncompleted reactions.

    Input:
    - `next_group`: A list of reaction SMILES with associated temperatures.

    Returns:
    - An array of bin edges representing temperature ranges, calculated on a logarithmic scale.
    """
    min_temp = min(sorted_next_group, key=lambda x: x[1])[1]
    max_temp = max(sorted_next_group, key=lambda x: x[1])[1]

    # +273 converts into Kelvin and -273 at end converts back into Celsius; bins get bigger
    # as temp increases since this is a logarithmic scale
    bin_edges = np.logspace(np.log10(min_temp+273), np.log10(max_temp+273), num_bins+1) - 273
    return bin_edges


def temperature_groups(next_group, rxns_top_conditions):
    """
    Creates logarithmic temperature bins and groups reactions
    based on their temperatures.

    Parameters:
    - next_group: List of reaction SMILES.
    - rxns_top_conditions: Dictionary mapping reaction SMILES to their top conditions, including temperature.

    Returns:
    - A dictionary where keys are temperature ranges (as strings) and values are sets of reactions within those ranges.
    """
    # Initialize set of reactions and sort them by temperature
    rxns_set = set(next_group)
    sorted_next_group = sort_rxns_by_temp(next_group, rxns_top_conditions)

    # Generate logarithmic bin edges
    bin_edges = log_bins(sorted_next_group)
    temp_groups = {}

    # Initialize temperature groups in the dict
    for i in range(len(bin_edges)-1):
        min_temp = int(float(bin_edges[i]))
        max_temp = int(float(bin_edges[i+1]))
        temp_groups[f'{min_temp}_{max_temp}'] = set()

    # Assign reactions to appropriate temperature groups
    for rxn, temp in sorted_next_group:
        for temp_range, rxns in temp_groups.items():
            min_temp = int(temp_range.split('_')[0])
            max_temp = int(temp_range.split('_')[1])
            if temp < max_temp and temp >= min_temp:
                rxns.add(rxn)
                rxns_set.remove(rxn)

    # Handle any remaining reactions not assigned to a bin
    if rxns_set != set():
        for rxn in rxns_set:
            rxn_no_id = rxn[2:]
            min_temp = int(rxns_top_conditions[rxn_no_id]["temperature"])
            max_temp = min_temp + 10
            temp_groups[f'{min_temp}_{max_temp}'] = {rxn}

    # Return temperature groups (may contain empty sets)
    return temp_groups


def largest_temp_groups(temp_groups, num=5):
    """
    Helper function that returns up to the top [num] (default 5)
    largest temperature groups.
    temp_groups is a dict where each key is a temperature range and each
    value is a set of reactions in that range.
    Returns a list of dicts.
    """
    """
    Helper function that returns up to the top `num` (default 5)
    largest temperature groups.

    Parameters:
    - `temp_groups`: Dictionary where keys are temperature ranges (strings) and values are sets of reactions within those ranges.
    - `num`: Maximum number of top temperature groups to return (default is 5).

    Returns:
    - A list of dictionaries representing the top temperature groups.

    """
    sorted_list = sorted(temp_groups.items(), key=lambda item: len(item[1]), reverse=True)

    # Select the top five dictionaries with the largest lists
    top_five = sorted_list[:5]

    return top_five


def update_step(recent_group, completed, uncompleted):
    """
    Helper function that processes a set of reaction SMILES.

    Parameters:
    - `to_be_performed`: Set of reaction SMILES to be performed.
    - `completed`: Set of completed reactions.
    - `uncompleted`: Set of uncompleted reactions.

    Returns:
    - A 2-element tuple containing updated sets: (new_completed, new_uncompleted).
    """
    new_completed = recent_group | completed
    new_uncompleted = uncompleted - recent_group
    return new_completed, new_uncompleted


def get_conditions(recent_group, rxns_top_conditions):
    """
    Helper function that processes a list of reaction SMILES and their top conditions.

    Parameters:
    - `reactions`: List of reaction SMILES.
    - `top_conditions`: Dictionary mapping reaction SMILES to their top reaction conditions.

    Returns:
    - A list of dictionaries, where each dictionary contains a single key-value pair:
    - Key: Reaction SMILES (string).
    - Value: Dictionary of optimal reaction conditions.
    """
    recent_group_with_conditions = []
    for reaction in recent_group:
        rxn_with_conditions = {}
        rxn_without_id = reaction[2:]
        rxn_with_conditions[reaction] = rxns_top_conditions[rxn_without_id]
        recent_group_with_conditions.append(rxn_with_conditions)
    return recent_group_with_conditions


def top_sequences(filtered_pathways, num):
    """
    Generates a list of possible wellplate sequences by branching on the top temperature groups.

    This function recursively explores different wellplate sequences by selecting the largest possible
    temperature groups at each iteration. It starts with the top `num` largest temperature groups in
    the first iteration, then branches into further groups in subsequent iterations.

    Parameters:
    - filtered_pathways: Dictionary where keys are product SMILES strings and values are lists of reaction pathways.
    - num: Integer specifying the number of top temperature groups to branch on at each iteration.

    Returns:
    - A list of dictionaries, where each dictionary represents a wellplate sequence. Each key in the
      dictionary is a wellplate ID (combining plate index and temperature range), and each value is a list
      of reactions with their corresponding conditions.
    """
    def generate_sequences(completed, uncompleted, rxns_top_conditions, rxns_to_pathways, num, depth=3):
        """
        Recursive helper function to generate wellplate sequences by branching on top temperature groups.
        """
        if depth == 0 or not uncompleted:
            return [{}]  # Base case: return an empty wellplate sequence

        next_group = make_next_group(completed, uncompleted, rxns_to_pathways)
        temp_groups = temperature_groups(next_group, rxns_top_conditions)
        top_temp_groups = largest_temp_groups(temp_groups, num)

        sequences = []
        for temp_range, recent_group in top_temp_groups:
            new_completed, new_uncompleted = update_step(recent_group, completed, uncompleted)
            recent_group_with_conditions = get_conditions(recent_group, rxns_top_conditions)
            sub_sequences = generate_sequences(new_completed, new_uncompleted, rxns_top_conditions, rxns_to_pathways, num, depth - 1)

            for sub_seq in sub_sequences:
                plate_id = len(sub_seq)
                wellplate_key = f'{plate_id}_{temp_range}'
                new_wellplate = sub_seq.copy()
                new_wellplate[wellplate_key] = recent_group_with_conditions
                sequences.append(new_wellplate)

        return sequences

    # Initialize variables
    rxns_top_conditions, uncompleted, rxns_to_pathways = transform_data(filtered_pathways)
    completed = set()

    # Generate sequences starting from the top temperature groups
    top_seqs = generate_sequences(completed, uncompleted, rxns_top_conditions, rxns_to_pathways, num)
    return top_seqs


def save_completed_paths(completed_paths, dir, filename):
    """
    Saves each completed path to a separate file.
    """
    id = 0
    for path in completed_paths:
        no_dir_filepath = f'{filename}_{id}'
        filepath = os.path.join(dir, no_dir_filepath)
        with open(filepath, 'w') as outfile:
            json.dump(path, outfile, indent=4)
        print(f'File saved to {filepath}')
        id += 1


# /////////////////////////////////////////
# EVALUATING GENERATED WELLPLATE SEQUENCES
# /////////////////////////////////////////

def iterate_reactions(rxns, products):
    """
    Helper function that counts the number of products
    in a given list of reactions.
    """
    update_count = 0
    # eaxh rxn is a dict with one k,v pair
    for rxn in rxns:
        for k in rxn.keys():
            if k.split('>>')[1] in products:
                update_count += 1
    return update_count


def count_products(wellplates, num, products):
    """
    Helper function that returns the number of products produced in the
    first "num" wellplates of a sequence.
    """
    count = 0
    for i in range(0, num):
        for k in wellplates.keys():
            if k.startswith(f'{str(i)}_'):
                rxns = wellplates[k]
                update_count = iterate_reactions(rxns, products)
                count += update_count
                break
    return count


def max_count_sequence(dir, filename, sequences, num, filtered_pathways):
    """
    Returns the ID of the wellplate sequence with the greatest number
    of products produced in the first `num` wellplates.
    """
    products = set(filtered_pathways.keys())
    product_counts = []

    for id in range(0, sequences+1):
        filepath = os.path.join(dir, filename)
        filepath += f'_{id}'
        with open(filepath, 'r') as jsonfile:
            wellplates = json.load(jsonfile)
        count = count_products(wellplates, num, products)
        product_counts.append(count)

    id_of_greatest = product_counts.index(max(product_counts)) + 1
    max_count = max(product_counts)
    return id_of_greatest, max_count, product_counts


if __name__ == "__main__":
    """
    Step 1: Filter candidate molecules
    -----------------------------------
    First, filter out candidate molecules whose synthesis pathways contain intermediate reactions
    with low probability of success or no valid conditions.
    """
    # Uncomment to generate filtered pathways
    # pathways_with_conditions_filepath = 'path/to/iter2_pathways_w_conditions.json'
    # with open(pathways_with_conditions_filepath, 'r') as jsonfile:
    #     pathways_with_conditions = json.load(jsonfile)

    # filtered_pathways = filter_data(pathways_with_conditions)
    filtered_pathways_filepath = 'path/to/iter2_filtered_pathways.json'
    # with open(filtered_pathways_filepath, 'w') as outfile:
    #     json.dump(filtered_pathways, outfile, indent=4)

    """
    Step 2: Generate wellplate sequences
    -------------------------------------
    Load filtered pathways and determine optimal wellplate sequences by branching on top temperature groups.
    """
    with open(filtered_pathways_filepath, 'r') as jsonfile:
        filtered_pathways = json.load(jsonfile)

    rxns_top_conditions, uncompleted, rxns_to_pathways = transform_data(filtered_pathways)

    # Uncomment to inspect intermediate steps
    # print(rxns_to_pathways)
    # completed = set()
    # next_group = make_next_group(completed, uncompleted, rxns_to_pathways)
    # print(sort_rxns_by_temp(next_group, rxns_top_conditions))
    # temp_groups = temperature_groups(next_group, rxns_top_conditions)
    # print(temp_groups)
    # next_temp_groups = largest_temp_groups(temp_groups, 5)
    # print(len(next_temp_groups))

    top_seqs = top_sequences(filtered_pathways, 5)

    # Save tagged reactions to pathways
    out_filepath = 'path/to/iter2_tagged_rxns_to_pathways.json'
    with open(out_filepath, 'w') as outfile:
        json.dump(rxns_to_pathways, outfile, indent=4)

    """
    Step 3: Evaluate wellplate sequences
    -------------------------------------
    Count the number of products covered by each wellplate sequence and select the best sequence.
    """
    counts = []
    max_count = 0
    best_seq = None

    for i, top_seq in enumerate(top_seqs):
        products = set(filtered_pathways.keys())

        # Check different numbers of wellplates here, as is reasonable within
        # the constraints of the autonomous liquid handler
        count = count_products(top_seq, 4, products)
        if count > max_count:
            max_count = count
            best_seq = top_seq
        counts.append(count)
        print(f'{i}: {count}')

    print(f'Maximum count: {max(counts)}')
    print(f'Best sequence count: {count_products(best_seq, 4, products)}')

    # Uncomment to save the best wellplate sequence
    best_seq_filepath = 'path/to/iter2_best_wellplate_seq.json'
    # with open(best_seq_filepath, 'w') as outfile:
    #     json.dump(best_seq, outfile, indent=4)
