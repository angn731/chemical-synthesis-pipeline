# imports
import json
import os

# input data is a dictionary where keys are product SMILES strings
# and values are a list of reactions, each of which is a dictionary
# mapping the reaction SMILES to a list of 5 reaction conditions,
# each of which is a dictionary

# first step is to transform the data: need to map each reaction
# to the set of conditions with the highest probability
# to do so, make a dict that matches eaches reaction to dict with optimal conditions
# but i also need to keep information about order
# when i'm transforming the data, i need to create another data structure
# that stores information about the reaction pathway associated with each reaction
# (actually can pull this from synthesize_pathway())
# somewhere along the way create a set of reactions that have been completed

def save_file(file_path, data):
    """
    Helper function that saves data into a json file.
    """
    with open(file_path, 'w') as outfile:
        json.dump(data, outfile, indent=4)
    print(f"File saved to {file_path}")


def check_empty_or_small_prob(rxn_dict):
    """
    Helper function that checks if the list of conditions associated with a
    reaction is empty.
    Input is a dict containing one k,v pair, where the key is the reaction SMILES
    and the value is a list of dicts, each of which stores a set of dicts specifying
    conditions that produce the desired product.
    Returns a Boolean.
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
    Helper function that takes as input a dict with product SMILES as key values
    and a list of reactions as values. Each reaction is a dict containing 1 k,v pair,
    where the key is the reaction SMILES and the value is a list of dicts, each of
    which stores a set of dicts specifying conditions that produce the desired product.
    Returns a new dictionary that has the same format as the original input dict but
    doesn't contain any products with reaction steps within the synthesis pathway
    that don't contain valid conditions.
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


def transform_data(filtered_pathways):
    """
    Helper function that takes in the final output from SynthesisPathway and
    returns a tuple of 3 elements.
    1st element: Dictionary where each key is a reaction SMILES and each value
    is a dict storing the set of conditions with the highest probability of success.
    2nd element: Set of all the reaction strings.
    3rd element: Dictionary where each key is a reaction SMILES and each value is
    a list of the synthesis pathway it's in.
    """
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
                all_rxns.add(reaction_smiles)
                current_pathway.append(reaction_smiles)
        # print(f'current pathway: {current_pathway}')
        for reaction in current_pathway:
            rxns_to_pathways[reaction] = current_pathway
    return rxns_top_conditions, all_rxns, rxns_to_pathways


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
    # reactions_list = [(reaction, rxns_top_conditions[reaction]["temperature"])
    #                   for reaction in next_group]
    reactions_list = []
    for reaction in next_group:
        # try:
        rxn_and_temp_tuple = (reaction, rxns_top_conditions[reaction]["temperature"])
        reactions_list.append(rxn_and_temp_tuple)
        # except:
        #     continue

    # Sort the list of tuples by temperature
    sorted_reactions = sorted(reactions_list, key=lambda x: x[1])
    return sorted_reactions


def temperature_groups(next_group, rxns_top_conditions):
    """
    Helper function that takes in a list of unsorted uncompleted reactions.
    Returns a list of dicts, where each dict contains one key-value pair
    mapping a temperature range to a set of reactions within that range.
    """
    sorted_next_group = sort_rxns_by_temp(next_group, rxns_top_conditions)
    temp_groups = {}
    current_group = set()
    current_temp = None

    for rxn, temp in sorted_next_group:
        if current_temp is None:
            current_group.add(rxn)
            current_temp = temp
        elif (temp - current_temp) < 10:
            current_group.add(rxn)
        else:
            temp_groups[f'{current_temp}_{current_temp + 10}'] = current_group
            # reinitialize current_group to only contain the current reaction
            current_group = {rxn}
            # reinitialize minimum temp at start of new group
            current_temp = temp
    temp_groups[f'{current_temp}_{current_temp + 10}'] = current_group

    return temp_groups


def largest_temp_groups(temp_groups, num=5):
    """
    Helper function that returns up to the top [num] (default 5)
    largest temperature groups.
    temp_groups is a dict where each key is a temperature range and each
    value is a set of reactions in that range.
    Returns a list of dicts.
    """
    # Sorting the list of dictionaries by the length of the lists in descending order
    # sorted_dict_values = sorted(my_dict.items(), key=lambda item: item[1])

    sorted_list = sorted(temp_groups.items(), key=lambda item: len(item[1]), reverse=True)

    # Selecting the top five dictionaries with the largest lists
    top_five = sorted_list[:5]

    return top_five


def update_step(temp_groups, completed, uncompleted):
    """
    Helper function that selects the biggest temperature group and mutates
    the completed and uncompleted sets to reflect that these reactions
    have been performed.
    Returns a 4-elem tuple: a string of the temp range; a list of the most
    recently completed reactions; an updated completed set;
    an updated uncompleted set.
    """
    temp_range = None
    largest_group = []

    for key, value_list in temp_groups.items():
        if len(value_list) > len(largest_group):
            temp_range = key
            largest_group = value_list

    # largest_group = max(temp_groups.values(), key=len)
    completed.update(largest_group)
    uncompleted.difference_update(largest_group)
    return temp_range, list(largest_group), completed, uncompleted


def update_step_2(recent_group, completed, uncompleted):
    """
    Helper function that takes in a set of reaction SMILES to be performed,
    a set of completed reactions, and a of uncompleted reactions.
    Returns a 2-elem tuple of two new (not mutated) completed and uncompleted sets.
    """
    new_completed = recent_group | completed
    new_uncompleted = uncompleted - recent_group
    return new_completed, new_uncompleted


def get_conditions(recent_group, rxns_top_conditions):
    """
    Helper function that takes as input a list of reactions SMILES
    and a dict of top reaction conditions.
    Returns a list where each element is a dict with one key-value pair,
    the key being the reaction SMILES and the value being another dict
    containing the optimal reaction conditions.
    """
    recent_group_with_conditions = []
    for reaction in recent_group:
        rxn_with_conditions = {}
        rxn_with_conditions[reaction] = rxns_top_conditions[reaction]
        recent_group_with_conditions.append(rxn_with_conditions)
    return recent_group_with_conditions


# initially uncompleted is the same as all_rxns
def wellplate_sequence(filtered_pathways):
    """
    Returns a dictionary mapping a well plate ID to a list of reactions,
    each of which is a dictionary mapping the reaction to another dictionary
    with the reaction conditions.
    """
    # initialize variables
    rxns_top_conditions, uncompleted, rxns_to_pathways = transform_data(filtered_pathways)
    # print(f'rxns_top_conditions: {rxns_top_conditions}')
    # print(f'uncompleted: {uncompleted}')
    # print(f'rxns_to_pathways: {rxns_to_pathways}')
    completed = set()
    wellplates = {}
    plate_id = 0

    # continue sorting reactions into wellplates until uncompleted is empty
    while uncompleted:
        next_group = make_next_group(completed, uncompleted, rxns_to_pathways)
        # print(f'next_group: {next_group}')
        # temp_groups is a list of dicts??
        temp_groups = temperature_groups(next_group, rxns_top_conditions)
        # print(f'temp_groups: {temp_groups}')
        temp_range, recent_group, completed, uncompleted = update_step(temp_groups, completed, uncompleted)
        recent_group_with_conditions = get_conditions(recent_group, rxns_top_conditions)
        wellplate_key = f'{plate_id}_{temp_range}'
        wellplates[wellplate_key] = recent_group_with_conditions
        plate_id += 1

    return wellplates


def top_sequences(filtered_pathways, num):
    """
    Returns a list of "num" possible wellplate sequences, where each sequence is a
    dict. Each wellplate sequence is generated from one of the top "num" largest
    temp_groups from the first iteration.
    """
    # initialize variables
    rxns_top_conditions, uncompleted, rxns_to_pathways = transform_data(filtered_pathways)
    completed = set()
    # wellplates = {}
    top_seqs = []

    # first find the top "num" largest temp_groups from the first iteration
    next_group = make_next_group(completed, uncompleted, rxns_to_pathways)
    temp_groups = temperature_groups(next_group, rxns_top_conditions)
    top_temp_groups = largest_temp_groups(temp_groups, num)

    # make a separate path for each of these groups
    for top_temp_group in top_temp_groups:
        # initialize new wellplates dict for each temperature group
        plate_id = 0
        wellplates = {}
        recent_group = set()
        # top_temp_group is a tuple where 1st elem is the temp range
        # and 2nd elem is a set of reactions
        temp_range, recent_group = top_temp_group
        # be careful of mutation!! update_step_2 doesn't mutate the original completed and uncompleted sets
        new_completed, new_uncompleted = update_step_2(recent_group, completed, uncompleted)
        recent_group_with_conditions = get_conditions(recent_group, rxns_top_conditions)
        wellplate_key = f'{plate_id}_{temp_range}'
        wellplates[wellplate_key] = recent_group_with_conditions
        plate_id += 1

        # continue sorting reactions into wellplates until uncompleted is empty
        while new_uncompleted:
            next_group = make_next_group(new_completed, new_uncompleted, rxns_to_pathways)
            # print(f'next_group: {next_group}')
            # temp_groups is a list of dicts??
            temp_groups = temperature_groups(next_group, rxns_top_conditions)
            # print(f'temp_groups: {temp_groups}')
            temp_range, recent_group, new_completed, new_uncompleted = update_step(temp_groups, new_completed,new_uncompleted)
            recent_group_with_conditions = get_conditions(recent_group, rxns_top_conditions)
            wellplate_key = f'{plate_id}_{temp_range}'
            wellplates[wellplate_key] = recent_group_with_conditions
            plate_id += 1

        top_seqs.append(wellplates)

    return top_seqs


def find_all_paths(filtered_pathways):
    """
    Returns a list of all possible wellplate sequences, given that the top specified num of
    largest groups are expanded into paths at each iteration.
    """
    # first transform data
    rxns_top_conditions, uncompleted, rxns_to_pathways = transform_data(filtered_pathways)

    # initialize variables and start state
    wellplates = {}
    completed = set()
    start = (wellplates, completed, uncompleted)

    # initialize list to store paths that are being extended
    agenda = [start]
    agenda_path = 0

    # initialize list to store paths that are completed
    completed_paths = []

    # keep looping while there are paths with a non-empty uncompleted set
    while agenda:
        this_path = agenda.pop(0)
        wellplates, completed, uncompleted = this_path
        next_group = make_next_group(completed, uncompleted, rxns_to_pathways)
        temp_groups = temperature_groups(next_group, rxns_top_conditions)
        # a list of 2-elem tuples of the format (temp_range, {set of rxns})
        neighbors = largest_temp_groups(temp_groups)

        # try building paths with each of the neighbors
        for neighbor in neighbors:
            temp_range, recent_group = neighbor
            # this function returns new completed and uncompleted sets to avoid mutating the original
            new_completed, new_uncompleted = update_step_2(recent_group, completed, uncompleted)
            # make a copy to avoid mutation
            updated_wellplates = {k:v for k,v in this_path[0].items()}
            # returns a list where each element is a dict mapping a reation SMILES to a dict
            # of its optimal conditions
            recent_group_with_conditions = get_conditions(recent_group, rxns_top_conditions)
            new_plate_id = len(updated_wellplates) + 1
            wellplate_key = f'{new_plate_id}_{temp_range}'
            updated_wellplates[wellplate_key] = recent_group_with_conditions
            new_state = (updated_wellplates, new_completed, new_uncompleted)

            agenda_path += 1
            # print(f'agenda path: {agenda_path}')

            if new_uncompleted == set():
                completed_paths.append(new_state[0])
            else:
                agenda.append(new_state)

    return completed_paths


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
    for i in range(1, num+1):
        for k in wellplates.keys():
            if k.startswith(str(i)):
                rxns = wellplates[k]
                update_count = iterate_reactions(rxns, products)
                count += update_count
                break
    return count

def max_count_sequence(dir, filename, sequences, num, filtered_pathways):
    """
    Returns the ID of the wellplate sequence with the greatest number
    of products produced in the first "num" wellplates.
    """
    products = set(filtered_pathways.keys())
    product_counts = []

    for id in range(0, sequences+1):
        filepath = os.path.join(dir, filename)
        filepath += f'_{id}'
        with open(filepath, 'r') as jsonfile:
            wellplates = json.load(jsonfile)
        # print(f'wellplates: {wellplates}')
        count = count_products(wellplates, num, products)
        product_counts.append(count)

    id_of_greatest = product_counts.index(max(product_counts)) + 1
    max_count = max(product_counts)
    return id_of_greatest, max_count, product_counts


if __name__ == "__main__":
    # helper function test cases

    # rxns_top_conditions, uncompleted, rxns_to_pathways = transform_data(data)
    # completed = {'B(C1=CC=CC=C1)(O)O.Nc1ncccn1>>c1ccc(Nc2ncccn2)cc1'}
    # uncompleted.remove('B(C1=CC=CC=C1)(O)O.Nc1ncccn1>>c1ccc(Nc2ncccn2)cc1')
    # next_group = make_next_group(completed, uncompleted, rxns_to_pathways)
    # temp_groups = temperature_groups(next_group, rxns_top_conditions)
    # print(temp_groups)
    # temp_range, largest_group, completed, uncompleted = update_step(temp_groups, completed, uncompleted)
    # print(temp_range)
    # print(largest_group)
    # print(completed)
    # print(uncompleted)

    # largest temp groups
    # temp_groups = {
    # '15_25': {'rxn1', 'rxn2', 'rxn3'},
    # '28_38': {'rxn4', 'rxn5'},
    # '45_55': {'rxn6', 'rxn7', 'rxn8', 'rxn9'},
    # '60_70': {'rxn10', 'rxn11'},
    # '75_85': {'rxn12'},
    # '90_100': {'rxn13', 'rxn14', 'rxn15', 'rxn16', 'rxn17'}
    # }

    # result = largest_temp_groups(temp_groups)
    # print(f'result: {result}')
    # result = [('90_100', {'rxn17', 'rxn14', 'rxn16', 'rxn15', 'rxn13'}),
    #         ('45_55', {'rxn9', 'rxn6', 'rxn8', 'rxn7'}),
    #         ('15_25', {'rxn1', 'rxn3', 'rxn2'}),
    #         ('28_38', {'rxn5', 'rxn4'}),
    #         ('60_70', {'rxn11', 'rxn10'})]

    # main function test case
    # wellplates = wellplate_sequence(data)
    # print(wellplates)

    # MFBO mols
    # inp_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/MFBO_contexts_and_preds/MFBO_selected_mols_contexts_and_preds.json'
    # inp_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/MFBO_selected_mols_pathways_w_conditions.json'
    # with open(inp_filepath, 'r') as jsonfile:
    #     pathways_with_conditions = json.load(jsonfile)
    # print(f'number of products initially: {len(list(pathways_with_conditions.keys()))}')
    # filtered_pathways = filter_data(pathways_with_conditions)
    # out_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/MFBO_selected_mols_filtered_pathways.json'
    # with open(out_filepath, 'w') as outfile:
    #     json.dump(filtered_pathways, outfile, indent=4)
    # print(f'number of products after: {len(list(filtered_pathways.keys()))}')

    # wellplates time!!

    inp_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/MFBO_selected_mols_filtered_pathways.json'
    with open(inp_filepath, 'r') as jsonfile:
        filtered_pathways = json.load(jsonfile)

    # inp_filepath = "/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/test_cases/pathfinding_small.json"
    # with open(inp_filepath, 'r') as jsonfile:
    #     filtered_pathways = json.load(jsonfile)

    top_seqs = top_sequences(filtered_pathways, 5)
    for i, top_seq in enumerate(top_seqs):
        products = set(filtered_pathways.keys())
        count = count_products(top_seq, 9, products)
        print(f'{i}: {count}')
    # print(f'top_seqs: {top_seqs}')

    # rxns_top_conditions, all_rxns, rxns_to_pathways = transform_data(filtered_pathways)
    # print(f'rxns_top_conditions: {rxns_top_conditions[]}')
    # print(f'all_rxns: {all_rxns}')
    # print(f'rxns_to_pathways: {rxns_to_pathways}')

    # print('starting wellplates')
    # wellplates_dict = wellplate_sequence(filtered_pathways)
    # out_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/MFBO_selected_mols_wellplates.json'
    # with open(out_filepath, 'w') as outfile:
    #     json.dump(wellplates_dict, outfile, indent=4)
    # print(f'File saved to {out_filepath}')

    # find all paths
    # completed_paths = find_all_paths(filtered_pathways)
    # print(completed_paths)


    # dir = "/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/test_cases/"
    # max = max_count_sequence(dir, "pathfinding_small_path", 44, 5, filtered_pathways)
    # print(max)


    # save_completed_paths(completed_paths, dir, "pathfinding_small_path")

    def total_elements(dictionary):
        total = 0
        for key, value in dictionary.items():
            total += len(value)
        return total

    # print(f'total rxns: {total_elements(wellplates_dict)}')
