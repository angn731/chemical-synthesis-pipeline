# imports


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

def transform_data(pathways_with_conditions):
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
    for product, pathway in pathways_with_conditions.items():
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
    reactions_list = [(reaction, rxns_top_conditions[reaction]["temperature"])
                      for reaction in next_group]

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
            temp_groups[f'{current_temp} to {current_temp + 10}'] = current_group
            # reinitialize current_group to only contain the current reaction
            current_group = {rxn}
            # reinitialize minimum temp at start of new group
            current_temp = temp
    temp_groups[f'{current_temp} to {current_temp + 10}'] = current_group

    return temp_groups


def update_step(temp_groups, completed, uncompleted):
    """
    Helper function that selects the biggest temperature group and mutates
    the completed and uncompleted sets to reflect that these reactions
    have been performed.
    Returns a 3-elem tuple: a list of the most recently completed reactions;
    updated completed set; updated uncompleted set.
    """
    largest_group = max(temp_groups.values(), key=len)
    completed.update(largest_group)
    uncompleted.difference_update(largest_group)
    return list(largest_group), completed, uncompleted


# initially uncompleted is the same as all_rxns
def wellplate_sequence(pathways_with_conditions):
    """
    Returns a dictionary mapping a well plate ID to a list of reactions,
    each of which is a dictionary mapping the reaction to another dictionary
    with the reaction conditions.
    """
    # initialize variables
    rxns_top_conditions, uncompleted, rxns_to_pathways = transform_data(pathways_with_conditions)
    completed = set()
    wellplates = {}
    plate_id = 0

    # continue sorting reactions into wellplates until uncompleted is empty
    while uncompleted:
        next_group = make_next_group(completed, uncompleted, rxns_to_pathways)
        # temp_groups is a dict
        temp_groups = temperature_groups(next_group, rxns_top_conditions)
        recent_group, completed, uncompleted = update_step(temp_groups, completed, uncompleted)
        recent_group_with_conditions = []
        for reaction in recent_group:
            rxn_with_conditions = {}
            rxn_with_conditions[reaction] = rxns_top_conditions[reaction]
            recent_group_with_conditions.append(rxn_with_conditions)
        wellplates[plate_id] = recent_group_with_conditions
        plate_id += 1

    return wellplates


if __name__ == "__main__":
    # transform data
    data = {
    "O=C(O)c1ccc(Nc2ncccc2Br)cc1": [
        {
            "C1=CC(=C(N=C1)N)Br.O=C(O)c1ccc(B(O)O)cc1>>O=C(O)c1ccc(Nc2ncccc2Br)cc1": [
                {
                    "temperature": 28.327255249023438,
                    "solvent": "ClCCl",
                    "reagent": "c1ccncc1",
                    "catalyst": "CC(=O)[O-].[Cu+2]",
                    "prob": 0.001115021812759695
                },
                {
                    "temperature": 41.79078674316406,
                    "solvent": "ClCCl",
                    "reagent": "",
                    "catalyst": "CC(=O)[O-].[Cu+]",
                    "prob": 4.193923258137252e-06
                },
                {
                    "temperature": 38.532230377197266,
                    "solvent": "ClCCl",
                    "reagent": "",
                    "catalyst": "CC(=O)[O-].[Cu+2]",
                    "prob": 4.193923258137252e-06
                },
                {
                    "temperature": 28.02731704711914,
                    "solvent": "ClCCl",
                    "reagent": "CCN(CC)CC",
                    "catalyst": "CC(=O)[O-].[Cu+2]",
                    "prob": 9.984051611924022e-06
                },
                {
                    "temperature": 35.62443161010742,
                    "solvent": "ClCCl",
                    "reagent": "c1ccncc1",
                    "catalyst": "CC(=O)[O-].[Cu+]",
                    "prob": 0.001115021812759695
                }
            ]
        }
    ],
    "O=S(=O)(N(c1ccccc1)c1ncccn1)C(F)(F)F": [
        {
            "B(C1=CC=CC=C1)(O)O.Nc1ncccn1>>c1ccc(Nc2ncccn2)cc1": [
                {
                    "temperature": 26.51129913330078,
                    "solvent": "ClCCl",
                    "reagent": "CCN(CC)CC",
                    "catalyst": "CC(=O)[O-].[Cu+2]",
                    "prob": 0.605522485929338
                },
                {
                    "temperature": 30.631725311279297,
                    "solvent": "ClCCl",
                    "reagent": "c1ccncc1",
                    "catalyst": "CC(=O)[O-].[Cu+2]",
                    "prob": 0.434369254770806
                },
                {
                    "temperature": 88.17671203613281,
                    "solvent": "Cc1ccccc1",
                    "reagent": "O=P([O-])([O-])[O-].[K+]",
                    "catalyst": "",
                    "prob": 0.019987005521394673
                },
                {
                    "temperature": 60.11564636230469,
                    "solvent": "O",
                    "reagent": "O=C([O-])[O-].[K+]",
                    "catalyst": "",
                    "prob": 0.015679980018138476
                },
                {
                    "temperature": 56.967620849609375,
                    "solvent": "",
                    "reagent": "c1ccncc1",
                    "catalyst": "CC(=O)[O-].[Cu+2]",
                    "prob": 0.16607293981531082
                }
            ]
        },
        {
            "c1ccc(Nc2ncccn2)cc1.O=S(=O)(Cl)C(F)(F)F>>O=S(=O)(N(c1ccccc1)c1ncccn1)C(F)(F)F": [
                {
                    "temperature": 13.223930358886719,
                    "solvent": "ClCCl",
                    "reagent": "CCN(CC)CC",
                    "catalyst": "",
                    "prob": 0.8094813695538945
                },
                {
                    "temperature": 17.88831329345703,
                    "solvent": "ClC(Cl)Cl",
                    "reagent": "CCN(CC)CC",
                    "catalyst": "",
                    "prob": 0.8533710279069024
                },
                {
                    "temperature": 21.842477798461914,
                    "solvent": "CC#N",
                    "reagent": "CCN(CC)CC",
                    "catalyst": "",
                    "prob": 0.8756587644695389
                },
                {
                    "temperature": 27.828472137451172,
                    "solvent": "CC#N",
                    "reagent": "O=C([O-])[O-].[K+]",
                    "catalyst": "",
                    "prob": 0.34261022868695973
                },
                {
                    "temperature": 16.287261962890625,
                    "solvent": "ClCCl",
                    "reagent": "",
                    "catalyst": "",
                    "prob": 0.29794359498509
                }
            ]
        }
    ]}

    # helper function test cases

    # rxns_top_conditions, uncompleted, rxns_to_pathways = transform_data(data)
    # completed = {'B(C1=CC=CC=C1)(O)O.Nc1ncccn1>>c1ccc(Nc2ncccn2)cc1'}
    # uncompleted.remove('B(C1=CC=CC=C1)(O)O.Nc1ncccn1>>c1ccc(Nc2ncccn2)cc1')
    # next_group = make_next_group(completed, uncompleted, rxns_to_pathways)
    # temp_groups = temperature_groups(next_group, rxns_top_conditions)
    # largest_group, completed, uncompleted = update_step(temp_groups, completed, uncompleted)
    # print(largest_group)
    # print(completed)
    # print(uncompleted)

    # main function test case
    # wellplates = wellplate_sequence(data)
    # print(wellplates)

    # load in file from synthesis_pathway(); need to store file in SynethsisPathway
    # pass
