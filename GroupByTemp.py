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
    returns a dictionary where each key is a reaction SMILES and each value
    is a dict storing the set of conditions with the highest probability of success.
    Also returns a set of all the reaction strings.
    """
    rxns_top_conditions = {}
    all_rxns = set()
    # pathway is a list
    for product, pathway in pathways_with_conditions.items():
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
    return rxns_top_conditions, all_rxns

# create a helper function that pulls the next group of reactions to examine
def temperature_groups(completed, uncompleted, rxns_top_conditions):
    """
    Helper function that takes in a set of completed reactions, a set of uncompleted
    reactions, and a dict storing the top reaction conditions for each reaction.
    Returns a list, where each element is a set of reactions that can be completed
    together next.
    """


if __name__ == "__main__":
    # transform data

    # load in file from synthesis_pathway(); need to store file in SynethsisPathway
    pass

{
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
    ]
