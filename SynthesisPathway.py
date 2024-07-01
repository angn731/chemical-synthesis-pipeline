from rdkit import Chem
from rdkit.Chem import AllChem
import json
import ForwardPred

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
        product = products[0][0]
        product_smiles = Chem.MolToSmiles(product)
    return product_smiles

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
    templates = data["templates"]
    reactant1 = data["initial_mol"]
    rxn_pathway = ()
    # loop through all the reaction templates, where each template is a dict
    for temp in templates:
        template_id = temp["_id"]
        # assume there is only one reactant
        reactant2 = temp["reactants"][0]
        intermediate = react(reactant1, reactant2, template_id)
        rxn_smiles = reactant1 + '.' + reactant2 + ">>" + intermediate
        rxn_pathway += (rxn_smiles,)
        # update reactant 1 to react the intermediate by the previous step
        # with the reactant in the next step
        reactant1 = intermediate
    # the final product should be the desired product
    if reactant1 != product:
        raise ValueError("The final product does not match the desired product.")
    return rxn_pathway

# forward pred took a json file that was a list of reaction SMILES strings, then output a dict
# where the keys were reaction SMILES strings and the values were a list of dicts with conditions and
# the probability of generating the desired product


# ////////////////////////////////////////////
# GENERATING SYNTHESIS PATHWAYS AND CONDITIONS
# ////////////////////////////////////////////



# can refactor code below into a separate function

# load json file with products and their associated data
products_path = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/iteration_4.json'

with open(products_path, 'r') as file:
    products_data = json.load(file)

# store synthesis pathways in a dict
reaction_pathways = {}

# simultaneously store a list of all the reactions
all_reactions = []
for prod, data in products_data.items():
    # rxn_pathway is a tuple of reaction strings
    reaction_pathways[prod] = synthesis_pathway(prod, data)
    rxns = list(reaction_pathways[prod])
    all_reactions.extend(rxns)

print(f'all_reactions: {all_reactions}')
print(f'reaction_pathways: {reaction_pathways}')

# feed the list of reactions into the forward pred
rxn_contexts_and_preds = ForwardPred.contexts_and_preds(all_reactions)
print('finished predicting contexts and top products')

file_path = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/iteration_4_contexts_and_preds.json'

with open(file_path, 'r') as file:
    rxn_contexts_and_preds = json.load(file)

# only keep conditions where the highest probability product matches the desired product
valid_conditions = ForwardPred.compare_products(rxn_contexts_and_preds)
print('finished comparing top products to desired product')
# print(valid_conditions)

# reaction_pathways is a dict with desired product mapped to a tuple of SMILES strings
# representing reactions
# loop through all of the tuples and grab a list of dicts of conditions from
# valid_conditions
pathways_with_conditions = {}
for product, pathway in reaction_pathways.items():
    pathway_with_conditions = tuple()
    for reaction in pathway:
        reaction_dict = {}
        reaction_dict[reaction] = valid_conditions[reaction]
        pathway_with_conditions += (reaction_dict,)
    pathways_with_conditions[product] = pathway_with_conditions

print('finished adding conditions to pathways')
print(pathways_with_conditions)

pathways_path = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/small_test_case_pathways_and_conditions.json'
with open(pathways_path, 'w') as outfile:
    json.dump(pathways_with_conditions, outfile, indent=4)
print(f"File saved to {pathways_path}")

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

    pass
