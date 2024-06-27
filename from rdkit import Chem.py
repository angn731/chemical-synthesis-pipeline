from rdkit import Chem
from rdkit.Chem import AllChem
import json

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
