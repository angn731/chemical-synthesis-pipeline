from rdkit import Chem
from rdkit.Chem import AllChem
import json

# Example reactant SMILES strings
reactant1_smiles = 'CC(=O)O'
reactant2_smiles = 'CC(C)(C)C#N'

# Convert SMILES to RDKit molecule objects
reactant1 = Chem.MolFromSmiles(reactant1_smiles)
reactant2 = Chem.MolFromSmiles(reactant2_smiles)

# Example SMARTS string for Chan-Lam coupling
reaction_smarts = '[C:1]-[N:2].[C:3]-[Br,I,F,Cl:4]>>[C:1]-[N:2]-[C:3]'

# Create an RDKit reaction object from SMARTS
rxn = AllChem.ReactionFromSmarts(reaction_smarts)

# Run the reaction to predict the product
products = rxn.RunReactants((reactant1, reactant2))
print(products)

# Assuming the reaction produces one product, convert it to SMILES
if products:
    product = products[0][0]  # Assuming one product for simplicity
    product_smiles = Chem.MolToSmiles(product)
    print(f'Predicted Product SMILES: {product_smiles}')
else:
    print('No product generated.')

# output: dictionary, where keys are SMILES strings representing reaction products and the values are tuples
# each element of the tuple is a dictionary produced by the forward prediction code.
# so first element would be the first reaction, second element would be the second reaction, etc.

templates_dict_path = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/templates_dict.json'

with open(templates_dict_path, 'r') as file:
    rxn_templates = json.load(file)

def react(reactant1_smiles, reactant2_smiles, template_id):
    """
    Helper function that returns a SMILES of the product,
    given two reactant SMILES and a reaction template ID.
    """
    reactant1 = Chem.MolFromSmiles(reactant1_smiles)
    reactant2 = Chem.MolFromSmiles(reactant2_smiles)
    reaction_smarts =
