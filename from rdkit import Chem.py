from rdkit import Chem
from rdkit.Chem import AllChem

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
