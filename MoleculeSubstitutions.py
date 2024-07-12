# imports
import json
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, DataStructs

def cas_to_smiles(cas_nums):
    """
    Takes in a list of CAS numbers as strings and returns a list
    of SMILES.
    """
    smiles_list = []
    for i, cas in enumerate(cas_nums):
        print(f'on cas {i}')
        compound = pcp.get_compounds(cas, 'name')
        if compound:
            smiles = compound[0].canonical_smiles
            smiles_list.append(smiles)
        else:
            smiles_list.append(None)
    return smiles_list


def morgan_fingerprints(molecules):
    """
    Given a list of molecule SMILES, returns a list of their Morgan fingerprints.
    """
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
    molecule_fps = []
    for smiles in molecules:
        molecule_fps.append(mfpgen.GetFingerprint(Chem.MolFromSmiles(smiles)))
    return molecule_fps


def tanimoto(missing_mols, reagents):
    """
    Given a list of SMILES of missing molecules and a list of SMILES of reagents,
    returns a dict where each key is a missing molecule and each value is a tuple where
    the first element is a list of all the tanimoto values between that molecule and
    each reagent and the second element is the SMILES of the most similar reagent.
    """
    reagent_fps = morgan_fingerprints(reagents)
    missing_fps = morgan_fingerprints(missing_mols)
    tanimoto_values = {}
    closest_reagents = {}

    for i, fp in enumerate(missing_fps):
        tanimoto_similarity = DataStructs.BulkTanimotoSimilarity(fp, reagent_fps)
        max_value = max(tanimoto_similarity)
        closest_reagent = reagents[tanimoto_similarity.index(max_value)]
        tanimoto_values[missing_mols[i]] = tanimoto_similarity
        closest_reagents[missing_mols[i]] = (closest_reagent, max_value)

    return tanimoto_values, closest_reagents


if __name__ == "__main__":

    filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/reagents_smiles.json'
    with open(filepath, 'r') as jsonfile:
        reagents = json.load(jsonfile)

    filepath_2 = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/missing_mols_smiles.json'
    with open(filepath_2, 'r') as jsonfile:
        missing_mols = json.load(jsonfile)

    filepath_3 = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/no_color_reagents_smiles.json'
    with open(filepath_3, 'r') as jsonfile:
        no_color_reagents = json.load(jsonfile)

    reagents = set(reagents)
    missing_mols = set(missing_mols)
    no_color_reagents = set(no_color_reagents)
    all_reagents = reagents | no_color_reagents
    all_reagents = all_reagents - missing_mols

    all_reagents = list(all_reagents)
    missing_mols = list(missing_mols)


    tanimoto_values, closest_reagents = tanimoto(missing_mols, all_reagents)

    out_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/missing_mols_tanimoto_all.json'
    with open(out_filepath, 'w') as outfile:
        json.dump(tanimoto_values, outfile, indent=4)

    out_filepath_2 = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/missing_mols_closest_reagents_all.json'
    with open(out_filepath_2, 'w') as outfile:
        json.dump(closest_reagents, outfile, indent=4)
