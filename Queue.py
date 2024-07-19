# imports
import json
import os

# goal: for each wellplate, generate a list where each
# element is a dict representing one rxn


def get_rxn_data(rxn_smiles, conditions, pathway):
    """
    Conditions is a dict containing information about the specified rxn.
    Pathways is a dict that maps each rxn SMILES to a tuple of the pathway
    it's in.
    Returns a dictionary containing information about the reaction
    for the queue.
    """
    # extract information from pathways
    # pathway = pathways[rxn_smiles]
    rxn_steps = len(pathway)
    rxn_step = None
    pathway_product = None
    for i, rxn_in_pathway in enumerate(pathway):
        if rxn_smiles == rxn_in_pathway:
            rxn_step = i + 1
        if i == len(pathway) - 1:
            pathway_product = rxn_in_pathway.split('>>')[1]

    # store information about rxn in dict
    collect_product = ''
    if conditions['final product']:
        collect_product = 'yes'
    else:
        collect_product = 'no'
    rxn_data = {}
    final_product = [rxn_step, collect_product,
                     pathway_product, rxn_steps]
    predicted_catalysts = [conditions['catalyst']]
    predicted_reactants = rxn_smiles.split('>>')[0].split('.')
    predicted_reagents = [conditions['reagent']]
    predicted_solvents = [conditions['solvent']]
    reaction_smiles = rxn_smiles
    target_product = [[rxn_smiles.split('>>')[1]],]
    total_volume = 500
    rxn_data = {
        'final_product': final_product,
        'predicted_catalysts': predicted_catalysts,
        'predicted_reactants': predicted_reactants,
        'predicted_reagents': predicted_reagents,
        'predicted_solvents': predicted_solvents,
        'reaction_smiles': reaction_smiles,
        'target_product': target_product,
        'total_volume': total_volume
    }
    return rxn_data


def get_multiple_rxn_data(rxn_smiles, conditions, pathways):
    """
    Returns a list of data for a particular reaction.
    """
    rxn_data = []
    for pathway in pathways[rxn_smiles]:
        pathway_data = get_rxn_data(rxn_smiles, conditions, pathway)
        rxn_data.append(pathway_data)
    return rxn_data


def get_plate_queue_data(rxns, pathways):
    """
    rxns is a list where each element is a dict with 1 k,v pair. Key is
    a rxn SMILES and value is a dict of conditions.
    Returns a list of dictionaries, containing data about each reaction
    for the queue.
    """
    plate_queue_data = []
    # rxn is a dict with 1 k,v pair
    for rxn in rxns:
        for rxn_smiles, conditions in rxn.items():
            rxn_data = get_multiple_rxn_data(rxn_smiles, conditions, pathways)
            # extend instead of append
            plate_queue_data.extend(rxn_data)
    return plate_queue_data

if __name__ == "__main__":
    # seqs is a list of dicts
    # rxns is a list
    filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/new_rxns_to_pathways.json'
    with open(filepath, 'r') as jsonfile:
        pathways = json.load(jsonfile)

    # print(f'pathways: {pathways['Nn1cnnc1.COc1ccc(B(O)O)cc1>>COc1ccc(Nn2cnnc2)cc1']}')
    # print(f'pathways: {pathways['COc1ccc(Nn2cnnc2)cc1.Cc1noc(C)c1B(O)O>>COc1ccc(N(c2c(C)noc2C)n2cnnc2)cc1']}')
    # first_two_items = {k: pathways[k] for k in list(pathways)[:2]}
    # print(first_two_items)

    seq_filepath = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/best_six_plate_seq/MFBO_selected_mols_six_plate_seq_product_key.json'
    with open(seq_filepath, 'r') as jsonfile:
        seqs = json.load(jsonfile)

    # print(seqs)
    plate_ids = {'0_', '1_', '2_', '3_', '4_', '5_'}
    for plate, rxns in seqs.items():
        if plate[0:2] in plate_ids:
            plate_queue_data = get_plate_queue_data(rxns, pathways)
            dir = '/Users/angelinaning/Downloads/jensen_lab_urop/reaction_pathways/reaction_pathways_code/MFBO_selected_mols/queue_data'
            filename = f'new_{plate}_queue.json'
            filepath = os.path.join(dir, filename)

        with open(filepath, 'w') as outfile:
            json.dump(plate_queue_data, outfile, indent=4)
