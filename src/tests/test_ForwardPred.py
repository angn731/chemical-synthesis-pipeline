from src.ForwardPred import react, synthesis_pathway, make_synthesis_pathways, read_file


def test_react():
    """
    Test case for react() with two reactants and a template ID.
    """
    reactant1 = "C1=CC(=C(N=C1)N)Br"
    reactant2 = "O=C(O)c1ccc(B(O)O)cc1"
    template_id = "ChanLam"
    result = react((reactant1, reactant2), template_id, check_prods=False)
    print(f"Result: {result}")


def test_synthesis_pathway():
    """
    Test case for synthesis_pathway() with a product and its data.
    """
    product = "O=S(=O)(N(c1ccccc1)c1ncccn1)C(F)(F)F"
    data = {
        "QED": ["0.8735248151923747", 0],
        "similarity": ["0.579330638082168", 0],
        "rank": 3,
        "initial_mol": "B(C1=CC=CC=C1)(O)O",
        "templates": [
            {"_id": "ChanLam", "reactants": ["Nc1ncccn1"]},
            {"_id": "Sulfonamide", "reactants": ["O=S(=O)(Cl)C(F)(F)F"]}
        ]
    }
    result = synthesis_pathway(product, data)
    print(f"Synthesis Pathway Result: {result}")


def test_make_synthesis_pathways():
    """
    Test case for make_synthesis_pathways() with a sample JSON file.
    """
    products_path = './gen_mols/gen_mols_Jun2024.json'
    products_data = read_file(products_path)
    all_reactions, reaction_pathways = make_synthesis_pathways(products_data)
    print("All Reactions:", all_reactions)
    print("Reaction Pathways:", reaction_pathways)


def test_intermediate_rxn():
    """
    Test case to ensure intermediate reaction in synthesis pathway returns None when it fails.
    """
    molecule = {
        "Cc1ccc(S(=O)(=O)N(c2ccc(-n3cc(CCCO)nn3)cc2)c2ccccc2C)cc1": {
            "HDAC_Docking": ["0.002662945346571803", "-1.26"],
            "LogP": ["4.39157551055326", "0.116244274465366"],
            "Toxicity": ["1.383194056838727", "0.032358073478579774"],
            "similarity": ["0.41583626050085826", 0],
            "rank": 4,
            "initial_mol": "Cc1ccccc1F",
            "templates": [
                {"_id": "SnAr_ForCl", "reactants": ["Nc1ccc(N)cc1"]},
                {"_id": "Sulfonamide", "reactants": ["Cc1ccc(S(=O)(=O)Cl)cc1"]},
                {"_id": "ClickChem_aryl_amine2azide", "reactants": ["C#CCCCO"]}
            ]
        }
    }
    product = "Cc1ccc(S(=O)(=O)N(c2ccc(-n3cc(CCCO)nn3)cc2)c2ccccc2C)cc1"
    data = molecule[product]
    pathway_result = synthesis_pathway(product, data)
    full_result = make_synthesis_pathways(molecule)
    full_expected = ([], {})
    assert pathway_result is None
    assert full_result == full_expected


if __name__ == "__main__":
    test_react()
    test_synthesis_pathway()
    test_make_synthesis_pathways()
    test_intermediate_rxn()
    print("All test cases completed.")
