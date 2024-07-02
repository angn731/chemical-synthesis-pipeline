import pytest
import SynthesisPathway

import sys
sys.setrecursionlimit(10000)

def test_intermediate_rxn():
    """
    Intermediate reaction in synthesis pathway doesn't work -> return None.
    """
    molecule = {"Cc1ccc(S(=O)(=O)N(c2ccc(-n3cc(CCCO)nn3)cc2)c2ccccc2C)cc1":
                {"HDAC_Docking": ["0.002662945346571803", "-1.26"],
                 "LogP": ["4.39157551055326", "0.116244274465366"],
                 "Toxicity": ["1.383194056838727", "0.032358073478579774"],
                 "similarity": ["0.41583626050085826", 0],
                 "rank": 4,
                 "initial_mol": "Cc1ccccc1F",
                 "templates": [{"_id": "SnAr_ForCl", "reactants": ["Nc1ccc(N)cc1"]},
                               {"_id": "Sulfonamide", "reactants": ["Cc1ccc(S(=O)(=O)Cl)cc1"]},
                               {"_id": "ClickChem_aryl_amine2azide", "reactants": ["C#CCCCO"]}]
                    }
                }
    lab_result = SynthesisPathway.make_synthesis_pathways(molecule)
    assert lab_result == None
