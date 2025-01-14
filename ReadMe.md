# Reaction Synthesis Pipeline

## Introduction

This project implements a pipeline for synthesizing chemical reactions by selecting optimal reaction conditions and grouping reactions into well plates based on temperature constraints.

## Pipeline

### Step I: Generating Synthesis Pathways

From a JSON file containing SMILES strings of candidate molecules and associated lists of template reactions, RDKit was used to determine the intermediates formed during the synthesis of each molecule. The output was a dictionary where:

- **Key**: SMILES of a candidate molecule.
- **Value**: A tuple of reaction SMILES leading to the final product.

Molecules were excluded from further analysis if:

1. An intermediate reaction did not yield a valid product.
2. The final product did not match the specified candidate molecule.

### Step II: Filtering Out Non-viable Synthesis Pathways

Using ASKCOS, predicted reaction conditions were generated for each reaction. I retained only sets of conditions where:

1. The predicted top product matched the specified product for the reaction.
2. The probability of producing the desired product exceeded 5%.

The output was a filtered list of synthesis pathways, each with viable reaction conditions.

### Step III: Creating Temperature Bins and Optimizing Wellplate Sequences

To maximize the number of candidate molecules produced within six well plates, reactions were organized based on temperature using logarithmically scaled bins. The algorithm ensured that:

- All reactions within a well plate were run at the same temperature.
- Reactions were assigned to well plates in the correct order, respecting dependencies within synthesis pathways.

The sequence of well plates that enabled the production of the most candidate molecules was selected.

## Technologies Used

- **RDKit**: Used for parsing SMILES strings and determining reaction intermediates.
- **ASKCOS**: Used for predicting reaction conditions and forward prediction.
- **Python**: Used for implementing the pipeline.
