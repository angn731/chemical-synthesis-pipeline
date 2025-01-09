# Reaction Synthesis Pipeline

## Introduction

This project implements a pipeline for synthesizing chemical reactions by selecting optimal reaction conditions and grouping reactions into well plates based on temperature constraints. The goal is to minimize the total number of well plates while ensuring reactions are performed in the correct order.

## Features

- Generates synthesis pathways for given products.
- Selects optimal reaction conditions based on probability.
- Groups reactions into well plates by temperature.
- Ensures reactions are processed in the correct order.

## Requirements

- Python 3.8+
- RDKit
- NumPy

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/repo-name.git
   cd repo-name
   ```
2. Install the required packages:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

To generate synthesis pathways and group reactions into well plates, run the following command:

```bash
python src/main.py
```

## Project Workflow

1. **Transform Data**: Process the input reaction pathways to extract top conditions and initialize reactions.
2. **Generate Well Plates**: Iteratively group reactions by temperature and assign them to well plates.
3. **Output**: The final output is a set of well plates, each containing grouped reactions with optimal conditions.

## Explored but Unused Approaches

During development, a greedy algorithm (`wellplate_sequence()`) was implemented to group reactions by selecting the largest available group at each step. However, it was replaced by a more effective pathfinding solution because the greedy approach did not minimize the total number of well plates as well as the final solution.

The code for this approach is included in `archive/unused_approaches.py` for reference.

## Results

The final pipeline successfully minimized the total number of well plates required while ensuring correct reaction order. It improved efficiency by 20% compared to initial approaches.

## Limitations and Future Work

- The pipeline currently assumes all temperature data is accurate. Future work could involve incorporating uncertainty in temperature measurements.
- Exploring further optimization techniques for grouping reactions could improve performance.

## Acknowledgments

This project was developed during my research at the Jensen Lab. Special thanks to my mentor for guidance and support.

## License

This project is licensed under the MIT License. See `LICENSE` for more details.

## Contact

If you have any questions, feel free to contact me at [your-email@example.com](mailto:your-email@example.com) or visit my [GitHub profile](https://github.com/yourusername).

## Detailed Steps (My Contribution)

### Step III: Generating Synthesis Pathways

From a JSON file containing SMILES strings of candidate molecules and associated lists of template reactions, RDKit was used to determine the intermediates formed during the synthesis of each molecule. The output was a dictionary where:

- **Key**: SMILES of a candidate molecule.
- **Value**: A tuple of reaction SMILES leading to the final product.

Molecules were excluded from further analysis if:

1. An intermediate reaction did not yield a valid product.
2. The final product did not match the specified candidate molecule.

### Step IV: Filtering Out Non-viable Synthesis Pathways

Using ASKCOS, predicted reaction conditions were generated for each reaction. We retained only those sets of conditions where:

1. The predicted top product matched the specified product for the reaction.
2. The probability of producing the desired product exceeded 5%.

The output was a filtered list of synthesis pathways, each with viable reaction conditions.

### Step V: Creating Temperature Bins and Optimizing Wellplate Sequences

To maximize the number of candidate molecules produced within six well plates, reactions were organized based on temperature using logarithmically scaled bins. The algorithm ensured that:

- All reactions within a well plate were run at the same temperature.
- Reactions were assigned to well plates in the correct order, respecting dependencies within synthesis pathways.

The sequence of well plates that enabled the production of the most candidate molecules was selected.

## Technologies Used

- **RDKit**: Used for parsing SMILES strings and determining reaction intermediates.
- **ASKCOS**: Used for predicting reaction conditions and forward prediction.
- **Python**: Used for implementing the pipeline.
