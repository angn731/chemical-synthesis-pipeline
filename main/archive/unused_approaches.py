"""
This file contains alternative approaches that were explored during the development of the reaction synthesis pipeline.
These approaches were ultimately replaced by more effective solutions but are included here for reference.
"""

# imports
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))
from GroupByTemp import transform_data, make_next_group, temperature_groups, update_step, get_conditions

# /////////////////////////////////////
# Greedy Approach: Wellplate Sequence
# /////////////////////////////////////

def update_step(temp_groups, completed, uncompleted):
    """
    Helper function that returns up to the top `num` (default 5) largest temperature groups.

    Parameters:
    - `temp_groups`: Dictionary where keys are temperature ranges (strings) and values are sets of reactions within those ranges.
    - `num`: Maximum number of top temperature groups to return (default is 5).

    Returns:
    - A list of dictionaries representing the top temperature groups.
    """
    temp_range = None
    largest_group = []

    for key, value_list in temp_groups.items():
        if len(value_list) > len(largest_group):
            temp_range = key
            largest_group = value_list

    completed.update(largest_group)
    uncompleted.difference_update(largest_group)
    return temp_range, list(largest_group), completed, uncompleted


def wellplate_sequence(filtered_pathways):
    """
    Returns a dictionary mapping a well plate ID to a list of reactions,
    each of which is a dictionary mapping the reaction to another dictionary
    with the reaction conditions.

    This approach uses a greedy algorithm to iteratively select the largest
    possible group of reactions that can be performed at each step, grouped by temperature.
    """
    # Initialize variables
    rxns_top_conditions, uncompleted, rxns_to_pathways = transform_data(filtered_pathways)
    completed = set()
    wellplates = {}
    plate_id = 0

    # Continue sorting reactions into wellplates until uncompleted is empty
    while uncompleted:
        next_group = make_next_group(completed, uncompleted, rxns_to_pathways)
        temp_groups = temperature_groups(next_group, rxns_top_conditions)
        temp_range, recent_group, completed, uncompleted = update_step(temp_groups, completed, uncompleted)
        recent_group_with_conditions = get_conditions(recent_group, rxns_top_conditions)
        wellplate_key = f'{plate_id}_{temp_range}'
        wellplates[wellplate_key] = recent_group_with_conditions
        plate_id += 1

    return wellplates

# Note: This approach was explored but not used in the final pipeline because it did not
# perform well at maximizing the number of products generated in as few plates as possible.
