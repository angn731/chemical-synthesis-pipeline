import GroupByTemp
from GroupByTemp import transform_data
# current algorithm: gets all the molecules that can be run together in the next step (basically a DAG situation)
# groups these rxns into temperature groups, then selects the biggest group; keeps doing this until all rxns completed

# new algorithm: gets all the molecules that can be run together in the next step and groups these rxns by temp. again
# but this time searches the top 5 largest temp. groups

# states: need to store all groups of rxns completed at each step, as well as completed and uncompleted sets from the
# previous step; could use a 3-elem tuple where the first element is a dict storing the wellplates up to that point
# the second element is the updated completed set and the third element is the updated uncompleted set

# neighbors: 5 next largest temp. groups

# goal test: uncompleted is empty

def find_all_paths(filtered_pathways, neighbors_function, start, goal_test):
    all_paths = []
    # first transform data
    rxns_top_conditions, uncompleted, rxns_to_pathways = transform_data(filtered_pathways)
    wellplates = {}
    start = ()

def find_path(neighbors_function, start, goal_test):
    if goal_test(start):
        return (start,)

    agenda = [(start,)]
    visited = {start}

    while agenda:
        this_path = agenda.pop(0)
        terminal_state = this_path[-1]

        for neighbor in neighbors_function(terminal_state):
            if neighbor not in visited:
                new_path = this_path + (neighbor,)
                if goal_test(neighbor):
                    return new_path

                agenda.append(new_path)
                visited.add(neighbor)

    return None
