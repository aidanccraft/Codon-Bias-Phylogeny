import numpy as np


class Node:
    def __init__(self):
        self.left = None
        self.right = None
        self.val = None


def differences(seq1, seq2):
    '''Counts the number of pairwise differences between
    two sequences'''
    count = 0

    # Count the number of pairwise mismatches
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            count += 1

    return count


def create_diff_dict(animals, nodes):
    '''Creates a matrix of differences between the given species.
    Stored in a dictionary of dictionary based on the indices of 
    the provided list.'''
    diff_dict = {}

    # Find the difference between each pair of animals
    for i, animal1 in enumerate(nodes):
        for j, animal2 in enumerate(nodes):
            # Don't check animals against themselves
            if(animal1 == animal2):
                continue

            diff = differences(animals[i], animals[j])

            # Update difference matrix with distance calculated
            if(not animal1 in diff_dict):
                diff_dict[animal1] = {}

            if(not animal2 in diff_dict):
                diff_dict[animal2] = {}

            diff_dict[animal1][animal2] = diff
            diff_dict[animal2][animal1] = diff

    return diff_dict


def pair_group(diff_dict):
    """ given a matrix of differences, returns the indices of the closest two related organisms"""
    min_val = np.inf
    index_1 = 0
    index_2 = 0

    # Find the lowest number of pairwise different
    for i in diff_dict.keys():
        for j in diff_dict[i].keys():
            if diff_dict[i][j] < min_val:
                min_val = diff_dict[i][j]
                index_1 = i
                index_2 = j

    # Return the pair of indices that have the lowest score
    return (index_1, index_2)


def calculate_weight(tup):
    '''Determine number of items in a nested tuple'''
    if(type(tup) != tuple):
        return 1

    num_values = 0

    # Recursively determine weight of the tuple by unpacking it based on its contents
    for i in tup:
        if type(i) == tuple:
            num_values += calculate_weight(i)
        else:
            num_values += 1

    return num_values


def get_distance_to_farthest_node(root):
    '''Gets the linear distance from a node to its furthest child'''
    if(root.right == None):
        return 0

    # Return the distance between this node and its child
    return root.right[1] + get_distance_to_farthest_node(root.right[0])


def generate_phylogenetic_tree(diff_dict):
    '''Creates a phylogentic tree from a list of sequences'''
    # Create the dictionary of nodes for the phylogentic tree creation
    node_dict = {}

    for species in diff_dict:
        node = Node()
        node.val = species
        node_dict[species] = node

    # Apply UPGMA until all nodes have been combined
    while(len(diff_dict.keys()) > 1):
        # Get the two most related nodes
        (animal_1, animal_2) = pair_group(diff_dict)

        # Add a new row for the new node
        new_group = (animal_1, animal_2)
        diff_dict[new_group] = {}

        # Get the distance between these two nodes
        distance = diff_dict[animal_1][animal_2]

        # Create a new Node representing the new group being formed. Include the edge weight between the nodes
        node = Node()
        node.val = new_group
        node.left = (node_dict[animal_1], distance / 2 -
                     get_distance_to_farthest_node(node_dict[animal_1]))
        node.right = (node_dict[animal_2], distance / 2 -
                      get_distance_to_farthest_node(node_dict[animal_2]))

        # Remove the old nodes and add the new ones to the node dictionary
        node_dict.pop(animal_1)
        node_dict.pop(animal_2)
        node_dict[new_group] = node

        # Delete node pairs from the distance matrix
        diff_dict[animal_1].pop(animal_2)
        diff_dict[animal_2].pop(animal_1)

        # Calculate the weights for the arithmetic mean
        weight1 = calculate_weight(animal_1)
        weight2 = calculate_weight(animal_2)

        # Populate new node's row
        for i in diff_dict[animal_1].keys():
            diff_dict[new_group][i] = (
                weight1 * diff_dict[animal_1][i] + weight2 * diff_dict[animal_2][i]) / (weight1 + weight2)

        # Delete old node's rows
        diff_dict.pop(animal_1)
        diff_dict.pop(animal_2)

        # Update remaining nodes' distances to the newly added node and delete the remaining
        # references to the old nodes
        for i in diff_dict.keys():
            if(i != new_group):
                diff_dict[i][new_group] = (
                    weight1 * diff_dict[i][animal_1] + weight2 * diff_dict[i][animal_2]) / (weight1 + weight2)

                diff_dict[i].pop(animal_1)
                diff_dict[i].pop(animal_2)

    # The last created node is the completed tree
    return node


def write_tree_xml(root, filename):
    '''Writes a weighted binary tree to a PhyloXML file'''
    with open('../data/' + filename + '.xml', 'w') as f:
        # Write the necessary header information
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns="http://www.phyloxml.org">\n')
        f.write('<phylogeny rooted="true">\n')
        f.write('<clade>')

        # Recursively write each node to the file
        write_tree_recursive(root, f)

        # Close off the file
        f.write('</clade>\n')
        f.write('</phylogeny>\n')
        f.write('</phyloxml>\n')


def write_tree_recursive(root, f):
    ''' Recursively adds clades to the file'''
    if(root.right == None):
        f.write('<name>' + root.val + '</name>\n')
        return

    # Add a new clade for the left child
    f.write('<clade branch_length="' + str(root.left[1]) + '">\n')
    write_tree_recursive(root.left[0], f)
    f.write('</clade>\n')

    # Add a new clade for the right child
    f.write('<clade branch_length="' + str(root.right[1]) + '">\n')
    write_tree_recursive(root.right[0], f)
    f.write('</clade>\n')
