import toytree  # for working with trees
from toytree.TreeNode import TreeNode  # make TreeNode directly available
import math


def euclidean_distance(p1, p2):
    """The Euclidean distance between two profiles."""
    return math.sqrt(sum((e1 - e2) ** 2 for e1, e2 in zip(p1, p2)))


def manhattan_distance(p1, p2):
    """The Manhattan distance between two profiles."""
    return sum(abs(e1 - e2) for e1, e2 in zip(p1, p2))


def cluster_bottom_up(profiles, profile_names, linkage="single", distance=euclidean_distance):
    """Performs a bottom-down hierarchical clustering of a list of profiles, returning
    a tree that has the given profile names labeling the leaves.
    Args:
        profiles: a list of profiles/points (each of which is represented as a tuple)
        profile_names: a list of the same length as profiles giving the names of the profiles
        linkage: a string indicating the linkage method to use, which should be one of
            "single", "complete", and "average"
        distance: a function that takes as input two profiles and returns a number giving
                  the distance between the two profile, i.e., distance(p1, p2) should return
                  the distance between profile p1 and profile p2.
    Returns:
        A TreeNode instance representing the root of the hierarchical clustering tree.
    """

    nodes = {}
    coordinates = {}
    dists = {}

    active_nodes = []

    for i in range(len(profiles)):
        new_node = TreeNode(name=profile_names[i])
        nodes[i] = new_node
        coordinates[i] = [profiles[i]]
        active_nodes.append(i)

    def get_distance(coorA, coorB, linkage, distance):
        max_dist = 0
        min_dist = 99999999
        avg_dist = 0
        for pointA in coorA:
            for pointB in coorB:
                dist = distance(pointA, pointB)
                avg_dist += dist
                if dist < min_dist:
                    min_dist = dist
                if dist > max_dist:
                    max_dist = dist
        avg_dist /= len(coorA) * len(coorB)
        if linkage == 'single':
            return min_dist
        if linkage == 'complete':
            return max_dist
        if linkage == 'average':
            return avg_dist

    j = len(profiles) - 1
    while len(active_nodes) > 1:
        j += 1
        min_dist = 99999999
        for left in range(len(active_nodes) - 1):
            for right in range(left + 1, len(active_nodes)):
                nodeA = active_nodes[left]
                nodeB = active_nodes[right]
                coorA = coordinates[nodeA]
                coorB = coordinates[nodeB]

                if (nodeA, nodeB) in dists:
                    dist = dists[(nodeA, nodeB)]
                else:
                    dist = get_distance(coorA, coorB, linkage, distance)
                    dists[(nodeA, nodeB)] = dist
                    dists[(nodeB, nodeA)] = dist

                if dist < min_dist:
                    mergeA = nodeA
                    mergeB = nodeB
                    min_dist = dist

        new_node = TreeNode()
        new_node.add_child(nodes[mergeA])
        new_node.add_child(nodes[mergeB])

        coordinates[j] = coordinates[mergeA] + coordinates[mergeB]
        # new_node.height(min_dist / 2)
        if nodes[mergeA].is_leaf():
            nodes[mergeA].add_feature("dist", min_dist / 2)
        else:
            nodes[mergeA].add_feature("dist", min_dist / 2 - nodes[mergeA].dist)
        if nodes[mergeB].is_leaf():
            nodes[mergeB].add_feature("dist", min_dist / 2)
        else:
            nodes[mergeB].add_feature("dist", min_dist / 2 - nodes[mergeB].dist)
        new_node.add_feature("dist", min_dist / 2)

        nodes[j] = new_node

        active_nodes.remove(mergeA)
        active_nodes.remove(mergeB)
        active_nodes.append(j)

    '''
    newnode = TreeNode(name = 'A')
    root.add_child(newnode)
    newnode = TreeNode(name = 'B')
    root.add_child(newnode)
    '''
    root = nodes[active_nodes[0]]

    # root.sort_descendants()
    # print(root.write(format=1))
    # print(profiles, profile_names)
    return root
