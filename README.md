# Bottom-up hierarchical clustering

Implement bottom-up hierarchical clustering with the function cluster_bottom_up below. The function takes as input a list of profiles and a list of corresponding names (which will be the names for the leaves of the resulting tree). In addition, the function will take as input a string specifying which linkage function (e.g., single) to use for the clustering as well as a distance function that computes the distance (e.g., euclidean distance) between a pair of profiles. The function will output a TreeNode object, representing the root of the hierarchical clustering tree.

Distances: The tree should have branch lengths computed in the same way as for the UPGMA algorithm for phylogenetic trees. That is, each node should have a "height", with the leaf nodes being at height zero and the root node being the highest node. After merging two nodes (clusters),  𝑖  and  𝑗 , into a new node  𝑘 , the height of node  𝑘  should be half the distance from cluster  𝑖  to cluster  𝑗 , i.e.,  ℎ𝑒𝑖𝑔ℎ𝑡(𝑘)=𝑑𝑖𝑗/2 .

Tie-breaking: For tie-breaking purposes, you should keep track of an index for each node (cluster) in the tree. The input profiles will correspond to the leaves of the tree, and should have indices 0 to  𝑛−1  where  𝑛  is the number of profiles. Each successive node that is created should have the next available integer index (e.g., the very first merge of the algorithm should produce the node with index  𝑛  and the following merge should produce the node with index  𝑛+1 ). When finding the next pair of clusters to merge, if two or more pairs have the same minimum distance, pick the pair with the lexicographically smallest pair of indices  (𝑖,𝑗) . For example, if the pairs of clusters (3, 8) and (5, 7) have the same minimum distance, you should choose the pair (3, 8) to merge next.

Efficiency: Your implementation should have runtime complexity of  𝑂(𝑛3) . You are welcome to implement the more efficient  𝑂(𝑛2)  (for single-link) and  𝑂(𝑛2log𝑛)  (for complete and average-link) algorithms, but this is not required.

Hierarchical clustering data structure: You are to use objects of the TreeNode class as we did in notebook 22 to build your hierarchical clustering structure.

Note: To run the program that should install the "toytree" module first.