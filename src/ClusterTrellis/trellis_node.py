
import numpy as np
import logging
from scipy.special import logsumexp, softmax

from .utils import get_logger
logger = get_logger(level=logging.WARNING)


class TrellisNode:
    """Class to define the nodes of the trellis and the node splitting likelihood """

    def __init__(self,
                 model_params,
                 elements = None,
                 children = None,
                 map_features = None):

        self.model_params = model_params
        self.map_features = map_features
        self.sample_map_features = None

        if elements is None:
            elements = []
        if children is None:
            children = []
        self.elements = elements
        self.parents = []
        self.children = children
        self.map_children = []
        self.map_tree_energy = None
        self.Z = None
        self.logZ = None
        self.N_trees = 1
        self.levelLHweight = []

    def __repr__(self):
        return "HierarchicalTrellisNode: " + str(self.elements)

    def is_leaf(self):
        return not self.children


    def get_energy_of_split(self, a_node, b_node):
        """ Define the energy function"""
        raise NotImplementedError


    def compute_map_features(self, a_node, b_node):
        """ Method to compute model features for inner vertices"""
        raise NotImplementedError


    def compute_map_tree_split(self, a_trellis):
        """Compute the momentum, delta, and energies over the trellis.
        Save the MLE(MAP) energy split, marginal (partition function Z) and count number of trees.
        Assumes the map_energy of all the descendent nodes has already been computed!!"""
        logger.debug(f"Computing MAP tree rooted at {self}")


        """If the node is a leaf, set it to some default likelihood =1"""
        if self.is_leaf():
            self.map_tree_energy = 0
            self.logZ = 0.

            logger.debug(f"Leaf Z = {self.Z}")

        elif self.map_tree_energy is None and self.logZ is None:
            """Choose the 1st element of the parent node as a special element. 
            Get all sub-clusters that contain this element and their complement. 
            As this element has to be in one of the 2 sub-clusters that are the children, 
            then this way we cover all possible clusterings.
            For each node, it's "energy" is the llh of the join between the node and its complement.
            Find the max energy of all possible pairs of children for this parent node."""
            special_element = list(self.elements)[0]
            logger.debug(f"special_element = {special_element}")

            """ Get nodes containing element within current subtree"""
            special_nodes = a_trellis.get_nodes_containing_element(self, special_element)
            logger.debug(f"special nodes = {special_nodes}")

            """Root node can't be MAP (its the parent node => can't be a children node), so remove self from special_nodes"""
            special_nodes.remove(self)
            self.map_tree_energy = -np.inf
            self.logZ = -np.inf

            N_trees = 0

            for node in special_nodes:
                complement = a_trellis.get_complement_node(node, self.elements)
                logger.debug(f" complement = {complement}")

                """Compute the pairing of nodes energy (llh)"""
                split_llh = self.get_energy_of_split(node, complement)

                """Add last splitting energy to subtree energy."""
                energy = split_llh + node.map_tree_energy + complement.map_tree_energy

                """Count number of trees, only if tree is allowed under the model"""
                if energy > -np.inf:
                    N_trees += node.N_trees * complement.N_trees

                """Compute partition function in a stable way
                logZ_i = log[LH(a, b)] + logZ_a + logZ_b
                logZ_p = scipy.misc.logsumexp(np.asarray([logZ_p, logZ_i]))."""

                partial_logZ = split_llh + node.logZ + complement.logZ
                self.logZ = logsumexp(np.asarray([self.logZ, partial_logZ]))

                self.levelLHweight.append(partial_logZ)

                logger.debug(f" Pair of nodes Z = { node.Z, complement.Z}")
                logger.debug(f" Pair of nodes Energy = { node.map_tree_energy, complement.map_tree_energy}")


                """Save if new value for energy is higher"""
                if energy > self.map_tree_energy:

                    self.map_features = self.compute_map_features(node, complement)
                    self.map_tree_energy = energy

                    """Save the 2 clusters of leaves that give the MAP (MLE). This is for each inner node. 
                    The we could start from the root of the tree and traverse down the tree."""
                    self.map_children = (node, complement)

            """Assign momentum and delta for nodes that are not possible under the model, i.e. give a deltaP<delta_min. 
            These nodes have a llh=-np.inf so they don't contribute to the MLE or Z, but we need them to have a complete trellis"""
            if self.map_tree_energy == - np.inf:
                complement = a_trellis.get_complement_node(special_nodes[0], self.elements)
                self.map_features = self.compute_map_features(special_nodes[0], complement)


            self.N_trees = N_trees

            """Normalize"""
            logger.debug(f"level LH soft before = {self.levelLHweight}")

            if self.logZ == -np.inf:
                """Assign 0 probability when Delta_parent is below the treshold not allowed under the model"""
                self.levelLHweight = np.zeros(len(self.levelLHweight))

            else:
                """ Normalize to 1"""
                self.levelLHweight = softmax(np.asarray(self.levelLHweight))

            logger.debug(f"---------" * 10)
            logger.debug(f"Special nodes= {special_nodes}")
            logger.debug(f"level LH soft after = {self.levelLHweight}")
            logger.debug("---------" * 10)
            logger.debug("---------" * 10)



    def sample_singleTreeLH(self, a_trellis, root, treeLLH):
        """ After we save the llh (energy) for each node in the trellis, we can sample trees following the posterior distribution for the likelihood of each of them. This follows a top down approach and can start from any subtree (of the root of the tree)"""

        if self.is_leaf():
            return

        """Choose the 1st element of the parent node as a special element.
        Get all sub-clusters that contain this element and their complement.
        As this element has to be in one of the 2 sub-clusters that are the children,
        then this way we cover all possible clusterings.
        For each node, it's "energy" is the llh of the join between the node and its complement.
        Find the max energy of all possible pairs of children for this parent node."""
        special_element = list(root.elements)[0]
        logger.debug(f"special_element = {special_element}")

        """ Get nodes containing element within current subtree"""
        special_nodes = a_trellis.get_nodes_containing_element(root, special_element)


        """Root node can't be MAP (its the parent node => can't be a children node), so remove self from special_nodes"""
        special_nodes.remove(root)

        logger.debug(f"root  = {root}")
        logger.debug(f"special nodes = {special_nodes}")

        logger.debug(f"root.levelLHweight = {root.levelLHweight}")


        """Sample a node at current level following the likelihhood of each of them"""
        node = np.random.choice(special_nodes, 1, p=root.levelLHweight)[0]

        complement = a_trellis.get_complement_node(node, root.elements)

        logger.debug(f"node = {node}")
        logger.debug(f"complement = {complement}")

        """ Get llh for the join of the pair {node,comlement} sampled"""
        split_llh = self.get_energy_of_split(node, complement)

        treeLLH.append(split_llh)

        self.sample_map_features = self.compute_map_features(node, complement)


        """ Recursively repeat for the next level"""
        node.sample_singleTreeLH(a_trellis,
                               node,
                               treeLLH)

        complement.sample_singleTreeLH(a_trellis,
                                     complement,
                                     treeLLH)



    def traverseMLEtree(self, parent_idx, is_left, tree, content, leaves, deltas):

        idx = len(tree) // 2
        if parent_idx >= 0:
            if is_left:
                tree[2 * parent_idx] = idx
            else:
                tree[2 * parent_idx + 1] = idx

        tree.append(-1)
        tree.append(-1)


        content.append(self.map_features[0])

        logger.debug(f"Tree = {tree}")
        logger.debug(f"Content = {content}")
        logger.debug(f"self.map_children = {self.map_children}")
        logger.debug(f"len(self.map_children)= {len(self.map_children)}")

        if len(self.map_children)<2:
            leaves.append(self.map_features[0])
            deltas.append(0.)

        else:

            deltas.append(self.map_features[1])

            self.map_children[0].traverseMLEtree(
                idx,
                True,
                tree,
                content,
                leaves,
                deltas,
            )

            self.map_children[1].traverseMLEtree(
                idx,
                False,
                tree,
                content,
                leaves,
                deltas,
            )



    def sample_singleTree(self, a_trellis, root,  parent_idx, is_left, tree, content, leaves, treeLLH, deltas):
        """ After we save the llh (energy) for each node in the trellis, we can sample trees following the posterior distribution for the likelihood of each of them. This follows a top down approach and can start from any subtree (of the root of the tree)"""


        """Build tree dictionary with content, tree, leaves, etc"""
        idx = len(tree) // 2
        if parent_idx >= 0:
            if is_left:
                tree[2 * parent_idx] = idx
            else:
                tree[2 * parent_idx + 1] = idx

        tree.append(-1)
        tree.append(-1)

        logger.debug(f"Tree = {tree}")
        logger.debug(f"Content = {content}")

        if self.is_leaf():
            leaves.append(self.map_features[0])
            deltas.append(0.)
            content.append(self.map_features[0])

            return

        """Choose the 1st element of the parent node as a special element. 
        Get all sub-clusters that contain this element and their complement. 
        As this element has to be in one of the 2 sub-clusters that are the children, 
        then this way we cover all possible clusterings.
        For each node, it's "energy" is the llh of the join between the node and its complement.
        Find the max energy of all possible pairs of children for this parent node."""
        special_element = list(root.elements)[0]
        logger.debug(f"special_element = {special_element}")

        """ Get nodes containing element within current subtree"""
        special_nodes = a_trellis.get_nodes_containing_element(root, special_element)


        """Root node can't be MAP (its the parent node => can't be a children node), so remove self from special_nodes"""
        special_nodes.remove(root)

        logger.debug(f"root  = {root}")
        logger.debug(f"special nodes = {special_nodes}")

        logger.debug(f"root.levelLHweight = {root.levelLHweight}")


        """Sample a node at current level following the likelihhood of each of them"""
        node = np.random.choice(special_nodes, 1, p=root.levelLHweight)[0]

        complement = a_trellis.get_complement_node(node, root.elements)

        logger.debug(f"node = {node}")
        logger.debug(f"complement = {complement}")

        """ Get llh for the join of the pair {node,comlement} sampled"""
        split_llh, = self.get_energy_of_split(node, complement)

        temp_map_features = self.compute_map_features(node, complement)
        content.append(temp_map_features[0])
        deltas.append(temp_map_features[1])
        treeLLH.append(split_llh)


        """ Recursively repeat for the next level"""
        node.sample_singleTree(a_trellis,
                               node,
                               idx,
                               True,
                               tree,
                               content,
                               leaves,
                               treeLLH,
                               deltas,
                               )

        complement.sample_singleTree(a_trellis,
                                     complement,
                                     idx,
                                     False,
                                     tree,
                                     content,
                                     leaves,
                                     treeLLH,
                                     deltas,
                                     )