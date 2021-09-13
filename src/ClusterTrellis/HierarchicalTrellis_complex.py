
import itertools
import string
import numpy as np
import logging
from .utils import get_logger
logger = get_logger(level=logging.WARNING)


class HierarchicalTrellis:
    """ Class to build the trellis given a set of leaves (jet constituents)"""

    def __init__(self):
        self.root = None
        self.elements = None
        self.elements_2_node = {}

    def leaves(self):
        return [self.elements_2_node[frozenset(element)] for element in self.elements]

    def get_complement_node(self, node, elements=None):
        if elements is None:
            elements = self.elements

        """Subtract elements in node"""
        elem = elements - node.elements
        # print("Elements after subtracting special element = ",elem)
        return self.elements_2_node.get(elem)

    def get_nodes_containing_element(self, root, element):
        ## could make an iterator...
        e_node = self.elements_2_node[frozenset(element)]
        nodes = {root}
        def get_parents(node):
            if node == root:
                return
            if set(node.elements).issubset(set(root.elements)):
                nodes.add(node)
                for parent in node.parents:
                    get_parents(parent)
        get_parents(e_node)
        nodes = list(nodes)
        return nodes

    @property
    def map_tree_energy(self):
        if self.root is None or self.root.map_tree_energy is None:
            raise ValueError
        return self.root.map_tree_energy


    def create_from_leaves(self,
                           map_values,
                           model_params,
                           NodeClass = None,
                           alphabet = string.ascii_lowercase):
        """ Create trellis nodes starting from leaves, and assign children and parents"""

        self.elements= frozenset(alphabet[:len(map_values)])
        """Create leaf nodes and assign values. e.g. momentum and delta"""
        for letter, leaf in zip(self.elements, map_values):
            element = frozenset(letter)
            leaf_node = NodeClass(model_params, element,
                                  None,leaf)

            """ Leaves have logLH= 0"""
            leaf_node.map_tree_energy = 0+0j

            """Fill dict with label as key and node as value"""
            self.elements_2_node[element] = leaf_node

        """ Create inner nodes: nodes with > 1 elements.
        Make clusters of every possible number of leaves, 
        from leaves to root, where root is a cluster of all leaves"""
        for i in range(2, len(self.elements) + 1):

            """Make all possible combinations of i leaves"""
            for elem in itertools.combinations(self.elements, i):
                logger.debug(f"ELEMENTS {elem}")
                logger.debug(f"---"*10)
                elem = frozenset(elem)

                """ Assign all clusters with i-1 leaves as children (where i enumerates the leaves of the parent subtree) """
                children = [self.elements_2_node[frozenset(child_elems)] for child_elems in
                            itertools.combinations(elem, i-1) if child_elems]

                """Create a node for a given cluster of i leaves. 
                Note: elem here denotes each cluster of i leaves"""
                a_node = NodeClass(model_params,elem, children)

                """Append parent node the the list of parents of each children node 
                (we created nodes for each children in the previous iteration (i-1)"""
                [child.parents.append(a_node) for child in children]

                """Append current node to the dictionary"""
                self.elements_2_node[elem] = a_node

        """ The last node created is the root"""
        self.root = a_node
        logger.debug(f"elements_2_node = {self.elements_2_node}")
        logger.info(f"Finished creating trellis nodes")
        return



    def compute_map_tree(self):
        """Finf the MLE (MAP) tree, and compute its momentums, deltas, and energies."""

        """Start from the leaf nodes"""
        for leaf in self.leaves():
            """ Compute map tree where each leaf is the root"""
            logger.debug(f"Leaf == {leaf}")
            leaf.compute_map_tree_split(self)

        """Set map tree for all inner nodes in order from smallest to largest number of elements.
        Start from the leaves"""
        elements = list(self.elements)
        logger.debug(f"elements = {elements}")
        for i in range(2, len(elements) + 1):
            for elem in itertools.combinations(elements, i):
                logger.debug(f"element = {elem}")
                elem = frozenset(elem)
                """Find node for each trellis label"""
                node = self.elements_2_node[elem]

                """ Compute map tree where the current node is the root"""
                node.compute_map_tree_split(self)

        logger.info(f"====="*10)
        logger.debug(f"Assigned MAP values for {self, self.root.map_tree_energy, self.root.map_children, self.root.map_features}")
        logger.info(f"---------"*10)
        logger.info(f"Partition function (logLH) = {self.root.logZ}")
        logger.info(f"MLE (LH, LLH) = {np.exp(self.map_tree_energy), self.map_tree_energy}")
        logger.info(f"Tot trees # =  {self.root.N_trees}")


        return self.map_tree_energy, self.root.logZ, self.root.N_trees




    def sample_trees(self,root, Ntrees):
        """ Sample trees LH following the posterior distribution over likelihoods"""
        treedist =[]
        for _ in range(Ntrees):
            treeLLH = []

            root.sample_singleTreeLH(self,root,treeLLH)
            treedist.append(np.sum(treeLLH))

            logger.debug(f"treeLLH = {treeLLH}")

        logger.info(f"treedist = {treedist}")

        return treedist



    def sampleTree(self,root, Ntrees):
        """ Sample trees (and build the tree dict) following the posterior distribution over likelihoods"""

        jetList = []
        treedist = []
        for _ in range(Ntrees):
            tree = []
            content = []
            leaves = []
            treeLLH = []
            deltas = []

            root.sample_singleTree(self,
                                     root,
                                     -1,
                                     False,
                                     tree,
                                     content,
                                     leaves,
                                     treeLLH,
                                     deltas,
                                     )

            jet ={}
            jet["tree"] = np.asarray(tree).reshape(-1, 2)
            jet["content"]= np.asarray(content)
            jet["leaves"]=np.asarray(leaves)
            jet["logLH"] = treeLLH
            jet["deltas"] = np.asarray(deltas)

            jetList.append(jet)

            treedist.append(np.sum(treeLLH))
            logger.debug(f"treeLLH = {treeLLH}")

        return treedist, jetList




    def traverseMAPtree( self, root):


        tree = []
        content = []
        leaves = []
        deltas = []

        root.traverseMLEtree(
            -1,
            False,
            tree,
            content,
            leaves,
            deltas,
        )


        jet ={}
        jet["tree"] = np.asarray(tree).reshape(-1, 2)
        jet["content"]= np.asarray(content)
        jet["leaves"]=np.asarray(leaves)
        jet["sumLogLH"] = root.map_tree_energy
        jet["deltas"] = np.asarray(deltas)
        jet["root_id"] = 0
        jet["Lambda"] = self.elements_2_node[frozenset({'a'})].lam
        jet["LambdaRoot"] = root.lam
        jet["pt_cut"] = root.delta_min

        # likelihood.enrich_jet_logLH(jet, dij=False)
        jet["totLogLH"]= sum(jet['logLH'])

        return jet

