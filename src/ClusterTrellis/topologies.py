import itertools


import pickle
import os
import logging
import numpy as np
import time
# import wandb
# import sys
#
#
# from absl import flags
# from absl import logging
# from absl import app

from ginkgo import likelihood_invM as likelihood
from scipy.special import logsumexp

from absl import logging
logging.set_verbosity(logging.WARNING)


class build_trees_topologies(object):

    def __init__(self, max_nodes, propagate_values_up, Ntopologies):
        pClass_st = time.time()
        # children[i] gives the kids of node i
        self.nodes_explored = 0
        self.topologies = [[] for _ in range(Ntopologies)]
        self.Ntrees = 0
        self.children = [[] for _ in range(max_nodes)]
        self.children_elems = [[] for _ in range(max_nodes)]
        self.mapv = np.zeros(max_nodes, dtype=np.float32)  # MAP vertices
        self.arg_mapv = [[] for _ in range(max_nodes)]  # Children of mapv
        self.clusters = [[] for _ in range(max_nodes)]  # list of frozensets
        # pq[i] = f, g, h, a, b
        self.pq = [[] for _ in range(max_nodes)]  # List of priority queues where each entry index is the node id.
        self.explored = np.zeros(max_nodes, dtype=np.bool_)
        self.elems2node = dict()  # Nodes are indices and not a class like in the standard full trellis
        self.root = None
        self.next_id = 0
        self.last_memoized = -1 * np.ones(max_nodes, dtype=np.float32)
        self.current_extension_num = -1
        self.propagate_values_up = propagate_values_up
        if self.propagate_values_up:
            # parents[i] gives parents of node i
            self.parents = [[] for _ in range(max_nodes)]
            self.up_to_date_fgh = [{} for _ in range(max_nodes)]  # [i] -> {(c_l, c_r) -> f,g,h}
            self.values_on_q = [{} for _ in range(max_nodes)]  # [i] -> {(c_l, c_r) -> queue([f_j,g_j,h_j,...])}

    def _find_Ntopologies(self,elements):
        # if len(elements)==1: # we consider that leaves have 0 trees below
        #     return 0
        # else:
        Tot = 1
        for k in range(2 * len(elements) - 3, 1, -2):
            Tot *= k
        return Tot

    def set_root(self, elems):
        self.root = self.record_node(elems)
        self.clusters[self.root] = elems
        self.nodes_explored += int(self.root)  # Add the number of elements to the nodes explored

    def is_leaf(self, i):
        """ Returns true if i is a leaf.

        True if i has no children

        :param i: node id
        :return: if there are no children.
        """
        return len(self.children[i]) == 0 or len(self.pq[i]) == 0

    def get_children(self, i, elems):
        """Gives the children of node i that has elements elems.

        In this version, it grabs all 2 partitions if they are not there
        and caches this in children[i].

        :param i: The node i
        :param elems: The elements associated with node i
        :return: the children of i
        """
        # if len(elems) == 1:
        #     return []
        # elif self.explored[i]:
        #     return self.children[i]
        # else:
        self.children[i], self.children_elems[i] = self._get_children(list(elems))  # all_two_partitions(list(elems))

            # self.update_from_children(i, (ch_l, ch_r))
        return self.children[i], self.children_elems[i]

    def _get_children(self, root):
        ##
        nodes = []
        nodes_elems = []
        root = frozenset(root)
        if len(root)>1:
            other_element = frozenset(list(root)[1::])


            for i in range(1, len(other_element) + 1):
                nodes += [(self.record_node(frozenset(elems)), self.record_node(root - frozenset(elems))) for elems in
                          itertools.combinations(other_element, i)]
                nodes_elems += [(frozenset(elems), root - frozenset(elems)) for elems in
                          itertools.combinations(other_element, i)]
        # else:
        #     nodes+=[self.elems2node[frozenset(root)]]
        #     nodes_elems+=[root]

        return nodes, nodes_elems

    def build_topologies(self, elements, all_nodes, m=0, uncle = None ):

        # logging.info("Elements = %s", elements)
        k= self.elems2node[elements]
        # k=self.root
        # elements = self.clusters[k]
        ch, ch_elems = self.get_children(k, elements)
        all_nodes.append([elements,ch_elems])
        # logging.info("ch = %s", ch_elems)
        # logging.info("m = %s", str(m))
        # logging.info("# trees = %s", str(self._find_Ntopologies(elements)))
        # if elements == self.clusters[self.root]:
        #     [self.topologies[i].append(elements) for i in range(self._find_Ntopologies(elements))]

        # if len(elements)==1:
        #     # logging.info("Elements = %s", elements)
        #     self.topologies[m].append(elements)

        # if len(ch)>0:
        for pairwise_split in ch:
            logging.debug("pairwise_split = %s ", pairwise_split)

            # # logging.info("pairwise_split_clusters = %s ", [(self.clusters[pairwise_split[0]], len(self.clusters[pairwise_split[0]])), (self.clusters[pairwise_split[1]], self.clusters[pairwise_split[1]])])
            # pairwise_split_clusters = sorted([(self.clusters[pairwise_split[0]], len(self.clusters[pairwise_split[0]])), (self.clusters[pairwise_split[1]], len(self.clusters[pairwise_split[1]]))], key=lambda x: x[1])
            # pairwise_split_clusters = [x for (x,y) in pairwise_split_clusters]
            #
            # logging.info("pairwise_split_clusters = %s ", pairwise_split_clusters)
            # self.build_topologies(pairwise_split_clusters[0], all_nodes, m)
            #
            # m+=self._find_Ntopologies(pairwise_split_clusters[0])
            # # if m==len(self.clusters[self.root]): m-=1
            # self.build_topologies(pairwise_split_clusters[1], all_nodes, m)
            # m += self._find_Ntopologies(pairwise_split_clusters[1])


            L = self.clusters[pairwise_split[0]]
            R = self.clusters[pairwise_split[1]]

            Nsubtrees = self._find_Ntopologies(L) * self._find_Ntopologies(R)

            N_uncle=1
            if uncle:
                N_uncle = self._find_Ntopologies(self.clusters[uncle])

            [self.topologies[m + i].append((L,R)) for i in range(Nsubtrees*N_uncle)]

            self.build_topologies(L, all_nodes, m, uncle=pairwise_split[1])
            self.build_topologies(R, all_nodes, m, uncle=pairwise_split[0])





            # L = self.clusters[pairwise_split[0]]
            # R = self.clusters[pairwise_split[1]]
            #
            # Nsubtrees = self._find_Ntopologies(L) * self._find_Ntopologies(R)
            # [self.topologies[m + i].append((L,R)) for i in range(Nsubtrees)]
            #
            # self.build_topologies(L, all_nodes, m, p=pairwise_split[0])

            # m+=self._find_Ntopologies(self.clusters[pairwise_split[0]])
            # if m==len(self.clusters[self.root]): m-=1
            # self.build_topologies(R, all_nodes, m)

            # Nsubtrees = self._find_Ntopologies(L) * self._find_Ntopologies(R)
            # [self.topologies[m + i].append((L,R)) for i in range(Nsubtrees)]
            m += Nsubtrees


        # logger.debug(f"element = {elem}")
        # elem = frozenset(elem)
        # """Find node for each trellis label"""
        # node = self.elements_2_node[elem]
        #
        # """ Compute map tree where the current node is the root"""
        # node.compute_map_tree_split(self)


    def tree_log_likelihood(self):
        pass



    def invariant_permutations(self):
        pass



class build_jetTrees(build_trees_topologies):

    """We want to find a maximum likelihood, so to implement the priority queue with min heap, we flip all the likelihood values. Thus, we will find a solution for the minimum value of (- log LH). Note: we work with (- log LH) for practical purposes """

    def __init__(self, propagate_values_up: bool, max_nodes: int, Ntopologies: int, leaves = None,  min_invM=None, Lambda= None, LambdaRoot=None):
        super(build_jetTrees, self).__init__( propagate_values_up=propagate_values_up, max_nodes=max_nodes, Ntopologies=Ntopologies)
        # self.graph = graph
        st = time.time()
        self.leaves_momentum = leaves
        # print("self.leaves_momentum = ", self.leaves_momentum)
        self.momentum = dict()
        self.min_invM = min_invM
        self.Lambda = Lambda.numpy()
        self.LambdaRoot = LambdaRoot.numpy()
        # self.invariant_mass = dict()
        # assert self.graph.shape[0] == self.graph.shape[1], "Input graph must be square matrix!"
        # initialize all singleton clusters to be nodes 0 to num points -1
        # print("self.leaves_momentum.shape[0] = ",self.leaves_momentum.shape[0])
        for i in range(self.leaves_momentum.shape[0]):
            self.momentum[frozenset({i})] = self.leaves_momentum[i]
            self.record_node(frozenset({i}))
            # elif len(elements)==1:

        # print("self.momentum=", self.momentum)
        self.set_root(frozenset(range(self.leaves_momentum.shape[0])))
        logging.info('Root node has elements: %s', self.clusters[self.root])
        logging.info('Root node is id: %s', self.root)

        en_t = time.time()- st
        logging.info('Loading class time: %s', en_t)





    def record_node(self, elements: frozenset) -> int:
        """Get the node corresponding to the given elements, create new id if needed.

        Creates a new id if needed.

        :param elements: the elements (FrozenSet of Integers)
        :return: the node for the given elements (this is an index id)
        """
        logging.debug('get node id from elements %s', str(elements))
        if elements not in self.elems2node:
            logging.debug('get node id from elements %s. new node! %s', str(elements), self.next_id)
            logging.debug('Clusters =%s ', str(self.clusters))
            self.elems2node[elements] = self.next_id
            self.clusters[self.next_id] = elements
            if len(elements)>1:
                # print('element in elements=', [element for element in elements])
                # print("momentum =", np.asarray([self.momentum[frozenset({elem})] for elem in elements]))
                self.momentum[elements]= sum(np.asarray([self.momentum[frozenset({elem})] for elem in elements])) # Add the momentum of the leaves that compose the node
                # self.invariant_mass[self.next_id] =
            # elif len(elements)==1:
            #     self.momentum[elements]= self.leaves_momentum[list(elements)[0]]

            self.next_id += 1
            return self.next_id - 1
        else:
            return self.elems2node[elements]



    def invariant_permutations(self,tree,tree_with_all_permutations, i=0):


        # for tree in trees:
        # tree_with_all_permutations = [tree]
        temp = tree[0:i] + [(tree[i][1], tree[i][0])] + tree[i + 1::]
        if temp not in tree_with_all_permutations:
            tree_with_all_permutations.append(temp)
        for k,siblings in enumerate(temp):
            # if k>i:
            temp2=temp[0:k]+[(siblings[1], siblings[0])]+temp[k+1::]
            if temp2 not in tree_with_all_permutations:
                tree_with_all_permutations.append(temp[0:k]+[(siblings[1], siblings[0])]+temp[k+1::])
            # tree_logLH.append(self.get_energy_of_split(siblings[0], siblings[1]))

        if i<len(tree)-1:
            i += 1
            self.invariant_permutations(tree, tree_with_all_permutations, i)


    # def invariant_permutations(self,tree,tree_with_all_permutations, i=0):
    #
    #
    #     # for tree in trees:
    #     # tree_with_all_permutations = [tree]
    #     temp = tree[0:i] + [(tree[i][1], tree[i][0])] + tree[i + 1::]
    #     tree_with_all_permutations.append(temp)
    #     for k,siblings in enumerate(temp):
    #         tree_with_all_permutations.append(tree[0:k]+[(siblings[1], siblings[0])]+tree[k+1::])
    #             # tree_logLH.append(self.get_energy_of_split(siblings[0], siblings[1]))
    #
    #     if i<len(tree)-1:
    #         i += 1
    #         self.invariant_permutations(tree, tree_with_all_permutations, i)

        # for i,siblings in enumerate(tree):
        #     tree_with_all_permutations.append(tree[0:i]+[(siblings[1], siblings[0])]+tree[i+1::])
        #         # tree_logLH.append(self.get_energy_of_split(siblings[0], siblings[1]))
        #     # total_logLH.append(sum(tree_logLH))
        # return tree_with_all_permutations

    # def invariant_permutations(self,tree):
    #
    #     # for tree in trees:
    #     tree_with_all_permutations = [tree]
    #     for i,siblings in enumerate(tree):
    #         tree_with_all_permutations.append(tree[0:i]+[(siblings[1], siblings[0])]+tree[i+1::])
    #             # tree_logLH.append(self.get_energy_of_split(siblings[0], siblings[1]))
    #         # total_logLH.append(sum(tree_logLH))
    #     return tree_with_all_permutations

    def tree_log_likelihood(self, trees):
        total_logLH = []
        for tree in trees:
            tree_logLH = []
            for siblings in tree:
                tree_logLH.append(self.get_energy_of_split(siblings[0], siblings[1]))
            total_logLH.append(sum(tree_logLH))

        Z = logsumexp(np.asarray(total_logLH))
        return Z, max(total_logLH)

    def get_energy_of_split(self, elem_l, elem_r):
        """Add last splitting llh to subtree llh.
        Assumes delta_min and lam are the same for both a_node and b_node"""

        # logging.debug(f"computing energy of split: {ch_l, ch_r}")

        # elem_l = self.clusters[ch_l]
        # elem_r = self.clusters[ch_r]

        # To follow the convention on Ginkgo, to get the correct result, we set t==0 if we have a leaf, i.e. t<t_cut
        l_node_invM = 0
        r_node_invM =0
        l_node_invM = self.momentum[elem_l][0] ** 2 - np.linalg.norm(self.momentum[elem_l][1::]) ** 2
        r_node_invM = self.momentum[elem_r][0] ** 2 - np.linalg.norm(self.momentum[elem_r][1::]) ** 2

        logging.debug(f" t_l ={l_node_invM}")
        logging.debug(f" t_R ={r_node_invM}")

        logging.debug(f" p_l ={self.momentum[elem_l]}")
        logging.debug(f" p_r ={self.momentum[elem_r]}")

        split_llh = likelihood.split_logLH_with_stop_nonstop_prob(self.momentum[elem_l],
                                                       self.momentum[elem_r],
                                                       self.min_invM,
                                                       self.Lambda)


        # logging.debug(f"split_llh = {split_llh}")

        # llh = split_llh + a_node.map_tree_energy + b_node.map_tree_energy
        logging.debug(f"split likelihood ={split_llh}")

        return split_llh


