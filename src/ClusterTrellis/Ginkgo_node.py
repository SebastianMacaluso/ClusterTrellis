
import numpy as np
import logging
from .trellis_node import TrellisNode
from . import Ginkgo_likelihood as likelihood
from ClusterTrellis.utils import get_logger
logger = get_logger(level=logging.WARNING)


class ModelNode(TrellisNode):
    """Class to define the nodes of the trellis and the node splitting likelihood """

    def __init__(self,
                 model_params,
                 elements = None,
                 children = None,
                 map_features = None):
        TrellisNode.__init__(self, model_params, elements, children, map_features)


    def get_energy_of_split(self, a_node, b_node):
        """ Assumes delta_min and lam are the same for both a_node and b_node"""

        logger.debug(f"computing energy of split: {a_node, b_node}")

        split_llh = likelihood.split_logLH(a_node.map_features[0],
                                           a_node.map_features[1],
                                           b_node.map_features[0],
                                           b_node.map_features[1],
                                           self.model_params["delta_min"],
                                           self.model_params["lam"])
        logger.debug(f"split_llh = {split_llh}")

        return split_llh


    def compute_map_features(self, a_node, b_node):
        """Auxiliary method to get the parent momentum value.
         The tree momentum is indepemdent of the tree latent structure"""
        momentum = a_node.map_features[0] + b_node.map_features[0]
        logger.debug(f"computing momentum for {a_node, b_node, momentum}")

        """Auxiliary method to get delta for the parent node."""
        logger.debug(f"computing delta parent {a_node, a_node.map_features[0], b_node, b_node.map_features[0]}")
        pP = a_node.map_features[0] + b_node.map_features[0]

        """Parent invariant mass squared"""
        tp1 = pP[0] ** 2 - np.linalg.norm(pP[1::]) ** 2
        logger.debug(f"tp =  {tp1}")

        return [momentum,tp1]
