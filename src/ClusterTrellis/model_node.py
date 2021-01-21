
import logging
from .trellis_node import TrellisNode
from .utils import get_logger
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
        """Model energy function
        Args: sibling nodes
        returns: model splitting energy
        """
        logger.debug(f"computing energy of split: {a_node, b_node}")

        """"Add function to compute model splitting energy"""

        """ split_energy = ...."""
        pass



    def compute_map_features(self, a_node, b_node):
        """Auxiliary method to get the parent vertex features. This is model dependent. In Ginkgo these are the momentum and parent invariant mass.
        Args: sibling nodes
        returns: list where each entry is a parent feature, e.g. [feature1,feature2,...]
        """
        pass
