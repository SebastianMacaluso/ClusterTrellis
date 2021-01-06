import time
import numpy as np
import HierarchicalTrellis
import logging
from utils import get_logger
logger = get_logger(level=logging.WARNING)


def compare_map_gt_and_bs_trees(gt_tree, NodeClass):
    """ Create a trellis given a set of leaves.
    (Originally written for Ginkgo: Toy Model for Particle Physics Jets)
    Args: Ground truth dictionary with leaves"""

    startTime = time.time()

    data_params = gt_tree
    N=len(data_params['leaves'])

    leaves_map_values =[ [data_params['leaves'][i],0] for i in range(N)]
    model_params ={}
    model_params["delta_min"] = float(data_params['pt_cut'])
    model_params["lam"]= float(data_params['Lambda'])

    """ Create trellis"""
    a_trellis = HierarchicalTrellis.HierarchicalTrellis()
    a_trellis.create_from_leaves(leaves_map_values, model_params, NodeClass = NodeClass)

    """Compute MAP (MLE), partition function Z"""
    map_energy, Z, Ntrees = np.float64(a_trellis.compute_map_tree())

    gt_energy = np.sum(gt_tree["logLH"])

    logger.debug({map_energy, gt_energy})

    endTime = time.time() - startTime

    return a_trellis, Z, map_energy, Ntrees, endTime




if __name__ == '__main__':
    main()