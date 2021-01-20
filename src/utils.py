import os
import pickle
import string
import time
import logging
import numpy as np



def get_logger(name=__file__, level=logging.INFO):
    logger = logging.getLogger(name)

    if getattr(logger, "_init_done__", None):
        logger.setLevel(level)
        return logger

    logger._init_done__ = True
    logger.propagate = False
    logger.setLevel(level)

    formatter = logging.Formatter("%(asctime)s:%(levelname)s::%(message)s")
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    handler.setLevel(0)

    del logger.handlers[:]
    logger.addHandler(handler)

    return logger


## Utils
def load_jets():
    root_dir = "data/"
    filename = os.path.join(root_dir, "TruthBS_10")
    with open(filename + ".pkl", "rb") as fd:
        Truth10, BS10 = pickle.load(fd, encoding='latin-1')
    return Truth10, BS10

def sumLogLH(jetList):
    for jet in jetList:
        jet["totLogLH"] = np.sum(jet["logLH"])

def getConstituents(jet, node_id, outers_list):
    """
    Recursive function to get a list of the tree leaves
    """
    if jet["tree"][node_id, 0] == -1:

        outers_list.append(jet["content"][node_id])

    else:
        getConstituents(
        jet,
        jet["tree"][node_id, 0],
        outers_list,)

        getConstituents(
        jet,
        jet["tree"][node_id, 1],
        outers_list,)

    return outers_list

def get_leaves(jet):
    return getConstituents(jet, jet["root_id"], [])