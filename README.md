#  Cluster Trellis for Exact Inference on Hierarchical Clustering in Particle Physics

### **Craig S. Greenberg\*, Sebastian Macaluso\*, Nicholas Monath, Ji-Ah Lee, Patrick  Flaherty, Kyle Cranmer, Andrew McGregor, Andrew McCallum.**
#### * Both authors contributed equally to this work.

Note that this is an early development version. 


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Introduction

Hierarchical clustering is a fundamental task often used to discover meaningful structures in data, such as phylogenetic trees, taxonomies of concepts, subtypes of cancer, and cascades of particle decays in particle physics. Typically approximate algorithms are used for inference due to the combinatorial number of possible hierarchical clusterings. In contrast to existing methods, we present novel dynamic-programming algorithms for **exact** inference in hierarchical clustering based on a novel trellis data structure, and we prove that we can exactly compute the partition function, maximum likelihood hierarchy, and marginal probabilities of sub-hierarchies and clusters. Our algorithms scale in time and space proportional to the powerset of N elements which is super-exponentially more efficient than explicitly considering each of the (2N-3)!! possible hierarchies. 

<!--This should be  a jpg file for the figure to be displayed-->
![Fig.1](plots/LatentStructure5.jpg)

##### Fig. 1: Schematic representation of a hierarchical clustering. H denotes the latent state and X the dataset.

## Publications

A more detailed description of this work and implementations to other areas such as genomics can be found in our papers [`Data Structures & Algorithms for Exact Inference in Hierarchical Clustering`](https://arxiv.org/abs/2002.11661) and [`Exact Inference on Hierarchical Clustering in Particle Physics and Cancer Genomics`](https://ml4physicalsciences.github.io/2020/files/NeurIPS_ML4PS_2020_70.pdf).


## Data Structure

<!--This should be  a jpg file for the figure to be displayed-->
![Fig.2](plots/hierarchical_trellis_diagram_v7.png)

##### Fig. 2: Computing the partition function for the dataset {a,b,c,d}. Left: exhaustive computation, consisting of the summation of (2x4-3)!! =15  energy equations. Right: computation using the trellis.  The sum for the partition function is over 2^{4-1} - 1 = 7 equations, each making use of a memoized partition function value. Colors indicate corresponding computations over siblings in the trellis.

## Implementation on Ginkgo Jets

Ginkgo is a toy model for jets physics that can be downloaded from [`Ginkgo`](https://github.com/SebastianMacaluso/ToyJetsShower). A detailed description of the model can be found in [`Ginkgo Notes`](https://www.overleaf.com/read/pmmcqhyfsctf). Also, details and examples on how to access the jet binary tree structure in [`binaryTreeStructure.ipynb`](src/binaryTreeStructure.ipynb).

### Partition function and maximum likelihood (MAP) hierarchy.

 The trellis implements a memoized dynamic program to compute the partition function and the MAP hierarchy. There are examples on how to run the algorithm and plot the results in [`Ginkgo_Trellis.ipynb`](src/Ginkgo_Trellis.ipynb).



<img src="https://github.com/SebastianMacaluso/HierarchicalTrellis-Ginkgo/blob/master/plots/ZvsMAPscatter.png" width="395" align="left"><img src="https://github.com/SebastianMacaluso/HierarchicalTrellis-Ginkgo/blob/master/plots/posteriorSample100000_5_6.png" width="395" align="right">

<pre>


</pre>
##### Fig. 3. Left: Scatter plot of the partition function Z  vs. the trellis MAP value for a dataset of Ginkgo jets, with up to 10 leaves (jet constituents). The color indicates the number of leaves of each hierarchical clustering. Right: Comparison of the posterior distribution for  a  jet with five leaves for sampling 10^5 hierarchies (black dots with small error bars) and expected posterior distribution (in green). The plot shows the discrete nature of the distribution. The log likelihood for the ground truth tree is a vertical dashed red line.

<pre>

</pre>

## Installation


1. Clone the repository
2. `make install` from main dir.




-------------------------------------------------------------------------
## Citations

If you use this package, please cite our paper

```
@article{Greenberg:2020heb,
    author = "Greenberg, Craig S. and Macaluso, Sebastian and Monath, Nicholas and Lee, Ji-Ah and Flaherty, Patrick and Cranmer, Kyle and McGregor, Andrew and McCallum, Andrew",
    title = "{Data Structures \textbackslash{}\& Algorithms for Exact Inference in Hierarchical Clustering}",
    eprint = "2002.11661",
    archivePrefix = "arXiv",
    primaryClass = "cs.DS",
    month = "2",
    year = "2020"
}
```



<img src="https://github.com/SebastianMacaluso/ToyJetsShower/blob/master/notes/plots/IRIS-HEP.png" width="180" align="left"> <img src="https://github.com/SebastianMacaluso/ToyJetsShower/blob/master/notes/plots/NYU.png" width="140" align="center"> 
<img src="https://github.com/SebastianMacaluso/HierarchicalTrellis-Ginkgo/blob/master/plots/IESL_logo.png" width="90" align="center">                   <img src="https://github.com/SebastianMacaluso/HierarchicalTrellis-Ginkgo/blob/master/plots/UMASS_logo.png" width="90" align="center">








