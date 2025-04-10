### Matching Augmentation Problem

The matching augmentation problem is defined as follows.

Given : A graph G=(V, E) and a matching M of G

Goal  : Find a set of edges F such that H=(V, F) is 2-edge-connected and |F\M| is minimized

This code is a C++ implementation of the primal LP-round scheme described in:

>Bamas, É., Drygala, M., Svensson, O. (2022). A Simple LP-Based Approximation Algorithm for the Matching Augmentation Problem. In: Aardal, K., Sanità, L. (eds) Integer Programming and Combinatorial Optimization. IPCO 2022. Lecture Notes in Computer Science, vol 13265. Springer, Cham. https://doi.org/10.1007/978-3-031-06901-7_5


## The implementation

The implementation is efficient (and, to the best of our knowledge, the fastest one) and parallelized using Open-MP.
The code depends on standard GCC libraries and the following external libraries:

Dependencies
* [Gurobi](www.gurobi.com)
* [Lemon](https://lemon.cs.elte.hu/trac/lemon)
* [Nauty](https://pallini.di.uniroma1.it/)

## Additional tools
There are two main sets of tools in this repository.

Int he ```BDS```folder, is the implementation of the BDS algorithm and a framework to test the integrality ratio of the 2-edge-conencted spanning subgraph polytope restricted to the Matching augmentation problem. 

In the ```conv_comb``` folder is an implementation of an algorithm to test the decomposition of half-integral solutions into 2-edge-connected spanning subgraphs.

Each implementation contains separate documentation of its readme file.

## Results obtained
We give a summary of the results derived from the code. The complete results are contained in the following:

>Azevedo, G. M. (2024). On rounding algorithms for the 2-edge-connected spanning subgraph problem. Master's Dissertation, Instituto de Matemática e Estatística, University of São Paulo, São Paulo. doi:10.11606/D.45.2024.tde-22012025-113608. Retrieved 2025-04-10, from www.teses.usp.br

We briefly summarize the computation results.

We present counterexamples to two conjectures posed by Bamas, Drygala, and Svensson, regarding the performance of the primal rounding algorithm in the half-integral case. Notably, we posed the first family of 'hard to approximate' instances for the problem, which is connected to non-planar graphs.

Moreover, we proposed the first infinity family of instances attaining the integrality gap (some instances were known before, but we were able to find new ones using exploratory research with nauty). Given the efficiency of the code, we tested a core family of graphs (quartic 4-edge-connected graphs) up to 18 vertices, even with the combinatorial explosion.

Finally, tested the conjectures posed in the thesis (see Chapters 5.5 and 3.4) in a vast set of instances. 

* Every 4-regular 4-edge-connected multigraph up to 15 vertices;
* 4-regular 4-edge-connected graphs from House of Graphs;
* 4-regular 4-edge-connected graphs obtained from snarks up to 34 vertices

