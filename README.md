# Self Similarity in Human Brain Connectomes

This repository contains Fortran programs developed for my Bachelor's Thesis in Physics at the University of Barcelona, focusing on Complex Networks. The programs explore the multiscale organization of human brain connectomes.

## Abstract

The human brain is a highly intricate system, characterized by interactions across multiple spatial scales. This project investigates the multiscale organization of human brain connectomes using two datasets from healthy subjects. Two symmetries were identified through the analysis.

1. **Multiscale Properties:** The study examined multiscale properties at five hierarchical resolutions, reaffirming the self-similarity of these properties as the resolution decreases.
   ![Multiscale Properties](multiscale_properties.pdf)

2. **Structural Features:** A quantitative analysis demonstrated that structural features of connectomes remain self-similar when applying a degree-thresholding renormalization method on the highest resolution network layer.
   ![Structural Features](structural_features.pdf)

## Fortran programs

The Fortran programs in this repository analyze the structural properties of graphs, specifically concentrating on:

- `clustering.f90`: Calculates the local clustering coefficient for each vertex in the graph.

- `edgelistHCP.f90`: Converts an edgelist to an adjacency list for a dataset with no repeated links.

- `edgelistUL.f90`: Converts an edgelist to an adjacency list for a dataset with repeated links.

- `knn.f90`: Determines the average nearest neighbor degree as a measure of the average connectivity of vertices in the graph.

- `probability.f90`: Calculates the degree distribution of the graph.

- `rewiring.f90`: Applies a random rewiring method to randomize the network structure.

- `rnd_threshold.f90`: Applies threshold renormalization to the randomized networks.

- `thresholdHCP.f90`: Applies a threshold renormalization to an edgelist with no repeated links.

- `thresholdUL.f90`: Applies a threshold renormalization to an edgelist.


## Note on Code Organization
It's important to note that the Fortran code may not be thoroughly organized. 
