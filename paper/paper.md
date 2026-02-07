---
title: "RSNet: A Resampling-Based Framework for Network Structure Learning in High-Dimensional Data"

authors:
  - name: Ziwei Huang
    orcid: 0009-0005-4879-8175
    affiliation: [1,5]
  
  - name: Zeyuan Song
    orcid: 0000-0002-7352-4177
    affiliation: [2,3]
  
  - name: Paola Sebastiani
    orcid: 0000-0001-6419-1545
    affiliation: [2,3,4]
  
  - name: Stefano Monti
    orcid: 0000-0002-9376-0660
    affiliation: [5,6,7]

affiliations:
  - name: Department of Physics, Boston University, Boston, MA
    index: 1
  - name: Institute for Clinical Research and Health Policy Studies, Tufts Medical Center, Boston, MA
    index: 2
  - name: Department of Medicine, School of Medicine, Tufts University, Boston, MA
    index: 3
  - name: Data Intensive Study Center, Tufts University, Boston, MA
    index: 4
  - name: Division of Computational Biomedicine, Boston University Chobanian & Avedisian School of Medicine, Boston, MA
    index: 5
  - name: Department of Biostatistics, Boston University School of Public Health, Boston, MA
    index: 6
  - name: Bioinformatics Program, Faculty of Computing and Data Science, Boston University, Boston, MA
    index: 7


date: 2026-02-06
bibliography: references.bib
---

# Summary

RSNet is an open-source R package that provides a resampling-based framework for robust and interpretable network inference, designed to address the limited-sample-size challenges common in high-dimensional (e.g., ‘omics’) data. It supports both the estimation of partial correlation networks modeled as Gaussian networks [@whittaker2009graphical] and conditional Gaussian Bayesian networks for mixed data types that combine continuous and discrete variables [@lauritzen2001stable]. The framework incorporates multiple resampling strategies, including bootstrap, subsampling, and cluster-based approaches, to accommodate both independent and correlated (e.g., family-based) observations. To enhance interpretability, RSNet integrates graphlet-based topology analysis that captures higher-order connectivity and edge sign information, enabling single-node and subnetwork-level insights. Notably, RSNet is the first R package to efficiently construct signed graphlet degree vector matrices (GDVMs) in near-constant time for sparse networks, providing scalable analysis of higher-order network structure. Collectively, RSNet offers a versatile tool for statistically reliable and interpretable network inference in high-dimensional data.

# Statement of need

Network inference methods are widely used to model dependencies among variables in high-dimensional data, supporting discovery and hypothesis generation in diverse research domains [@federico2023structure; @fan2016overview; @friedman2008sparse; @friedman2000using]. Commonly applied approaches such as correlation or co-expression networks [@zhang2005general; @yang2014gene; @ruan2010general] are easy to implement but cannot distinguish direct from indirect dependencies [@federico2023structure]. In contrast, Gaussian networks, also known as partial correlation neworks [@federico2023structure; @fan2016overview; @friedman2008sparse], and conditional Gaussian Bayesian networks (CGBNs)[@mcgeachie2014cgbayesnets; @bottcher2005learning] estimate conditional dependencies, offering a higher-resolution representation of complex systems. However, the reliability of inferred network structures is often compromised by the limited sample sizes in high-dimensional data, a challenge commonly referred to as the “small n, large p” problem, where n denotes the number of samples and p the number of variables [@federico2023structure; @friedman2008sparse; @friedman2000using; @kalisch2007estimating]. RSNet addresses this limitation by introducing a resampling-based framework that quantifies edge-level uncertainty and integrates information across multiple inferred networks to construct a robust consensus network. The framework supports both Gaussian networks for continuous data and CGBNs for mixed data types and can accommodate correlated or family-based observations [@song2025learning]. This design provides empirical confidence intervals, adjusted p-values, and edge-selection frequencies, offering a fine-grained assessment of network reliability and structure [@friedman2000using; @zhang2018silggm; @jankova2017honest].

In addition to improving reliability, RSNet enhances interpretability through graphlet-based topology analysis, which captures higher-order local connectivity patterns and incorporates edge sign information [@milenkovic2008uncovering; @das2019signed; @doria2020probabilistic]. These functionalities enable detailed examination of node-level structural roles and facilitate comparative analyses between networks inferred under different conditions [@tu2021differential; @gill2010statistical; @higgins2018differential]. Existing R packages for network inference do not support the construction of signed graphlet degree vector matrices (GDVMs), where brute-force enumeration has complexity greater than $O(p^{3})$, with $p$ denoting the network size. RSNet overcomes this barrier by combining state-of-the-art graphlet counting algorithms [@das2020efficient; @das2019algorithm; @hovcevar2016computation] with parallelization, establishing an efficient method to construct signed GDVMs in $O(∣d∣)$, where $∣d∣$ is the average degree, resulting in near-constant time complexity for sparse networks.

# Implementation

## Overview

RSNet provides a unified, resampling-based framework for network inference and analysis that supports both independent and correlated datasets. For independent observations, resampling is performed at the sample level, while for correlated or family-based data, it is conducted at the cluster level to avoid spurious edge discoveries driven by within-cluster dependencies [@song2025learning]. Resampling yields an ensemble of inferred networks, which are subsequently integrated to construct a consensus network that assigns edge-level reliability estimates.

For Gaussian networks, these estimates include empirical confidence intervals and nominal or adjusted p-values; for conditional Gaussian Bayesian networks, they correspond to edge-selection frequencies. The resulting consensus network forms the basis for downstream analyses such as centrality analysis, community detection, graphlet-based topology analysis, and differential connectivity analysis (Figure 1).

![**Overview of RSNet.**  RSNet accepts a sample-by-feature dataset and allows users to specify whether observations are independent or correlated. The resampling framework generates an ensemble of weighted or binary adjacency matrices over m iterations, which are subsequently integrated into a consensus network that serves as the basis for downstream analyses.](figures/Figure1.pdf)

RSNet is organized into modular functions that enable flexible, end-to-end workflows from network inference to higher-order structural analysis. The framework supports parallel computing to accelerate large-scale analyses, making RSNet a reproducible, scalable, and efficient platform for high-dimensional network inference.

## Resampling-based framework

RSNet implements multiple resampling strategies to enhance the stability and reliability of inferred network structures. For both Gaussian networks and conditional Gaussian Bayesian networks (CGBNs), users can choose among four general approaches: (1) unstratified bootstrap, (2) unstratified subsampling, (3) stratified bootstrap, and (4) stratified subsampling, depending on data characteristics and study design. 

For Gaussian networks, RSNet additionally supports cluster-based resampling methods designed for correlated or family-based datasets, including (1) cluster bootstrap, which samples entire clusters with replacement to preserve intra-cluster dependencies, and (2) fractional cluster bootstrap, which samples a subset of clusters with replacement. These procedures are implemented in the function “ensemble_ggm()”,  which leverages inference algorithms from the *SILGGM* package [@zhang2018silggm]. For CGBNs, the function “ensemble_cgbn()”, provides analogous resampling-based network inference using algorithms from the *RHugin* package [@kalisch2007estimating].

The resampling module supports parallel computing, enabling efficient large-scale network inference and ensuring scalability across high-dimensional datasets.

## Consensus network construction

The ensemble of networks generated from resampling iterations is integrated to construct a consensus network, thereby mitigating noise and reducing numerical instability. The Gaussian consensus network is constructed using the function “consensus_net_ggm()”, which computes confidence intervals and both nominal and adjusted p-values for each edge, while supporting thresholding based on effect size and statistical significance. The conditional Gaussian Bayesian consensus network is generated using “consensus_net_cgbn()”, which calculates selection frequencies, the proportion of resampling iterations in which an edge is identified, and outputs a directed acyclic graph (DAG). The method also supports computing Markov-blanket [@friedman2000using; @gao2016efficient] frequencies between incident nodes, which can be used to generate an undirected network representation.

## Graphlet degree vector matrix and graphlet correlation matrix

Graphlets are small connected non-isomorphic induced subgraphs that capture higher-order connectivity patterns within a network [@milenkovic2008uncovering]. For graphlets containing between two and four nodes, there are eight unique structures and fifteen corresponding automorphism orbits (Figure 2A). The graphlet degree vector (GDV) of a node counts how many times the node participates in each orbit, generalizing the concept of node degree to higher-order structures [@milenkovic2008uncovering]. Collectively, these vectors form a graphlet degree vector matrix (GDVM) of dimension $p \times O$, where $p$ denotes the number of nodes and $O$ the number of orbits. 

![**Unsigned and signed graphlets.** (A) The 2-, 3-, and 4-node graphlets G_0,G_1,…, G_8 and their automorphism orbits 0, 1, 2, …, 14. (B) The 2- and 3-node signed graphlets G_0,G_1,…, G_8 and their automorphism orbits 0, 1, 2, …, 14. Positive relationships are represented by black lines, and negative relationships are represented by red lines.](figures/Figure2.pdf)


RSNet extends this framework to signed networks, incorporating edge sign information (positive or negative partial correlations) into graphlet enumeration [@das2020efficient]. Notably, the number of graphlets of size up to four and their orbits in unsigned networks is equivalent to the number of graphlets of size up to three in signed networks (Figure 2B) [@das2019algorithm], allowing a concise but informative representation of signed topologies [@das2020efficient; @das2019algorithm].

The graphlet correlation matrix (GCM) provides a compact, embedding-based representation of network structure [@yaverouglu2014revealing]. It is constructed by computing pairwise Spearman correlations between the columns of the GDVM, yielding an $O \times O$ matrix that captures inter-orbit relationships [@yaverouglu2014revealing]. This embedding enables meaningful comparisons between networks of different sizes or densities on a common structural scale. The functions “gdvm_gcm()” and “signed_gdvm_gcm()” implement the computation of both GDVMs and GCMs for unsigned and signed networks, respectively.

## Differential connectivity analysis

Differential connectivity analysis is an essential application in network-based studies, designed to characterize the rewiring patterns of interactive entities across different conditions [@tu2021differential; @gill2010statistical; @zhang2018diffgraph]. Unlike univariate analyses that evaluate each feature independently, this approach captures integrative, system-level signals arising from the collective behavior of individual entities, thereby providing complementary insights into the underlying mechanisms. 

Accurate estimation of the null distribution is critical for assessing statistical significance [@gill2010statistical]. However, for most network-based statistics, such as average shortest path length, modularity, and centrality measures, the null distribution is challenging to derive analytically [@shojaie2021differential]. Consequently, resampling-based approaches are indispensable for empirically estimating the null distribution [@pollard2004choice; @vavsa2022null]. RSNet supports permutation [@pesarin2010permutation; @van2023comparing]- and bootstrap-based [@karlsson2009bootstrap; @efron1982jackknife] null distribution generation for a broad range of network statistics, including standard metrics (e.g., centralities, modularity, clustering coefficient) and GDV distances to enable differential comparisons at the single-node level. These procedures are implemented in the functions “null_ggm()”, “diff_centrality()”and “diff_gdv()”.


## Other functionalities

In addition to network inference, consensus construction, and graphlet-based analysis, RSNet provides a suite of complementary tools for network analysis. Standard centrality measures, including degree, strength, eigenvector, betweenness, closeness, and PageRank centralities [@barabasi2004network], are implemented in the function “centrality()”, enabling assessment of node importance. RSNet also supports community detection using the Louvain [@guillaume2008fast], Leiden [@traag2019louvain], and Walktrap [pons2005computing] algorithms, implemented in the function “community_detection()”, to identify modular structures within networks.

# Discussion

RSNet provides a versatile and scalable R package for resampling-based network inference, designed to address the challenges of limited sample size in high dimensional datasets. By integrating both Gaussian networks and conditional Gaussian Bayesian networks, RSNet supports structure learning for continuous and mixed data types within a unified framework.
The package enhances interpretability by integrating standard network analysis tools with graphlet-based methods for higher-order topological characterization. To the best of our knowledge, RSNet is the first R package to implement the construction of a signed GDVM in approximately constant time for sparse networks.
To conclude, RSNet facilitates reproducible, statistically robust, and interpretable network analysis. Its modular and parallelized design supports large-scale applications while maintaining transparency and flexibility, making it a user-friendly open-source resource for high-dimensional network inference and comparative structural analysis.


# Availability

The latest version of the RSNet package along with additional information on the installation process can be found on github.com/montilab/RSNet. 

# Funding

This work was supported in part by the National Institutes of Health, NIA cooperative agreements U19 AG023122-16 and UH3 AG064704, and NIDCR R01 R01DE031831.  The findings and conclusions presented in this paper are those of the author(s) and do not necessarily reflect the views of the NIH.

# Reference

