# Phylogenetics-Project

## Introduction
This is a GitHub repository for my 3rd year phylogenetics research project. In this project I am comapring sequence data to analyse whether the seqeunces can be expressed as a phylogenetic tree.

## Analysis
In this project I am using a NEXUS file format to either store sequence data or distance matricies. From this file, I am computing the overall delta value, for every possible quartet of organisms inputted (Holland et al., 2002). I can also alter this method slighltly to be able to calculate the mean delta value for each taxa, and therefore find the taxa that are causing the non-tree like data. This is due to processes of recombination, horizontal gene transfer and parallel evolution, which can therefore effect the typical 'tree-like' evolution (Holland et al., 2002) .

## Methods
The plotting of delta plots will be carried out on python and the plotting of phylogenetic networks will be
carried out using the SplitsTree package (Huson and Bryant, 2005). I will need to locate past sequence comparison data when
forming the distance matrix or find and compare sequences from different species manually online. Where I
will also need to produce a method for random sampling of the sequence data into quartets. After carrying
out the analysis, I will also try and find a new scaling method for the delta values therefore producing more
accurate delta plots.

## Refrences
Huson, D.H. and Bryant, D. (2005). Application of Phylogenetic Networks in Evolutionary Studies. Molecular Biology and Evolution, 23(2), pp.254–267. doi:https://doi.org/10.1093/molbev/msj030.

Holland, B.R., Huber, K.T., Dress, A. and Moulton, V. (2002). δ Plots: A Tool for Analyzing Phylogenetic Distance Data. Molecular Biology and Evolution, 19(12), pp.2051–2059. doi:https://doi.org/10.1093/oxfordjournals.molbev.a004030.
