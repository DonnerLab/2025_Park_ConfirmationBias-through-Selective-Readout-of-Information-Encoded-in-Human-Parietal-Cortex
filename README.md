# 2025_Park_Confirmation_Bias_through_Selective_Use_of_Evidence_in_Human_Cortex

MATLAB codes used to compute information measures and behavioral metrics in article: <link to paper, DOI>
MATLAB version: MATLAB2020b, MATLAB2024

The Figure indices corresponding to figures in the manuscript are included in the script name. 
Script Figure4_*.m plots for main Figure4, and also Supplementary Figures S8, S9 and S10. 
The source data for plotting the figures can be found here at the University of Hamburg 
Center for Sustainable Data Management: https://doi.org/10.25592/uhhfdm.16918

plot_maps folder contains the necessary data and code for plotting brain maps with 180 cortical parcels, with the Glasser atlas (Glasser et al. A multi-modal parcellation of human cerebral cortex, Nature. 2016 Aug 11;536(7615):171–178. doi: 10.1038/nature18933)


Computation of the information measures which are central to this project was done with the MINT:
https://github.com/panzerilab/MINT
BIORXIV article: 
MINT: a toolbox for the analysis of multivariate neural information coding and transmission
Gabriel Matías Lorenz, Nicola M. Engel, Marco Celotto, Loren Kocillari, Sebastiano Curreli, Tommaso Fellin, Stefano Panzeri
https://doi.org/10.1101/2024.11.04.621910, PLOS Computational Biology, in press. 

Scripts which compute the information measures using MINT with source data (principal components, can be found here: https://www.fdr.uni-hamburg.de/deposit/xxxx):
  Intersection information: confbias_MINTII_int12_cic.m 
  Mutual information: confbias_MINTMI_int12_SE.m

Script which computes the intersection information linear proxy, corr(s^hat, E) in Figure S11:
  confbias_regress_sample_int12_cic.m
 
