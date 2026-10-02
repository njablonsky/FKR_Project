# FKR_Project
Associated code for "Genetic exchange drives relatedness of large linear and circular extrachromosomal elements of archaea" by Jablonsky et al. (2026). ChatGPT and Claude were both used to assist in the writing of this code. 

For development of the structural clustering matrix (Fig 5):
- tmalign_matrix.py = Generates matrix of TM-scores from all-against-all pairwise protein structure comparisons. Written by NJJ
- plot_matrix.R = Visualises the resulting matrix in R. Written by NJJ

For Kendall-Tau analysis:
- Kendall_Tau_conservation.py = Generates Kendall-Tau scores for genome pairs. Written by L-DS.
- Order_conservation_plots.R = Generates matrix of Kendall-Tau scores and linear map of SCG position in R. Written by L-DS.

For subfamily enrichment analysis / Fisher test:
- Fisher_test.R = Identifies subfamilies that are enriched in certain ECE groups over others and generates matrix of significant subfamilies. Written by NJJ
