#run the workflow using this command
Rscript run_WGCNA.R


#Input File
#expression_matrix_cortex.csv	#RNA-seq gene expression matrix downloaded from Brainspan

#Out Files
#WC2_4.txt	#Filtered gene expression matrix for period 2 to 4
#WC2_4_power.txt	#Power selectionto reach Scale Free Topology
#WC2_4_module.txt	#Module assignment for genes
#WC2_4_all_con.txt	#Connectivity scores of genes
#WC2_4_cmd.txt	#Multidimensional Scaling of genes
#WC2_4_module3_cor0.95_network.txt	#Network of genes by correlation cutoff of 0.95 for module 3
#WC2_4_eigengene.txt	#Module eigengene expression
#WC2_4_cmd_formodules.txt	#Multidimensional Scaling of modules based on eigengene expression
#WC2_4_cor0.3_modulecor.txt	#Network of modules by correlation cutoff of 0.5

