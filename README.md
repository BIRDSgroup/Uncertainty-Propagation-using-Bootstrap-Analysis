# Uncertainty-aware Gene Rankings for Network Measures

This repository contains scripts for implementation of the systematic workflow to generate uncertainty-aware gene rankings for downstream network (degree and PageRank centrality) measures, and their analyses.



## Overview of the Project

This project develops a systematic workflow that estimates and propagates uncertainty about the coexpression network to degree and PageRank centrality and produces robust (uncertainty-aware) gene rankings. 
The workflow uses ‘*bootstrapping*’ to generate multiple measures for the centrality – mean $(\mu)$ and variance $(\sigma^2)$ of which were used to propose different uncertainty-aware gene rankings (BCMs) of the form $\mu-c\sigma$.

The project is implemented using the following programming languages: 

-  *Functions for computations and analyses*: Python versions 3.10 and 3.12. 
   Packages like `numpy, pandas, scipy, random, corals.correlation` and `statsmodel` have been used repeatedly.
-  *Covariate adjustment of the gene expression data and GSEA analysis using WebGestalt*: R versions 4.2 and 4.4
   The WebGestalt package for R `WebGestaltR` has been used primarily.
-  *Workflows combining several scripts*: Bourne Again SHell (BASH) scripting language

Note: Some of the codes have been designed to run parallelly on a server with $50$ cores, however this is only done in the BASH scripts and can be easily modified to run sequentially on a local system.

### Repository Structure

Each folder contains the scripts and results for the same analyses. The codes have been implemented such that the functions will read the files from the corresponding directory structure directly, so that manual duplication of the output files will not be necessary. Hence, we suggest to maintain the same directory structure to run the codes without any hassle. 

The overall directory structure is given below: Modify to put it in a collapsible way.

```python
+---Bias	# Scripts to estimate the Estimator and Bootstrap Bias

+---Centrality Computation # scripts to generate centrality values. 
	'''
	For each tissue T, a folder is generated containing 7 files. 
	Directory structure of the same is given below as an example.
	/---T
	|	error.txt 	# errors during the running of the code.
	|	original_degree.csv # the degree cetrality of the genes (in the same order as in the 'genes.csv' file) by obs.
	|	original_pagerank.csv # the PageRank cetrality of the genes (in the same order as in the 'genes.csv' file) by obs.
	|	output.txt # outputs (if any) pertaining to the computation of centrality for the same tissue.
	|	sample_degree.csv # the degree centrality values of the genes (in the same order as in the 'gnes.csv' file) in each bootstrapped coexpression network.
    |	sample_pagerank.csv # the PageRank centrality values of the genes (in the same order as in the 'genes.csv' file) in each bootstrapped coexpression network.
    |	seed.txt # the seed used to perform bootstrapping on the same tissue (necessary to replicate the results exactly).
	'''

+---Gene Ordering # scripts to generate a random order of the genes
	# For each tissue, a csv file <tissue_name>_order.csv is generated in the same folder

+---Original Dataset
|   +---Covariates 
	# extracted covariate files for each tissue as provided by GTEx.
|   +---Gene Expression Matrices 
	# extracted gene expression files for each tissue as provided by GTEx.
|   +---Genes
	# extracted list of genes for each tissue as provided by GTEx.
|   +---Preprocessed Files # the covariate adjusted files for each tissue.
	'''
	For each tissue T, a folder with 3 files are generated.
	- T.csv: the covariate adjusted coexpression matrix for T
	- T_orig_filtered.csv: the original gene expression matrix for T but only for the protein-coding genes.
	- genes.csv: information related to all protein-coding genes in T. Both degree/PageRank centrality values of all the genes are generated in the same order as in this file. 
	'''

+---Parameter Analysis
	'''
	For every tissue T and for both the centrality measures, it contains scripts to generate:
	1. the ranks for each metric in T.
	2. the Spearman correlation among all genes for T.
	'''
    
+---Replication Analysis on RW Datasets
|   +---Centrality Computation 
	# scripts for centrality comptation of the discovery and real-world datasets
	# similar to the centrality computation on the GTEx tissues.
|   +---Original Dataset
	'''
	Scripts to preprocess and generate the final datasets. Preprocessing involves
	- covariate adjustment of replication dataset
	- extraction of protein-coding genes to both datasets
	'''
|   |   +---Muscle_Skeletal_gtex
|   |   +---Muscle_Skeletal_recount
|   |   +---Preprocessed Files
|   +---Results 
	# scripts to compute the gene ranking by each metric and generate the plots.

+---Simulations and Validation Analysis
	'''
	This contains scripts for two analyses:
	1. Generation of simlated data and all analyses on the same.
		- all of these scripts are named as 'SIM_*.py'.
	2. Validation analysis of centrality measures on GTEx tissues.
		- these scripts compute the gene ranking by each metric and generate the plots.
		- all of these scripts are named as 'RW_*.py'
	'''

+---Tissue Specificity Analysis - 0.05 FDR
	# Specificity Analysis with 5% FDR
|   +---Elevated Genes # list of specific genes as downloaded from Human Protein Atlas
|   +---Parameter Analysis 
	# scripts to generate correlations of specific genes and other genes. 
|   +---Pathway Analysis
	'''
	Stores results of the GSEA analysis by WebGestalt. 
	E.g., the output directory structure for centrality C and tissue T is given below.
	+---Elevated_C
	|	+---T
	|	|	+---Project_deg	# results from gene ranking of obs
	|	|	+---Project_mu
	|	|	+---Project_mu_2sigma
	|	|	+---Project_mu_sigma
	+---Pathway Results # scripts to generate the heatmaps and boxplots.
	'''
    
+---Tissue Specificity Analysis - 0.05 FDR
	# Specificity Analysis with 10% FDR
|   +---Elevated Genes # list of specific genes as downloaded from Human Protein Atlas
|   +---Parameter Analysis 
	# scripts to generate correlations of specific genes and other genes. 
|   +---Pathway Analysis
	'''
	Stores results of the GSEA analysis by WebGestalt. 
	E.g., the output directory structure for centrality C and tissue T is given below.
	+---Elevated_C
	|	+---T
	|	|	+---Project_deg	# results from gene ranking of obs
	|	|	+---Project_mu
	|	|	+---Project_mu_2sigma
	|	|	+---Project_mu_sigma
	+---Pathway Results # scripts to generate the heatmaps and boxplots.
	'''
```

## Input Files and Executing the Pipelines

### Generating the centrality measures

Before generating the centrality files, please ensure that the gene expression matrix ($\textrm{genes} \times \textrm{samples}$) is stored in a  comma-separated `.csv` file with the gene id’s as the first column and the sample ids as the first row.

To generate the centrality measures for both degree and PageRank centrality and to replicate our the results, run `Centrality Computation/run.sh`

1. Ensure the working directory is `Centrality Computation`.
2. Set the value of `B` with the number of bootstrap coexpression networks you want to generate. 
   Ensure `B` is set to a multiple of the number of cores you are using. For you are running the codes sequentially, then this constraint is irrelevant.  
3. Set `tissues` with the list of tissues for which you want to generate – each tissue name should be enclosed within quotes and separated by space.
4. Run the file from the terminal or GIT Bash (in a windows system) as `./run.sh`.

To run the program on a local system, edit `Centrality Computation/main.sh` as follows:

1. Comment line number $32$ `((counter++)); ((counter % cores == 0)) && wait`.
2. Remove `&` from line number $31$ `time python compute.py $tissue $((b-1))`.

Edit `Centrality Computation/compute.py` (and `Centrality Computation/compute.py` for values of the $\mathrm{obs}$ metric) under the following scenarios:

-  *You want to generate the centrality values for a specific centrality measure* – comment / remove the block corresponding to the centrality measure that you don’t want to generate.
-  *Your input/preprocessed file is in a different location* – update the `in_file` variable with the location of your file.
-  *Your input/preprocessed file has a different structure as mentioned above* – please update the code to read your gene expression matrix and enter the code to update the structure of the data frame before invoking the `bootstrap.BootstrapSample()` function.

### Validation Analysis on GTEx Datasets

To generate the plots for the validation analysis, run `Simulations and Validation Analysis/RW_validation.py` with the necessary inputs. 

For e.g., to generate the plots for centrality `C` ( `degree` or `pagerank` ) for `Lung`  tissue on the observed dataset with $237$ samples, do the following

1. Ensure the working directory is `Simulations and Validation Analysis`.
2. Ensure the folder `Semi-Simulated Data/Lung/C` exists in the same directory.
3. Set `klim` – the x-axis limit of the CAT and recall plots.
4. Run the code from the terminal as: `time python RW_validation.py Lung C 237`

To edit the code for the following scenarios:

-  *If you want to generate only the CAT or recall plot* – comment / remove the corresponding blocks for the same.
-  *If you want the output results in the CAT and/or recall plots in a separate location* – update the `filepath` variable for the corresponding plot with your destination location.
-  *If your input files are in a separate location* – then the codes to read the `.csv` files using `read_csv()` function from `pandas` in lines $27, 30, 33$ and $35$ may be updated. 
-  Please *update the parameters* of the `Plot()` function in the `plots` file as per your requirements.
-  To generate the POG values, instead of the plots, please run `RW_pog_vals.py` from the same directory with the necessary modifications as already mentioned.

### Replication Analysis on Real-World Datasets

To perform the procedure is similar as on the GTEx datasets, however, here, the main directory should be `Replication Analysis on RW Datasets`.

1. *For covariate adjustment of the replication dataset*:
   1. Ensure the current working directory is `Replication Analysis on RW Datasets/Original Dataset\Muscle_Skeletal_recount`.
   2. To extract the covariates: run `time python cov extract.py` from the terminal.
   3. To perform the adjustment run the `covariates.R` script.
2. *For extraction of genes common to the datasets* – run `genes.py` from `Replication Analysis on RW Datasets`.
3. Following this, the *centralities may be similarly computed* as mentioned above from `Replication Analysis on RW Datasets/Centrality Computation`.
4. To *plot the results* (or to generate the POG values), please run `Replication Analysis on RW Datasets/Results/results.py` (or `pog_vals.py`) using the same procedure as in validation analysis.

##  Simulated Dataset

For simulated datasets, the codes for each analysis is independent of each other. Each analysis generates the population dataset, and performs the analysis as required. Nevertheless, if the seed is the same, the same population and respective observed datasets will be generated irrespective of the script/analysis. Please ensure the current working directory is `Simulations and Validation Analysis`.

-  To generate the heatmaps of the adjacency matrix, the degree distribution and the $\mu$ vs. $\sigma$ plots, run `SIM_graphs.py`.
-  To perform validation analysis and plot the CAT and recall plots, run `SIM_validation.py`.
-  To perform replication analysis and plot the CAT and recall plots, run `SIM_replication.py`.
-  To generate the POG values in the validation and replication analysis, run `SIM_pog_vals_1.py` and `SIM_pog_vals_2.py` respectively.

The input number of genes, the number of bootstrap samples, sample sizes of the observed dataset and the x-axis limit of the CAT and recall plots can be modified with `n, B, SampleSizes` and `klim` variables respectively. In the scripts for replication analysis, set the value of `s` to modify the size of the discovery dataset.

All of modifications, if required, can be performed as in the scripts for the real-world datasets.

## License

This project is licensed under the Creative Commons License — see [Link to License file] for details.

## Citation, Contact

To be updated following paper acceptance.

## Acknowledgements

We thank Dr. Sanga Mitra for her invaluable suggestions and contribution towards the completion of this work.
