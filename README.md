# LINNAEUS
LINNAEUS is a single-cell organism wide lineage tracing method relying on CRISPR/Cas9-induced genetic marking. This repository contains all scripts used in LINNAEUS as described in "Simultaneous lineage tracing and cell type identification using CRISPR/Cas9 induced genetic scars".

# Tree-building workflow
A description of the full tree-building workflow, starting from sequenced scars and sequenced and (Cell Ranger) mapped mRNA data.

1. __Extract scars.__ Run the extract_scar_P5.sh UNIX shell script in the scar_extraction folder. This script calls scar_CIGAR_sc_10X_v2.pl and scar_filter.py, present in the same folder. The script scar_CIGAR_sc_10X_v2.pl requires [bwa](http://bio-bwa.sourceforge.net) to be installed, and the (Python 2.7) script scar_filter.py assumes the existence of a virtual environment to run in, and packages optparse, pandas and distance to be installed in that environment.
2. __Determine cell types using Seurat.__ Run the Larvae_Seurat_mt.R script in the R_scripts folder. This script combines the Cell Ranger-generated count tables of all sequenced larval 10X libraries, filters, clusters and determines differentially expressed genes for the cells sequenced. Uses [Seurat](http://satijalab.org/seurat/).
3. __Filter scars within the library and determine creation probability of scars.__ Run the Preprocess_droplet_scar_X_allreads.R scripts (one for each organism) and the Scar_comparison_between_fish.R script in the library_scar_filtering folder to combine scar libraries with mRNA libraries and perform further filtering steps.
4. __Run the tree-building algorithm.__ Run Iterative_tree_building.R (for larvae) or Iterative_tree_building_adults.R (for adults) to build and plot tree using cell type and filtered scar data. The data files provided in the directory Test_data can be used to build and display a larval tree.

# Data description
We provide scar and cell type data for one dataset in the Test_data folder, together with colors needed in Iterative_tree_building.R for tree visualization. Other datasets can be downloaded from GEO.

* Larvae_Seurat_batch_r_out_cells_2.csv - Cell type and tSNE coordinates for all larval cells.
* Z2_scars_compared.csv - Filtered scars from one organism.
* color_table_larva.csv - Colors used for subsets of larval cell types to plot subtrees.
* color_table_larvae.csv - Colors used for all larval cell types together.

# Script description
A brief description of all scripts in the repository, ordered by folder.

* __R_scripts__ - Scripts for data analysis and visualization.
	* F1_detection_rates.R - Determination of scar detection rates through F1 experiments.
	* Iterative_tree_building.R - Iterative tree building for larvae.
	* Iterative_tree_building_adults.R - Iterative tree building for adults.
	* Larvae_Seurat_mt.R - Single cell mRNA analysis to determine cell types. Uses [Seurat](http://satijalab.org/seurat/).
	* Link_enrichment.R - Calculations of scar connection enrichment of cell types.
	* Scar_expression_rates.R - Calculations of scar detection efficiency.
	* Scar_simulation_only.R - Simulation of developmental trees with scarring; simulation of a scar-sequencing experiment of cells from these trees.
	* scar_helper_functions.R - Functions written for the analysis scripts.
	* treeviz.r - Visualization of simulated developmental and scar trees, as well as reconstructed simulated LINNAEUS and PHYLIP trees.
* __collapsibleTree__ - A (with author permission) modified version of the [collapsibleTree](https://cran.r-project.org/web/packages/collapsibleTree/index.html) package. Install locally using devtools: ```install_local("./Scripts/LINNAEUS/collapsibleTree/")```.
* __library_scar_filtering__ - Scripts to filter the scars combining all libraries from one organism together, and a script to compare scars found in different organisms.
	* Preprocess_droplet_scar_X_allreads.R - Combine all libraries from one organism together, assess and filter suspicious-looking scars, filter scars with low number of reads, and filter cells with an unusually high number of scars.
	* Scar_comparison_between_fish.R - Count how in how many different organisms a certain scar is observed; this is a proxy for the creation probability of the scars.
* __old__ - Previously used iteration of collapsibleTree
* __scar_extraction__ - Scripts to extract single-cell and bulk scars.
	* scar_CIGAR_bulk.pl - Align bulk scar sequencing data to reference and return a text file with all scars observed. Requires [bwa](http://bio-bwa.sourceforge.net) to be installed.
	* scar_CIGAR_10X_v2.pl - Align 10X single-cell scar sequencing data to reference and return a text file with all scars observed. Requires [bwa](http://bio-bwa.sourceforge.net) to be installed.
	* scar_filter.py - Filter scars for sequencing errors. Assumes the existence of a virtual environment to run in, and packages optparse, pandas and distance to be installed in that environment.
	* extract_scar_P5.sh - UNIX shell script that extracts and filters scars by determining cellular barcodes, aligning sequencing data and filtering scars.
	* extract_scars_bulk1.sh - UNIX shell script that extracts and filters scars by determining cellular barcodes, aligning sequencing data and filtering scars.
	* count_correct_scars.sh - Collapse extracted bulk scars to used barcodes.
