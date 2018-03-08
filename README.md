# LINNAEUS
LINNAEUS is a single-cell organism wide lineage tracing method relying on CRISPR/Cas9-induced genetic marking. This repository contains all scripts used in LINNAEUS as described in "Simultaneous lineage tracing and cell type identification using CRISPR/Cas9 induced genetic scars".

# Tree-building workflow
A description of the full tree-building workflow, starting from sequenced scars and sequenced and (Cell Ranger) mapped mRNA data.

1. Extract scars using the extract_scar_P5.sh UNIX shell script in the scar_extraction folder. This script calls scar_CIGAR_sc_10X_v2.pl and scar_filter.py, present in the same folder. The script scar_CIGAR_sc_10X_v2.pl requires [bwa](http://bio-bwa.sourceforge.net) to be installed, and the (Python 2.7) script scar_filter.py assumes the existence of a virtual environment to run in, and packages optparse, pandas and distance to be installed in that environment.
2. Determine cell types using Seurat.
3. Filter scars within the library and determine creation probability of scars.
4. Run the tree-building algorithm.

# Script description
A brief description of all scripts in the repository, ordered by folder.

* R_scripts
* collapsibleTree - a (with author permission) modified version of the [collapsibleTree](https://cran.r-project.org/web/packages/collapsibleTree/index.html) package. Install locally using devtools: ```install_local("./Scripts/LINNAEUS/collapsibleTree/")```.
* library_scar_filtering - scripts to filter the scars combining all libraries from one organism together, and a script to compare scars found in different organisms.
	* Preprocess_droplet_scar_X_allreads.R - Combine all libraries from one organism together, assess and filter suspicious-looking scars, filter scars with low number of reads, and filter cells with an unusually high number of scars.
	* Scar_comparison_between_fish.R - Count how in how many different organisms a certain scar is observed; this is a proxy for the creation probability of the scars.
* old - previously used iteration of collapsibleTree
* scar_extraction - scripts to extract single-cell and bulk scars.
