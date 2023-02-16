A code pool for RNA expression analysis in Woodman Lab.

## Installation

Without SSH key: 

    devtools::install_git(
        url='http://gitlab.mdanderson.edu/XLi23/expr.git',
    )

With SSH key set up:

    devtools::install_git(
        url='git@gitlab.mdanderson.edu:XLi23/expr.git',
        quiet=FALSE
    )
 
  - [How to set up SSH key?](https://docs.gitlab.com/ee/user/ssh.html)

## Usage
View package vignette with `browseVignettes("expr")` in Rstudio console, or [here](http://127.0.0.1:10041/library/expr/doc/expr_vignette.html)

## Development
- This package is documented with [roxygen2](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html). Before pushing a commit, 'Cmd + Shift + B' to generate document, then 'Cmd + Shift + B' to build the package. Keep new functions going!
- Package tests are as good as how you write it. Please kindly report issues under "issue" tab every time an error is encountered.

## Notes 11/13/2023
**functions**
  1) prepare_clean_RNA_sample_info_and_protein_expressions
  2) prepare_unsupervised_data
  3) unsupervised_analysis

- Added line in `prepare_clean_RNA_sample_info_and_protein_expressions` function to ensure expression column names are identical to sample info sample IDs so that otherwise function error out.
- Param name in `prepare_clean_RNA_sample_info_and_protein_expressions` changed `tolerent_library_size_factor` to `tolerant_library_size_factor`.
- `Match.arg` renders error when arg has a length of 1 (when running function line by line). Changed to select first element in a vector.
- gg3D (needed for `unsupervised_analysis`) seems to be difficult to install (XQuarts needed for MAC), used plotly instead.
- The `guide` argument in `scale_*()` cannot be `FALSE`. This was deprecated in ggplot2 3.3.4. Fixed by adding `legend.position="none"`
- Fixed the error that setting any of analysis to `FALSE` fails the `unsupervised_analysis` function.
  1) grid.arrange line for organizing figures was removed. See function manual for arranging plots.
  2) Results were sperated into plots and analysis responses in the function output.

**vignette**
- The package vignette showcase exprClean in making RNA analysis a pipeline.  Added a k-mean clustering measure to auto-detect batch effect as an alternative for manual visualization. Manual confirmation still recommended.
- The package vignette can be easily customized to generate one-click html or pdf report for user datasets.

## Notes 2/15/2023
**Added exported functions**
- NbClust
- consensus_immunedeconvolute
- estimate_bestNumberofClusters
- map_clusters
- consensusCluster
- gsea
- MRNsurr
- getHeatMapAnnotation
- reassignNA
- table_org
- zscoreData
- changeColNames
- getFill
- t_test2

**Edits**
1. getHeatmapAnnotations: added param track - for selecting meta data of interest as track input.
2. Moved “protein coding ensemble to symbol” file and “hg19 ncbi protein coding gene info” to `data/ ` as RData to reduce package size.

**Fixed**
1.  In function `estimate_bestNumberofClusters()`, index_for_NbClust<- c(…,“rubbin”,..)should be “rubin”.
2. In function `map_clusters()`: used `apply(cluster_table,1,max)` instead of `rowMaxs(cluster_table), ` to reduce package dependencies.
3. In function `multiCluster()`: fixed below error by individual indices (adding param “allow1=T” in sourced `mcl_clusters()` ).
    1. `Error in if ((res[ncP - min_nc + 1, 15] <= resCritical[ncP - min_nc +  : missing value where TRUE/FALSE needed” `
