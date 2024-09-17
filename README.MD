# Normalization effects on dimensionality reduction and clustering efficiency

This repository investigates how various normalization techniques impact the efficiency of dimensionality reduction and clustering in single-cell RNA sequencing (scRNA-seq) data analysis. We explore and compare multiple normalization methods, apply dimensionality reduction techniques (t-SNE and UMAP), and leverage clustering algorithms to evaluate their effects.

## Initialization 

All the package/library dependencies are stored in `requirements.R` file. 

Enviromental variables are stored within config.yml file

## Repository Structure

The repository is organized as follows:

* `src/data/`: Contains scRNA-seq datasets, including liver, breast cancer, and peripheral blood mononuclear cell (PBMC) datasets.
    * `src/data/normalized`: Consists of normalized data.
    * `src/data/reduced`: Consists of normalized data after dimensionality reduction.
    * `src/data/clustered`: Consists of processed data, after clustering.
* `src/normalization/`: Scripts for applying various normalization methods. Each dataset has its corresponding normalization script.
* `src/dim_reduction/`: Scripts for performing dimensionality reduction using both t-SNE and UMAP techniques.
* `src/clustering/`: Contains scripts for clustering the data using different clustering algorithms.
* `src/functions/`: Contains utility functions used across normalization, dimensionality reduction, and clustering scripts.

## Details

### Normalization Methods Compared

| Normalization | R package | url |
| ------------- | --------- | --- |
| log2 | base      | -   |
| Freeman-Tukey | base | - |
| Seurat | Seurat: Tools for Single Cell Genomics | [CRAN database](https://cran.r-project.org/web/packages/Seurat/index.html) |
| scran | Methods for Single-Cell RNA-Seq Data Analysis | [Bioconductor database](https://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html) |
| Dino | Normalization of Single-Cell mRNA Sequencing Data | [Bioconductor database](https://www.bioconductor.org/packages/release/bioc/html/Dino.html) |
| SCnorm | Robust normalization of single-cell RNA-seq data | [GitHub repository](https://github.com/rhondabacher/SCnorm) |

### Dimensionality reduction techniques used

| Dimensionality reduction | 
| ------------------------ |
| t-Distributed Stochastic Neighbor Embedding |
| Uniform Manifold Approximation and Projection for Dimension Reduction | 

### Leveraged clustering algorithms 

| Clustering Algorithm | R Package |
| -------------------- | --------- |
| Hierarchical Clustering | stats | [R doc](https://www.rdocumentation.org/packages/stats/versions/3.6.2) |
| K-means | stats |  [R doc](https://www.rdocumentation.org/packages/stats/versions/3.6.2) |
| Louvain Community Detection | Seurat | [CRAN database](https://cran.r-project.org/web/packages/Seurat/index.html) |
| scCAN | scCAN | [CRAN](https://cran.r-project.org/web/packages/scCAN/index.html) |


### How to Use

1. Preprocessing the Data:

* Place scRNA-seq datasets into the src/data/ directory.
* Perform normalization within src/normalization/ directory.

Normalization:

* Run the normalization scripts to apply the selected normalization methods to your datasets.

3. Dimensionality Reduction:

* Use the scripts in src/dim_reduction/ to perform t-SNE or UMAP on the normalized data.

4. Clustering:

* Apply the clustering algorithms by running the scripts in src/clustering/.

5. Analysis:

* TODO: Evaluate the impact of each normalization method on clustering and dimensionality reduction performance using ARI, SI, Jaccard Index...

**_NOTE:_** There is [scCAN](https://www.nature.com/articles/s41598-022-14218-6) clustering method missing, that will soon updated... 
