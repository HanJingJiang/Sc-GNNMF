# Sc-GNNMF
Single-cell RNA sequencing (scRNA-seq) technology's advent has brought forth fresh perspectives on intricate biological processes, revealing the nuances and divergences present among distinct cells.
Accurate single-cell analysis is a crucial prerequisite for in-depth investigation into the underlying mechanisms of heterogeneity. 
Due to various technical noises, such as the impact of dropout values, the analysis of scRNA-seq data remains challenging.
In this work, we propose an unsupervised learning framework for scRNA-seq data analysis, named Sc-GNNMF. 
Based on the non-negativity and sparsity of scRNA-seq data, we propose employing graph-regularized non-negative matrix factorization (GNMF) algorithm for dimensionality reduction of scRNA-seq data, which involves estimating cell-to-cell similarity and gene-to-gene similarity through Laplacian kernels and p-nearest neighbor graphs (p-NNG).
By assuming intrinsic geometric local invariance, we use a weighted p-nearest known neighbors (p-NKN) of cell-cell interactions to guide the matrix decomposition process, promoting the closeness of cells with similar types in the cell-gene data space and determining a more suitable embedding space for clustering.
Experiments on 9 real scRNA-seq datasets demonstrate that Sc-GNNMF outperforms the other seven methods and maintains satisfactory compatibility and robustness. 
Furthermore, Sc-GNNMF yields excellent results in  clustering tasks, extracting useful genetic markers, and pseudo-temporal analysis.
