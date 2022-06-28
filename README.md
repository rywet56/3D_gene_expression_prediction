# 3D_gene_expression_prediction

### Enviornment Setup

### Predicting 3D gene expression profiles

### Predicting cell-to-cell mappings

### Parameters
``` verbose = True # or False ```  
Whether to print messages about the progress of the reconstruction.

``` keep_genes = ["AT2G17950", "AT2G27250"] ```  
TAIR IDs of genes that should be forced to remain in the scRNA-seq dataset despite being potentially removed during the pre-processing of the data. This will also enforce the prediction of 3D gene expression profiles of those genes.

``` pre_mode = "expr_pre" # or pre_mode = "cell_mapping" ```  
Setting the prediction mode influences which scRNA-seq cells are used in the cell-to-cell mapping during Optial Transport (OT). Setting this paramter to "cell_mapping" will contain all cells of the scRNA-seq dataset. This makes sense if we are mainly interested in the original location of cells in the scRNA-seq dataset. Setting the parameter to "expr_pre" will allow the selection and/or enrichment of scRNA-seq cells and is reccomended if the main goal is to reconstruct gene expression profiles.

``` method_dist = "euclidean" ```  
The distance measure that should be used in the context of OT to calculate distances between cells.  

``` top_sccells = 50 ```  
The number of cells in which a gene has be expressed in the scRNA-seq DGE in order to be kept during the pre-processing of the scRNA-seq dataset.  

``` enrich = False # or True ```  
If pre_mode = "expr_pre", this option can be used to not only select cells in the scRNA-seq dataset, but also to enrich some of those selected cells in such a way that the scRNA-seq dataset best represents the composition of the cells in the spatial expression dataset.  

``` genes_remove_list = [] ```  
Setting this parameter to a list of TAIR IDs allows the user to remove genes from the spatial gene expression dataset. The removed reference genes will not be used in the cell-to-cell mapping and prediction of 3D gene expression profiles. This may be useful if the user wants to screen the effect of gene removal from the spatial reference dataset on gene expression prediction performance.  

``` num_neighbors_source ```  
The number of nearest neighbors (cells) that should be considered for constructing a nearest neighbor graph of the cells in the scRNA-seq based on gene expression. If this number is small, the nearest neighbor graph describes more the local structure of cells instead of the global structure of cells.  

``` num_neighbors_target ```  
The number of nearest neighbors (cells) that should be considered for constructing a nearest neighbor graph of the cells in the spatial expression dataset based on their positional information. If this number is small, the nearest neighbor graph describes more the local structure of cells instead of the global structure of cells.

``` alpha = 0.1 ```  
This parameter expresses how much the cell-to-cell mapping (and therefore the 3D gene expression prediction) should be driven by either the information about the spatial position or the reference gene expression of cells in the spatial reference dataset. A value of 0 uses only information from the spatial position of cells to perform the cell-to-cell mapping between cells in the scRNA-seq dataset to cells in the spatial expression dataset. Setting alpha to 1, uses only the information about the gene expression of marker genes to obtain cell-to-cell mappings. A value for alpha of 0.5 balances those two sources of information.  

``` epsilon = 0.05 ```  
A regularization constant. The larger this value, the lower the reconstruction performance.

``` max_iter = 500```  
The maximum number of iterations of the underlying OT-based reconstruction algorithm. The larger this value, the longer the alogrithm takes to converge.

``` tol = 1e-9```  
The tolerance for termination of the underlying OT-based reconstruction algorithm. The smaller this value, the longer the alogrithm takes to converge.  


