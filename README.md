# 3D_gene_expression_prediction

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



num_neighbors_source
num_neighbors_target
alpha
epsilon

max_iter
tol
