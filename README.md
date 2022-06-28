# 3D_gene_expression_prediction

### Parameters
```diff
+ verbose=True  
--> whether to print messages about the progress of the reconstruction

+ keep_genes=["AT2G17950", "AT2G27250"]  
--> TAIR IDs of genes that should be forced to remain in the scRNA-seq dataset despite being potentially removed during the pre-processing of the data. This will also enforce the prediction of 3D gene expression profiles of those genes.

+ pre_mode = "expr_pre" # or pre_mode = "cell_mapping"  
Setting the prediction mode influences which scRNA-seq cells are used in the cell-to-cell mapping during Optial Transport (OT). Setting this paramter to "cell_mapping" will contain all cells of the scRNA-seq dataset. This makes sense if we are mainly interested in the original location of cells in the scRNA-seq dataset. Setting the parameter to "expr_pre" will allow the selection and/or enrichment of scRNA-seq cells and is reccomended if the main goal is to reconstruct gene expression profiles.



```
