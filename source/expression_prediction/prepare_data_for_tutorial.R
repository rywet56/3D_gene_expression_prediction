library(Seurat)

# preparation for tutorial (simplify by extracting Seurat data to .csv file)
# ---------------
# preparing metadata of cells
path <- "/home/storage/data/analyzed_data/novosparc_paper/day4/so_day4.rds"
so <- readRDS(path)
md <- so@meta.data
umap <- as.data.frame(so[['umap']]@"cell.embeddings")

cell_meta_data <- data.frame(cell_id = rownames(md), 
                                UMAP_1 = umap$UMAP_1, 
                                UMAP_2 = umap$UMAP_2, 
                                cluster_no = md$seurat_clusters)
# create cluster number - annotation mapping
clusters_annotations = c(
  "epidermis", "(pro-)cambium", "cortex", "dividing cells",
  "floral meristem", "S-phase cells", "dividing cells", "mesophyll",
  "xylem parenchyma", "epidermis", "epidermis/dividing", "epidermis/dividing", "phloem"
)
cluster_nos <- paste0("cluster_", 0:(length(clusters_annotations) - 1))
clusters <- paste0(cluster_nos, "___", clusters_annotations)
names(clusters) <-  0:(length(clusters_annotations) - 1)
# assign cells cluster annotations
cell_meta_data$cluster_annotation <- as.character(clusters[as.character(cell_meta_data$cluster_no)])
# save
path <- "/home/storage/scripts/high_level_analysis/novosparc_plus/gene_expression_prediction_3D/data/seurat_cell_meta_data.csv"
write.csv(x = cell_meta_data, file = path)
# ---------------






















