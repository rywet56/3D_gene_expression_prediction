import numpy as np
import novosparc as nc
from scipy.spatial.distance import cdist
import pandas as pd


def load_scDGE(dge_path):
    # load DGE
    dge = pd.read_csv(filepath_or_buffer=dge_path, sep=',', index_col=0)
    # get gene names (rows)
    gene_names = list(dge.index)
    gene_names = np.array(gene_names)
    # get cell ids (columns)
    cell_ids = list(dge.columns)
    cell_ids = np.array(cell_ids)
    # convert pandas.DataFrame to numpy.array
    dge = dge.to_numpy(dtype=np.float32)

    return dge, gene_names, cell_ids

def load_spatial_expression_matrix(coord_path):
    coor = pd.read_csv(coord_path, sep=",", index_col=0, decimal=".")
    return coor

def pre_process_scDGE(dge, gene_names, min_sccells_gene_expressed, keep_genes):
    keep = np.sum(dge > 0, axis=1) # get number of cells in which each gene is expressed
    keep = np.argwhere(keep > min_sccells_gene_expressed) # indices of genes that are at least expressed in "min_sccells_gene_expressed" cells
    keep = np.array([k[0] for k in keep]) # convert to 1D np.array
    genes_name_keep = np.unique(np.append(gene_names[keep], keep_genes)) # get gene names to be kept also considering genes in "keep_genes" that user forces to be kept
    genes_index_keep = []
    for gene in genes_name_keep:
        if gene in gene_names:
            genes_index_keep.append(np.where(gene == gene_names)[0][0])
    dge = np.transpose(dge) # transpose "dge" to now have shape (cells x genes)
    dge = dge[:, genes_index_keep] # only keep genes in "dge" that where selected
    
    return dge, genes_name_keep

def pre_process_spatial_expression_matrix(sem, genes_name_keep, genes_remove_list):
    sel_genes_1 = np.intersect1d(sem.columns.values, genes_name_keep) # get genes in ref. dataset that were selected to be kept in scRNA-seq DGE
    genes_keep_index = np.logical_not(np.isin(sel_genes_1, genes_remove_list)) # get boolean list indicating which genes in ref. dataset should be kept (considering a user specified list of genes that should be removed) during enrichment and prediction
    sel_genes_2 = sel_genes_1[genes_keep_index]  # use boolean list (genes_keep_index) to get final set of genes in ref. dataset that should be retained
    
    sel_genes = sel_genes_2
    
    insitu_matrix = sem[sel_genes] # generate new ref. dataset only containing selected ref. genes
    insitu_matrix = insitu_matrix.to_numpy(dtype=float) # convert "insitu_matrix" to numpy array
    coord = sem.to_numpy() # create another copy of the ref. dataset as a numpy array that will (after the next step) serve as the pure coordinate matrix
    coord = coord[:, [0, 1, 2]] # only keep the first three columns of the ref. dataset to obtain a ref. dataset with only the coordinates of cells
    
    return insitu_matrix, coord, sel_genes

def get_idx_ref_genes_in_scDGE(dge, genes_name_keep, sel_genes):
    index_genes = []
    for x in sel_genes:
        index_genes = np.append(index_genes, np.nonzero(np.in1d(genes_name_keep, x) + 0)[0])
    index_genes = index_genes.astype(int) # convert float to int
    return(index_genes)

def get_best_cells(dge, index_genes, insitu_matrix, method_dist, top_sccells, enrich):
    '''
    dge: pandas scDGE (cells x genes)
    index_genes: indices of sel ref. genes in pandas scDGE
    insitu_matrix: 2D numpy array containing expression of reference genes (cells x genes)
    return: list of indices with respect to orig. scDGE that should be kept
    '''
    # indices of original scDGE
    l = np.arange(0, dge.shape[0], 1)

    # Get indices of cells in scRNA-seq DGE which express at least one ref. gene
    # while only considering genes that are in the ref. dataset (after selecting genes in the ref. dataset)
    dge_selgenes = (dge[:, index_genes] > 0) + 0 # get scRNA-seq DGE that contains sel. ref. genes and then mark cells that have an 
                                                 # expression value for a particular gene of at least one
    keep_pattern = np.nonzero((np.sum(dge_selgenes, axis=1) > 0) + 0)[0] # get indices of cells that express at least one of the
                                                                         # ref. genes
    l_1 = l[keep_pattern]  # original indices of cells to be retained
    
    dge_bin_norm = transform_2D_array(array2d=dge, transform="binary", normalize=True) # gives binarized scRNA-seq DGE, normalize=True has no effect here
    insitu_matrix_bin = transform_2D_array(array2d=insitu_matrix, transform="binary", normalize=False) # gives binarized scRNA-seq DGE, only changes if insitu_matrix is not already binary

    # Get indices of cells in scRNA-seq DGE wich are most similar to cells in insitu matrix
    topcells = []
    dist_cells = cdist(insitu_matrix_bin, dge_bin_norm[np.ix_(keep_pattern, index_genes)], metric=method_dist) # calculate distances between all cells
    if method_dist in ["hamming", "euclidean"]:
        # if we use hamming or euclidean distance, negate the values since we will order them from 
        dist_cells = -dist_cells
    for x in dist_cells:
        # go over each cell in the ref. dataset and it's distances to cells in the scRNA-seq DGE
        # only keep top_sccells cells with lowest distances or highest similarity, respectively
        if(dist_cells.shape[1] < top_sccells):
            # it is possible that there are less scRNA-seq cells remaining than the user-specified top_sccells to be kept
            # in this case just use as many scRNA-seq cells as there are left.
            cells = np.argsort(x, kind="stable")[range(dist_cells.shape[1])]
        else:
            cells = np.argsort(x, kind="stable")[range(top_sccells)] # order from small to big and only keep "top_sccells" cells
        topcells = np.append(topcells, cells)
    if not enrich:
        print("no enrichment")
        topcells = np.unique(topcells.astype(int))
        print(len(topcells))
    else:
        print("cell enrichment performed")
        topcells = topcells.astype(int)
        print(len(topcells))
    l_2 = l_1[topcells]  # indices of cells with respect to the original matrix that should be kept

    return l_2

def transform_2D_array(array2d, transform, normalize):
    if transform == "binary":
        array2d = (array2d > 0) + 0
    if normalize:
        array2d = array2d / np.amax(array2d) # normalize by highest value in array2d (so the whole matrix)
    return array2d

def subset_dge_to_informative_features(dge, topcells, genes_name_keep, sel_genes, top_cor_genes, pre_mode):
    dge_enr = dge[topcells, :] # only use topcells for calculating finiding relevant features
    dd = np.corrcoef(dge_enr, rowvar=False) # calculate correlation of genes
    dd = dd[np.in1d(genes_name_keep, sel_genes), :]
    hvg = []
    for x in dd:
            # go through each row in correlation matrix "dd"
        temp = np.argsort(-x, kind="stable")[range(top_cor_genes)] # sort from small to big to get top "top_cor_genes" correlated genes
        hvg = np.append(hvg, temp)
    hvg = np.unique(hvg.astype(int)) # make list unique
    # depending on the pre_mode, produce matrix with all or only topcells
    if pre_mode == "cell_mapping":
        dge_hvg = dge[:, hvg]
    elif pre_mode == "expr_pre":
        dge_hvg = dge_enr[:, hvg]
    return dge_hvg


def get_expression_location_cost(dge_hvg, coord, ns, nt):
    cost_expression, cost_locations = nc.rc.setup_for_OT_reconstruction(dge_hvg,
                                                                        coord,
                                                                        num_neighbors_source=ns,
                                                                        num_neighbors_target=nt)
    return cost_expression, cost_locations

def get_distributions_over_expression_location(dge_hvg, coord):
    num_locations = coord.shape[0]
    num_cells = dge_hvg.shape[0]
    p_locations, p_expression = nc.rc.create_space_distributions(num_locations, num_cells)
    return p_locations, p_expression

def get_marker_gene_cost(dge, insitu_matrix, pre_mode, tr, index_genes, method_dist, topcells=None):
    dge_bin_norm = transform_2D_array(array2d=dge, transform=tr, normalize=True) # binarize and normalize DGE if specified
    if pre_mode == "cell_mapping":
        dge_enr_selgenes = dge_bin_norm[:,index_genes]
    elif pre_mode == "expr_pre":
        dge_enr_selgenes = dge_bin_norm[np.ix_(topcells, index_genes)]
    insitu_matrix_bin_norm = transform_2D_array(array2d=insitu_matrix, transform=tr, normalize=True)
    cost_marker_genes = cdist(dge_enr_selgenes, insitu_matrix_bin_norm, metric=method_dist)
    return cost_marker_genes

def predict_cell_to_cell_mappings(cost_marker_genes, cost_expression, cost_locations, p_expression, p_locations, 
                    alpha, epsilon, max_iter, tol, verbose=False):
    gw = nc.rc._GWadjusted.gromov_wasserstein_adjusted_norm(cost_mat=cost_marker_genes,
                                                            C1=cost_expression, C2=cost_locations,
                                                            alpha_linear=alpha,
                                                            p=p_expression, q=p_locations,
                                                            loss_fun='square_loss', epsilon=epsilon,
                                                            max_iter=max_iter, tol=tol,
                                                            verbose=verbose, log=False)
    return gw * 10e5

def predict_spatial_gene_expression(gw, dge, topcells, pre_mode):
    if pre_mode == "cell_mapping":
        dge_use = dge
    elif pre_mode == "expr_pre":
        dge_use = dge[topcells, :]
    sdge = np.dot(dge_use.transpose(), gw) # actual reconstruction of spatial expression profiles
    return sdge

def prepare_reconstruction(pre_mode,
                              min_sccells_gene_expressed, keep_genes, genes_remove_list,
                              method_dist, top_sccells, enrich, top_cor_genes,
                              num_neighbors_target, num_neighbors_source, transform,
                              dge, gene_names, sem):
    
    '''
    Prepares all data needed to perform OT-based mapping of cells and gene expression reconstruction.
    '''
    
    dge, genes_name_keep = pre_process_scDGE(dge=dge, gene_names=gene_names, 
                                             min_sccells_gene_expressed=min_sccells_gene_expressed, 
                                             keep_genes=keep_genes)
    
    insitu_matrix, coord, sel_genes = pre_process_spatial_expression_matrix(sem=sem, genes_name_keep=genes_name_keep, 
                                                                        genes_remove_list=genes_remove_list)
    
    index_genes = get_idx_ref_genes_in_scDGE(dge=dge, genes_name_keep=genes_name_keep, sel_genes=sel_genes)
    
    topcells = get_best_cells(dge=dge, index_genes=index_genes, insitu_matrix=insitu_matrix,
                              method_dist=method_dist, top_sccells=top_sccells, enrich=enrich)
    
    dge_hvg = subset_dge_to_informative_features(dge=dge, topcells=topcells, sel_genes=sel_genes, 
                                                 genes_name_keep=genes_name_keep, top_cor_genes=top_cor_genes, 
                                                 pre_mode=pre_mode)
    
    cost_expression, cost_locations = get_expression_location_cost(dge_hvg=dge_hvg, coord=coord, 
                                                                   ns=num_neighbors_source, nt=num_neighbors_target)
    
    cost_marker_genes = get_marker_gene_cost(dge=dge, insitu_matrix=insitu_matrix, pre_mode=pre_mode, 
                                         tr=transform, index_genes=index_genes, method_dist=method_dist, topcells=topcells)
    
    p_locations, p_expression = get_distributions_over_expression_location(dge_hvg=dge_hvg, coord=coord)
    
    return(cost_expression, cost_locations, cost_marker_genes, p_locations, p_expression, topcells, genes_name_keep, dge, coord)

def predict_3D_gene_expression(
                              min_sccells_gene_expressed, keep_genes, genes_remove_list,
                              method_dist, top_sccells, enrich, top_cor_genes,
                              num_neighbors_target, num_neighbors_source,
                              transform, alpha, epsilon, max_iter, tol, verbose,
                              dge=None, gene_names=None, dge_path=None, sem=None, sem_path=None):
    
    # Load data if not provided
    if dge is None:
        dge, gene_names, cell_ids = load_scDGE(dge_path)
    if sem is None:
        sem = load_spatial_expression_matrix(coord_path=sem_path)
        
    # prepare data for OT-based reconstruction
    cost_expression, cost_locations, cost_marker_genes, p_locations, p_expression, topcells, genes_name_keep, dge, coord = prepare_reconstruction(
                        pre_mode="expr_pre",
                        dge=dge, gene_names=gene_names, sem=sem, 
                            min_sccells_gene_expressed=min_sccells_gene_expressed, keep_genes=keep_genes, genes_remove_list=genes_remove_list,
                              method_dist=method_dist, top_sccells=top_sccells, enrich=enrich,
                              num_neighbors_target=num_neighbors_target, num_neighbors_source=num_neighbors_source, top_cor_genes=top_cor_genes,
                              transform=transform
                                )
    
    # predict cell-to-cell mapping
    gw = predict_cell_to_cell_mappings(cost_marker_genes=cost_marker_genes, cost_expression=cost_expression, 
                              cost_locations=cost_locations, p_expression=p_expression, p_locations=p_locations, 
                            alpha=alpha, epsilon=epsilon, max_iter=max_iter, tol=tol, verbose=verbose)
    
    # reconstruct spatial gene expression
    sdge = predict_spatial_gene_expression(gw=gw, dge=dge, topcells=topcells, pre_mode="expr_pre")
    
    # prepare data for output
    sdge = np.transpose(sdge)
    sdge = np.concatenate((coord, sdge), axis=1)
    col_names = np.concatenate((np.array(["x", "y", "z"]), genes_name_keep), axis=0)
    sdge = round(pd.DataFrame(sdge, columns=col_names), 3)
    
    return sdge

def predict_cell_mappings(
                              min_sccells_gene_expressed, keep_genes, genes_remove_list,
                              method_dist, top_sccells, enrich, top_cor_genes,
                              num_neighbors_target, num_neighbors_source,
                              transform, alpha, epsilon, max_iter, tol, verbose,
                              dge=None, gene_names=None, cell_ids=None, dge_path=None, sem=None, sem_path=None):
    
    pre_mode="cell_mapping"
    
    # Load data if not provided
    if dge is None:
        dge, gene_names, cell_ids = load_scDGE(dge_path)
    if sem is None:
        sem = load_spatial_expression_matrix(coord_path=sem_path)
        
    # prepare data for OT-based reconstruction
    cost_expression, cost_locations, cost_marker_genes, p_locations, p_expression, topcells, genes_name_keep, dge, coord = prepare_reconstruction(
                        pre_mode=pre_mode,
                        dge=dge, gene_names=gene_names, sem=sem, 
                            min_sccells_gene_expressed=min_sccells_gene_expressed, keep_genes=keep_genes, genes_remove_list=genes_remove_list,
                              method_dist=method_dist, top_sccells=top_sccells, enrich=enrich,
                              num_neighbors_target=num_neighbors_target, num_neighbors_source=num_neighbors_source, top_cor_genes=top_cor_genes,
                              transform=transform
                                )
    
    # predict cell-to-cell mapping
    gw = predict_cell_to_cell_mappings(cost_marker_genes=cost_marker_genes, cost_expression=cost_expression, 
                              cost_locations=cost_locations, p_expression=p_expression, p_locations=p_locations, 
                            alpha=alpha, epsilon=epsilon, max_iter=max_iter, tol=tol, verbose=verbose)
    
    # prepare cell-to-cell mapping matrix
    # cell_ids_sel = cell_ids[topcells]
    cell_ids_sel = cell_ids
    gw = round(pd.DataFrame(gw, index=cell_ids_sel), 3)
    
    return gw