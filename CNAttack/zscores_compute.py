def compute_all_cell_zscores_to_adata_optimized(smoothed: pd.DataFrame,
                                                celltype_avg: dict,
                                                global_avg: pd.Series,
                                                cluster_avg: pd.Series,
                                                adata: AnnData,
                                                group_by: str = 'cell_type') -> AnnData:
    cells = smoothed.index
    genes = smoothed.columns

    #Precompute arrays for faster access
    smoothed_np = smoothed.values
    global_avg_np = global_avg.values
    cluster_avg_np = cluster_avg.values

    #Build a matrix of per-cell-type averages aligned to cells
    group = adata.obs.loc[cells, group_by].values
    celltype_avg_np = np.stack([celltype_avg[ct].loc[genes].values for ct in group])

    #Helper function that calculates z-svore along gene axis
    def zscore_rowwise(matrix):
        mean = matrix.mean(axis=1, keepdims=True)
        std = matrix.std(axis=1, keepdims=True)
        std[std == 0] = 1
        return (matrix - mean) / std

    #Allocate arrays
    delta_global = smoothed_np - global_avg_np
    delta_celltype = smoothed_np - celltype_avg_np
    delta_cluster = smoothed_np - cluster_avg_np

    #Store results in adata.obs
    adata.obsm['delta_global_z'] = zscore_rowwise(delta_global)
    adata.obsm['delta_celltype_z'] = zscore_rowwise(delta_celltype)
    adata.obsm['delta_cluster_z'] = zscore_rowwise(delta_cluster)
    adata.uns['delta_zscore_genes'] = list(genes)

    return adata

#Using argparse to set argument for GitHub use
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute smoothed profiles and z-score deltas from an AnnData file.")
    parser.add_argument("h5ad_file", type=str, help="Path to the .h5ad file.")
    parser.add_argument("--layer", type=str, default=None, help="Layer name to use for expression.")
    parser.add_argument("--group_by", type=str, default="cell_type", help="Column in .obs to group cells.")
    parser.add_argument("--window_size", type=int, default=25, help="Window size for smoothing.")
    parser.add_argument("--ref", type=str, default="Monocyte", help="Reference cell type.")
    parser.add_argument("--out_prefix", type=str, default="output", help="Prefix for output files.")
    parser.add_argument("--compute_z", action="store_true", help="Whether to compute delta z-scores and store in AnnData.")

    args = parser.parse_args()

    #Load the annData object
    adata = sc.read_h5ad(args.h5ad_file)

    #Call and use the function
    smoothed_expr, celltype_profiles, global_avg, ref_avg = compute_smoothed_profiles_from_adata(
        adata,
        use_layer=args.layer,
        group_by=args.group_by,
        window_size=args.window_size,
        ref=args.ref
    )
