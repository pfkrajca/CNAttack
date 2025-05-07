def filter_and_count_zscores(
    zscore_array,
    upper_thresh: float = 1.0,
    lower_thresh: float = -1.0,
    min_cells: int = 600,
    adata: AnnData = None,
    obsm_key: str = None
) -> Tuple[pd.DataFrame, pd.Series, pd.Series]:

    #Validate dimensions
    if adata is not None and zscore_array.shape[0] != adata.n_obs:
        raise ValueError(f"zscore_array has {zscore_array.shape[0]} cells, but adata has {adata.n_obs}")

    #Convert to DataFrame with proper gene names
    zscore_df = pd.DataFrame(
        zscore_array,
        index=adata.obs_names if adata else None,
        columns=adata.var_names[:zscore_array.shape[1]] if adata else None
    )

    #Create significant msaks
    sig_mask = (zscore_df > upper_thresh) | (zscore_df < lower_thresh)
    sig_counts = sig_mask.sum(axis=0)

    #Filter genes
    keep_genes = sig_counts >= min_cells
    filtered_df = zscore_df.where(sig_mask & keep_genes, None)

    #Store in annData if requests
    if adata is not None and obsm_key is not None:
        filtered_matrix = np.zeros((adata.n_obs, adata.n_vars))
        filtered_matrix[:, :zscore_array.shape[1]] = filtered_df.values
        adata.obsm[obsm_key] = filtered_matrix
        adata.uns[f"{obsm_key}_genes"] = filtered_df.columns.tolist()

    #Count significant events
    up_counts = (zscore_df > upper_thresh).sum(axis=0)
    down_counts = (zscore_df < lower_thresh).sum(axis=0)

    print(f"Filtered to {keep_genes.sum()} genes with changes in ≥{min_cells} cells")
    return filtered_df.loc[:, keep_genes], up_counts, down_counts

#Using argparse to set argument for GitHub use
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter z-scores and count significantly changed genes.")
    parser.add_argument("h5ad_file", type=str, help="Input AnnData .h5ad file with z-scores in obsm")
    parser.add_argument("--zscore_key", type=str, default="delta_cluster_z", help="Key in obsm with z-score matrix")
    parser.add_argument("--obsm_key_out", type=str, default="filtered_zscores", help="Key to store filtered results in obsm")
    parser.add_argument("--upper_thresh", type=float, default=1.0, help="Upper z-score threshold")
    parser.add_argument("--lower_thresh", type=float, default=-1.0, help="Lower z-score threshold")
    parser.add_argument("--min_cells", type=int, default=600, help="Minimum cells with significant change")
    parser.add_argument("--output", type=str, default="filtered.h5ad", help="Output AnnData file")

    args = parser.parse_args()

    #Load the annData object
    adata = sc.read_h5ad(args.h5ad_file)

    if args.zscore_key not in adata.obsm:
        raise ValueError(f"Key '{args.zscore_key}' not found in adata.obsm.")

    #Call and use the function
    filtered_df, up_counts, down_counts = filter_and_count_zscores(
        zscore_array=adata.obsm[args.zscore_key],
        upper_thresh=args.upper_thresh,
        lower_thresh=args.lower_thresh,
        min_cells=args.min_cells,
        adata=adata,
        obsm_key=args.obsm_key_out
    )
