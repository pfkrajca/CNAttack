#Orders the variables (genes) in the AnnData object by chromosome and start.
def order_genes_by_position(adata: AnnData) -> AnnData:
    gene_order = (
        adata.var
        .sort_values(['chromosome', 'start'])
        .index
    )
    return adata[:, gene_order]

#Applies a sliding window average across genes.
def smooth_expression_matrix(matrix: pd.DataFrame, window_size: int = 25) -> pd.DataFrame:
    smoothed = []
    gene_names = matrix.columns.to_list()

    for i in range(len(gene_names) - window_size + 1):
        window = gene_names[i:i + window_size]
        avg_expr = matrix[window].mean(axis=1)
        smoothed.append(avg_expr.rename(f"{window[0]}_to_{window[-1]}"))

    return pd.concat(smoothed, axis=1)

def compute_smoothed_profiles_from_adata(
    adata: AnnData,
    use_layer: str = None,
    group_by: str = 'cell_type',
    window_size: int = 25,
    ref: str = 'Monocyte'
) -> Tuple[pd.DataFrame, Dict[str, pd.Series], pd.Series, pd.Series]:

    #Order the gene
    adata_ordered = order_genes_by_position(adata)

    #Get expression matrix
    if use_layer:
        expr = pd.DataFrame(
            adata_ordered.layers[use_layer].toarray(),
            index=adata_ordered.obs_names,
            columns=adata_ordered.var_names
        )
    else:
        expr = pd.DataFrame(
            adata_ordered.X.toarray() if hasattr(adata_ordered.X, "toarray") else adata_ordered.X,
            index=adata_ordered.obs_names,
            columns=adata_ordered.var_names
        )

    #Smooth once for all cells
    smoothed_expr = smooth_expression_matrix(expr, window_size)

    #Global average
    global_avg = smoothed_expr.mean(axis=0)

    #Per-cell type acerages
    celltypes = adata_ordered.obs[group_by]
    celltype_profiles = {
        ct: smoothed_expr.loc[celltypes == ct].mean(axis=0)
        for ct in celltypes.unique()
    }

    #Reference cluster average
    if ref in celltype_profiles:
        cluster_avg = celltype_profiles[ref]
    else:
        raise ValueError(f"Reference cell type '{ref}' not found in adata.obs['{group_by}'].")

    return smoothed_expr, celltype_profiles, global_avg, cluster_avg

#Using argparse to set argument for GitHub use
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute smoothed expression profiles from an AnnData file.")
    parser.add_argument("h5ad_file", type=str, help="Path to the .h5ad file.")
    parser.add_argument("--layer", type=str, default=None, help="Layer name to use for expression.")
    parser.add_argument("--group_by", type=str, default="cell_type", help="Column in .obs to group cells.")
    parser.add_argument("--window_size", type=int, default=25, help="Window size for smoothing.")
    parser.add_argument("--ref", type=str, default="Monocyte", help="Reference cell type.")
    parser.add_argument("--out_prefix", type=str, default="output", help="Prefix for output CSV files.")

    args = parser.parse_args()

    #Load the annData object
    adata = sc.read_h5ad(args.h5ad_file)

    #Call and use the function
    smoothed_expr, celltype_profiles, global_avg, cluster_avg = compute_smoothed_profiles_from_adata(
        adata,
        use_layer=args.layer,
        group_by=args.group_by,
        window_size=args.window_size,
        ref=args.ref
    )
