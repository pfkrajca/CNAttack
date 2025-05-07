def format_detected_cnvs_with_cell_counts(
    adata,
    state_col='hmm_cnv_state',
    output_col='detected_cnvs',
    chrom_col='chromosome',
    start_col='start',
    end_col='end',
    group_by_col='cell_type',
    max_gap=1e6,
    min_region_size=1000,
    z_threshold=1.5
):
    """
    Formats detected CNVs:
    - Excludes neutral regions
    - Merges overlapping CNVs
    - Counts affected cells per region and cell type
    - Annotates CNV type (gain/loss)
    """

    #Only consider CNV states
    valid_states = [0, 2]
    valid_features = adata.var[state_col].isin(valid_states)
    sorted_vars = adata.var[valid_features].sort_values([chrom_col, start_col])

    if len(sorted_vars) == 0:
        print("No CNV regions found!")
        adata.obs[output_col] = "neutral"
        return adata, pd.DataFrame()

    #Map state to CN type
    state_to_cn = {0: '0', 2: '4'}
    cn_to_type = {'0': 'loss', '4': 'gain'}

    #Build regions list
    regions = []
    for _, row in sorted_vars.iterrows():
        cn = state_to_cn[row[state_col]]
        regions.append({
            'chrom': row[chrom_col],
            'start': float(row[start_col]),
            'end': float(row[end_col]),
            'cn': cn,
            'cnv_type': cn_to_type[cn],
            'genes': [row.name]
        })

    #Merge overlapping CNVs
    merged = []
    current = None
    for region in sorted(regions, key=lambda x: (x['chrom'], x['start'])):
        if current is None:
            current = region.copy()
            continue

        if (
            region['chrom'] == current['chrom']
            and region['cn'] == current['cn']
            and region['start'] <= current['end'] + max_gap
        ):
            current['end'] = max(current['end'], region['end'])
            current['genes'].extend(region['genes'])
        else:
            if current['end'] - current['start'] >= min_region_size:
                merged.append(current)
            current = region.copy()

    if current and current['end'] - current['start'] >= min_region_size:
        merged.append(current)

    if not merged:
        print("No valid CNV regions after merging.")
        adata.obs[output_col] = "neutral"
        return adata, pd.DataFrame()

    #Create CNV matrix and collect stats
    results = []
    cell_cnv_matrix = np.zeros((adata.n_obs, len(merged)), dtype=bool)

    for i, region in enumerate(merged):
        expr = adata[:, region['genes']].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            z_scores = zscore(expr, axis=0, nan_policy='omit')
            region_z = np.nanmean(z_scores, axis=1)

        #Determine affected cells
        if region['cn'] == '0':
            cnv_cells = region_z < -z_threshold
        elif region['cn'] == '4':
            cnv_cells = region_z > z_threshold
        else:
            cnv_cells = np.zeros(adata.n_obs, dtype=bool)

        cell_cnv_matrix[:, i] = cnv_cells
        n_cells = np.sum(cnv_cells)

        #Count by cell type
        if group_by_col in adata.obs.columns:
            cell_types = adata.obs.loc[cnv_cells, group_by_col]
            cell_type_counts = cell_types.value_counts().to_dict()
        else:
            cell_type_counts = {}

        region_str = f"{region['chrom']}:{int(region['start'])}-{int(region['end'])} (CN {region['cn']})"

        results.append({
            'region': region_str,
            'n_cells': n_cells,
            'genes': ','.join(region['genes']),
            'cnv_type': region['cnv_type'],
            'group_by_counts': cell_type_counts
        })

    #Final DataFrame
    cnv_df = pd.DataFrame(results).sort_values('n_cells', ascending=False)
    cnv_df['percent_cells'] = (cnv_df['n_cells'] / adata.n_obs * 100).round(1)

    #Assign CNVs to cells
    adata.obs[output_col] = "neutral"
    for i in range(len(merged)):
        region_str = cnv_df.iloc[i]['region']
        adata.obs.loc[cell_cnv_matrix[:, i], output_col] = region_str

    #Handle multiple CNVs per cell
    multi_cnv_cells = np.sum(cell_cnv_matrix, axis=1) > 1
    if np.any(multi_cnv_cells):
        for cell_idx in np.where(multi_cnv_cells)[0]:
            patterns = [
                cnv_df.iloc[i]['region']
                for i in np.where(cell_cnv_matrix[cell_idx])[0]
            ]
            adata.obs.iloc[cell_idx, adata.obs.columns.get_loc(output_col)] = "; ".join(patterns)

    print(f"Detected {len(merged)} CNV regions across {np.sum(cell_cnv_matrix)} cell-region pairs")
    print("\nTop CNV regions by cell count:")
    print(cnv_df)

    #Store extra info
    adata.uns['cnv_stats'] = {
        'total_regions': len(merged),
        'total_cell_cnv_pairs': np.sum(cell_cnv_matrix),
        'cells_with_cnvs': np.sum(adata.obs[output_col] != "neutral"),
        'cell_cnv_matrix': cell_cnv_matrix,
        'merged_regions': merged
    }

    return adata, cnv_df

#Using argparse to set argument for GitHub use
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format detected CNVs with cell counts and cell type information.")
    parser.add_argument("h5ad_file", type=str, help="Input AnnData .h5ad file with CNV states in adata.var")
    parser.add_argument("--state_col", type=str, default="hmm_cnv_state", help="Column in var with CNV state")
    parser.add_argument("--output_col", type=str, default="detected_cnvs", help="Output column in obs to store CNV regions")
    parser.add_argument("--chrom_col", type=str, default="chromosome", help="Column in var for chromosome names")
    parser.add_argument("--start_col", type=str, default="start", help="Column in var for CNV start position")
    parser.add_argument("--end_col", type=str, default="end", help="Column in var for CNV end position")
    parser.add_argument("--group_by_col", type=str, default="cell_type", help="Column in obs for grouping cells (e.g., cell type)")
    parser.add_argument("--max_gap", type=float, default=1e6, help="Maximum gap for merging overlapping CNVs")
    parser.add_argument("--min_region_size", type=int, default=1000, help="Minimum region size to consider as valid CNV")
    parser.add_argument("--z_threshold", type=float, default=1.5, help="Z-score threshold for determining affected cells")
    parser.add_argument("--output", type=str, default="formatted_cnvs.h5ad", help="Output AnnData file")

    args = parser.parse_args()

    #Load the annData object
    adata = sc.read_h5ad(args.h5ad_file)

    #Call and use the function
    adata, cnv_df = format_detected_cnvs_with_cell_counts(
        adata=adata,
        state_col=args.state_col,
        output_col=args.output_col,
        chrom_col=args.chrom_col,
        start_col=args.start_col,
        end_col=args.end_col,
        group_by_col=args.group_by_col,
        max_gap=args.max_gap,
        min_region_size=args.min_region_size,
        z_threshold=args.z_threshold
    )
