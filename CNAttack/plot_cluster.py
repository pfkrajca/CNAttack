def plot_cluster_means_by_genomic_position(adata, layer="counts", group_by="cell_type", downsample=100):
    gene_indices = np.arange(0, adata.n_vars, downsample)

    #Get chromosome information
    chromosomes = adata.var["chromosome"].values[gene_indices] if "chromosome" in adata.var else None

    if chromosomes is not None:
        valid_chroms = [str(i) for i in range(1, 23)] + ['X', 'chrX']
        chrom_mask = np.isin(chromosomes, valid_chroms)
        chrom_changes = np.where((chromosomes[:-1] != chromosomes[1:]) &
                                 (chrom_mask[:-1] | chrom_mask[1:]))[0] + 1
        chrom_boundaries = [0] + chrom_changes.tolist() + [len(chromosomes)]
        chrom_midpoints = [(chrom_boundaries[i] + chrom_boundaries[i+1]) // 2
                           for i in range(len(chrom_boundaries)-1)]
        chrom_labels = [chromosomes[i] for i in chrom_boundaries[:-1]]

    #Calculate global y-axis limits
    all_means = []
    for group in adata.obs[group_by].unique():
        group_mask = adata.obs[group_by] == group
        group_data = adata[group_mask].layers[layer][:, gene_indices]
        if issparse(group_data):
            group_data = group_data.toarray()
        all_means.append(np.mean(group_data, axis=0))

    global_ymin = min(np.min(means) for means in all_means)
    global_ymax = max(np.max(means) for means in all_means)
    y_padding = (global_ymax - global_ymin) * 0.1

    #Combined plot
    plt.figure(figsize=(14, 6))
    for group, means in zip(adata.obs[group_by].unique(), all_means):
        plt.plot(means, label=f"{group}", linewidth=1)

    if chromosomes is not None:
        for x in chrom_boundaries[1:-1]:
            if chrom_mask[x]:
                plt.axvline(x, color='gray', linestyle='--', alpha=0.5, linewidth=0.5)
        plt.xticks(
            [mp for mp, lbl in zip(chrom_midpoints, chrom_labels) if lbl in valid_chroms],
            [lbl for lbl in chrom_labels if lbl in valid_chroms],
            rotation=90
        )

    plt.ylim(global_ymin - y_padding, global_ymax + y_padding)
    plt.title(f"All Groups (Downsampled {downsample}x)")
    plt.xlabel("Genomic Position")
    plt.ylabel("Mean Counts")
    plt.legend(bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.show()

    #Individual group plots
    for group, means in zip(adata.obs[group_by].unique(), all_means):
        plt.figure(figsize=(14, 4))
        plt.plot(means, color='steelblue', linewidth=1)

        if chromosomes is not None:
            for x in chrom_boundaries[1:-1]:
                if chrom_mask[x]:
                    plt.axvline(x, color='gray', linestyle='--', alpha=0.5, linewidth=0.5)
            plt.xticks(
                [mp for mp, lbl in zip(chrom_midpoints, chrom_labels) if lbl in valid_chroms],
                [lbl for lbl in chrom_labels if lbl in valid_chroms],
                rotation=90
            )

        plt.ylim(global_ymin - y_padding, global_ymax + y_padding)
        plt.title(f"Group: {group}\n(Mean Expression, {downsample}x downsampled)")
        plt.xlabel("Genomic Position")
        plt.ylabel("Mean Counts")
        plt.tight_layout()
        plt.show()

#Using argparse to set argument for GitHub use
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plots mean expression across genomic positions for each group in annData.")
    parser.add_argument("h5ad_file", type=str, help="Path to the AnnData .h5ad file.")
    parser.add_argument("--layer", type=str, default="counts", help="Layer name to use (default: counts).")
    parser.add_argument("--group_by", type=str, default="cell_type", help="Column in obs to group cells by (default: cell_type).")
    parser.add_argument("--downsample", type=int, default=100, help="Downsample factor for plotting, plot every nth gene to reduce visual clutter (default: 100).")

    args = parser.parse_args()

    #Load the annData object
    adata = sc.read_h5ad(args.h5ad_file)

    #Call and use the function
    plot_cluster_means_by_genomic_position(
        adata,
        layer=args.layer,
        group_by=args.group_by,
        downsample=args.downsample
    )
