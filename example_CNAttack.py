#!/usr/bin/env python
"""
CNAttack - A CNA detection choose your own adventure!

This example script demonstrates the full CNAttack workflow.

Usage:
    python example_CNAttack.py input.h5ad --output output_prefix [options]

"""
pip install -r requirements.txt

# Import CNAttack modules
from CNAttack.plot_cluster import plot_cluster_means_by_genomic_position
from CNAttack.smooth_profiles import compute_smoothed_profiles_from_adata
from CNAttack.zscores_compute import compute_all_cell_zscores_to_adata_optimized
from CNAttack.filter_zscores import filter_and_count_zscores
from CNAttack.hmm_for_CNAs import detect_cnvs_with_hmm_final
from CNAttack.format_CNVs import format_detected_cnvs_with_cell_counts


def main():
    #Parse command line arguments
    parser = argparse.ArgumentParser(
        description="CNAttack - Complete Copy Number Variation Analysis from scRNA-seq data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    #Required arguments
    parser.add_argument("input_file", type=str, help="Path to the AnnData (.h5ad) file")
    parser.add_argument("--output", type=str, default="cnattack_output", 
                        help="Output prefix for results")
    
    #Data parameters
    parser.add_argument("--layer", type=str, default=None, 
                        help="Layer in AnnData to use (e.g., 'counts', 'normalized')")
    parser.add_argument("--group_by", type=str, default="cell_type", 
                        help="Column in obs to group cells by")
    
    #Plotting options
    parser.add_argument("--plot", action="store_true", help="Generate genomic position plots")
    parser.add_argument("--downsample", type=int, default=100, 
                        help="Downsample factor for plotting")
    
    #Smoothing options
    parser.add_argument("--window_size", type=int, default=25, 
                        help="Window size for smoothing gene expression")
    parser.add_argument("--ref_group", type=str, default="Monocyte", 
                        help="Reference cell type/group for comparisons")
    
    #Z-score filtering options
    parser.add_argument("--upper_thresh", type=float, default=1.0, 
                        help="Upper z-score threshold")
    parser.add_argument("--lower_thresh", type=float, default=-1.0, 
                        help="Lower z-score threshold")
    parser.add_argument("--min_cells", type=int, default=600, 
                        help="Minimum cells with significant change")
    
    #HMM options
    parser.add_argument("--n_components", type=int, default=3, 
                        help="Number of hidden states for HMM")
    parser.add_argument("--n_iter", type=int, default=50, 
                        help="Number of iterations for HMM fitting")
    parser.add_argument("--min_non_nan", type=int, default=20, 
                        help="Minimum number of non-NaN cells for HMM")
    parser.add_argument("--chunk_size", type=int, default=100, 
                        help="Chunk size for processing features in HMM")
    
    #CNV formatting options
    parser.add_argument("--max_gap", type=float, default=1e6, 
                        help="Maximum gap for merging overlapping CNVs")
    parser.add_argument("--min_region_size", type=int, default=1000, 
                        help="Minimum region size to consider as valid CNV")
    parser.add_argument("--z_threshold", type=float, default=1.5, 
                        help="Z-score threshold for determining affected cells")
    
    #Workflow control
    parser.add_argument("--start_step", type=int, default=1, choices=range(1, 7),
                        help="Start workflow from step (1-6)")
    parser.add_argument("--end_step", type=int, default=6, choices=range(1, 7),
                        help="End workflow at step (1-6)")
    
    args = parser.parse_args()
    
    #Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output) if os.path.dirname(args.output) else '.', exist_ok=True)
    
    #Define intermediate file paths [checks work]
    smoothed_file = f"{args.output}_smoothed.h5ad"
    zscores_file = f"{args.output}_zscores.h5ad"
    filtered_file = f"{args.output}_filtered.h5ad"
    hmm_file = f"{args.output}_hmm.h5ad"
    cnv_file = f"{args.output}_cnvs.h5ad"
    
    print(f"CNAttack - A CNA detection choose your own adventure!")
    print(f"Input: {args.input_file}")
    print(f"Output prefix: {args.output}")
    
    #Load input data
    print(f"\nLoading input data from {args.input_file}...")
    adata = sc.read_h5ad(args.input_file)
    print(f"Data loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")
    
    #Execute workflow steps
    if 1 >= args.start_step and 1 <= args.end_step:
        print(f"\n Plotting mean expression across genomic positions")
        if args.plot:
            print(f"Generating genomic position plots...")
            plot_cluster_means_by_genomic_position(
                adata, 
                layer=args.layer, 
                group_by=args.group_by, 
                downsample=args.downsample
            )
            print(f"Plots generated.")
        else:
            print(f"Skipping plot generation (use --plot to enable)")
    
    if 2 >= args.start_step and 2 <= args.end_step:
        print(f"\n Computing smoothed expression profiles")
        print(f"Window size: {args.window_size}, Reference group: {args.ref_group}")
        smoothed_expr, celltype_profiles, global_avg, ref_avg = compute_smoothed_profiles_from_adata(
            adata,
            use_layer=args.layer,
            group_by=args.group_by,
            window_size=args.window_size,
            ref=args.ref_group
        )
        print(f"Smoothed profiles computed.")
        print(f"Shape of smoothed expression matrix: {smoothed_expr.shape}")
        print(f"Number of cell types found: {len(celltype_profiles)}")
        adata.uns['smoothed_genes'] = smoothed_expr.columns.tolist()
        adata.write(smoothed_file)
        print(f"Saved to {smoothed_file}")
    
    if 3 >= args.start_step and 3 <= args.end_step:
        print(f"\n Computing z-scores")
        # Load data if starting from this step
        if args.start_step > 2:
            print(f"Loading from {smoothed_file}...")
            adata = sc.read_h5ad(smoothed_file)
            smoothed_expr = None  # Will need to recompute from step 2
            celltype_profiles = None
            global_avg = None
            ref_avg = None
            
            # Recompute if needed
            if smoothed_expr is None:
                print(f"Recomputing smoothed profiles...")
                smoothed_expr, celltype_profiles, global_avg, ref_avg = compute_smoothed_profiles_from_adata(
                    adata,
                    use_layer=args.layer,
                    group_by=args.group_by,
                    window_size=args.window_size,
                    ref=args.ref_group
                )
        
        print(f"Computing z-scores...")
        adata = compute_all_cell_zscores_to_adata_optimized(
            smoothed=smoothed_expr,
            celltype_avg=celltype_profiles,
            global_avg=global_avg,
            cluster_avg=ref_avg,
            adata=adata,
            group_by=args.group_by
        )
        print(f"Z-scores computed. Results stored in:")
        print(f"  - adata.obsm['delta_global_z']")
        print(f"  - adata.obsm['delta_celltype_z']")
        print(f"  - adata.obsm['delta_cluster_z']")
        adata.write(zscores_file)
        print(f"Saved to {zscores_file}")
    
    if 4 >= args.start_step and 4 <= args.end_step:
        print(f"\n Filtering z-scores")
        # Load data if starting from this step
        if args.start_step > 3:
            print(f"Loading from {zscores_file}...")
            adata = sc.read_h5ad(zscores_file)
        
        print(f"Filtering z-scores (threshold: {args.upper_thresh}/{args.lower_thresh}, min cells: {args.min_cells})...")
        filtered_df, up_counts, down_counts = filter_and_count_zscores(
            zscore_array=adata.obsm['delta_cluster_z'],
            upper_thresh=args.upper_thresh,
            lower_thresh=args.lower_thresh,
            min_cells=args.min_cells,
            adata=adata,
            obsm_key='filtered_z'
        )
        print(f"Filtered to {filtered_df.shape[1]} significant genes.")
        adata.write(filtered_file)
        print(f"Saved to {filtered_file}")
    
    if 5 >= args.start_step and 5 <= args.end_step:
        print(f"\n Detecting CNVs using HMM")
        # Load data if starting from this step
        if args.start_step > 4:
            print(f"Loading from {filtered_file}...")
            adata = sc.read_h5ad(filtered_file)
        
        print(f"Running HMM with {args.n_components} states...")
        adata = detect_cnvs_with_hmm_final(
            adata=adata,
            matrix_name='filtered_z',
            n_components=args.n_components,
            n_iter=args.n_iter,
            random_state=42,
            output_prefix='hmm_cnv',
            min_non_nan=args.min_non_nan,
            chunk_size=args.chunk_size
        )
        print(f"HMM states detected and stored in adata.var['hmm_cnv_state']")
        adata.write(hmm_file)
        print(f"Saved to {hmm_file}")
    
    if 6 >= args.start_step and 6 <= args.end_step:
        print(f"\n Formatting and annotating CNVs")
        # Load data if starting from this step
        if args.start_step > 5:
            print(f"Loading from {hmm_file}...")
            adata = sc.read_h5ad(hmm_file)
        
        print(f"Formatting CNVs...")
        adata, cnv_df = format_detected_cnvs_with_cell_counts(
            adata=adata,
            state_col='hmm_cnv_state',
            output_col='detected_cnvs',
            chrom_col='chromosome',
            start_col='start',
            end_col='end',
            group_by_col=args.group_by,
            max_gap=args.max_gap,
            min_region_size=args.min_region_size,
            z_threshold=args.z_threshold
        )
        
        # Save CNV results
        adata.write(cnv_file)
        cnv_csv_file = f"{args.output}_cnv_regions.csv"
        cnv_df.to_csv(cnv_csv_file)
        print(f"CNV annotations stored in adata.obs['detected_cnvs']")
        print(f"Full AnnData saved to: {cnv_file}")
        print(f"CNV regions saved to: {cnv_csv_file}")
    
    print(f"Congratulations on choosing your own CNA detection adventure with CNAttack!")
    print(f"Final output file: {cnv_file}")


if __name__ == "__main__":
    main()
