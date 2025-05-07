#!/bin/bash
# CNV Detection Pipeline Usage Example
# This script demonstrates how to use the CNV detection functions in sequence

# Path to your input data
INPUT_FILE="your_data.h5ad"

# Step 1: Order genes by genomic position and create smoothed expression profiles
# This step prepares the data by organizing genes by chromosome position and smoothing expression
echo "Step 1: Creating smoothed expression profiles..."
python3 smooth_profiles.py "$INPUT_FILE" \
    --layer "counts" \
    --group_by "cell_type" \
    --window_size 25 \
    --ref "Monocyte" \
    --out_prefix "smoothed_data"

# Step 2: Compute Z-scores relative to global, cell type, and reference averages
# This identifies cells with unusual expression patterns
echo "Step 2: Computing Z-scores..."
python3 compute_zscores.py "$INPUT_FILE" \
    --layer "counts" \
    --group_by "cell_type" \
    --window_size 25 \
    --ref "Monocyte" \
    --compute_z \
    --out_prefix "zscores"

# Step 3: Filter Z-scores to focus on significant gene expression changes
# This reduces noise and focuses on regions with strong signals
echo "Step 3: Filtering Z-scores..."
python3 filter_zscores.py "$INPUT_FILE" \
    --zscore_key "delta_cluster_z" \
    --obsm_key_out "filtered_z" \
    --upper_thresh 1.0 \
    --lower_thresh -1.0 \
    --min_cells 600 \
    --output "filtered_data.h5ad"

# Step 4: Detect CNVs using Hidden Markov Model
# This identifies patterns consistent with copy number variations
echo "Step 4: Detecting CNVs with HMM..."
python3 detect_cnvs.py "filtered_data.h5ad" \
    --matrix_name "filtered_z" \
    --n_components 3 \
    --n_iter 50 \
    --output_prefix "hmm_cnv" \
    --min_non_nan 20 \
    --chunk_size 100 \
    --output "hmm_results.h5ad"

# Step 5: Format and summarize the detected CNVs
# This provides a final report of CNV regions and affected cells
echo "Step 5: Formatting detected CNVs..."
python3 format_cnvs.py "hmm_results.h5ad" \
    --state_col "hmm_cnv_state" \
    --output_col "detected_cnvs" \
    --chrom_col "chromosome" \
    --start_col "start" \
    --end_col "end" \
    --group_by_col "cell_type" \
    --max_gap 1000000 \
    --min_region_size 1000 \
    --z_threshold 1.5 \
    --output "final_cnv_results.h5ad"

# Optional visualization step
echo "Step 6: Visualizing genomic expression profiles..."
python3 plot_genomic_profiles.py "$INPUT_FILE" \
    --layer "counts" \
    --group_by "cell_type" \
    --downsample 100

echo "CNV detection pipeline completed!"
