def detect_cnvs_with_hmm_final(adata, matrix_name="filtered_z", n_components=3,
                             n_iter=50, random_state=42, output_prefix="hmm_cnv",
                             min_non_nan=20, chunk_size=100):
    """
    Final robust HMM implementation that:
    1. Uses pre-initialized model parameters
    2. Handles sparse data gracefully
    3. Provides detailed error reporting
    4. Stores per-gene HMM state and posterior in adata.var
    """

    if matrix_name not in adata.obsm:
        raise ValueError(f"Matrix '{matrix_name}' not found in adata.obsm")

    data = adata.obsm[matrix_name]
    n_cells, n_features = data.shape

    all_states = np.full(n_features, -1, dtype=np.int8)
    all_posteriors = np.full((n_features, n_components), np.nan, dtype=np.float16)

    print("Identifying valid features...")
    valid_features = [f for f in range(n_features) if np.sum(~np.isnan(data[:, f])) >= min_non_nan]

    if len(valid_features) < n_components:
        raise ValueError(f"Only {len(valid_features)} valid features found - need at least {n_components} for HMM")

    #Initialize HMM model
    model = hmm.GaussianHMM(
        n_components=n_components,
        covariance_type="diag",
        n_iter=0,
        init_params="",
        means_prior=np.linspace(-2, 2, n_components).reshape(-1, 1),
        covars_prior=np.ones((n_components, 1)),
        transmat_prior=np.full((n_components, n_components), 1/n_components),
        random_state=random_state
    )

    model_fitted = False

    print("Processing in chunks...")
    for chunk_start in tqdm(range(0, len(valid_features), chunk_size)):
        chunk_features = valid_features[chunk_start:chunk_start + chunk_size]
        chunk = data[:, chunk_features]

        valid_cells = ~np.isnan(chunk).any(axis=1)
        clean_chunk = chunk[valid_cells, :]

        if clean_chunk.shape[0] < min_non_nan:
            continue

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            scaled_chunk = StandardScaler().fit_transform(clean_chunk)

        if not model_fitted:
            try:
                model.n_iter = n_iter
                model.fit(scaled_chunk.T)
                model_fitted = True
                print("Model successfully fitted!")
            except Exception as e:
                print(f"Fit failed: {str(e)}")
                continue

        try:
            states = model.predict(scaled_chunk.T)
            posteriors = model.predict_proba(scaled_chunk.T)

            for i, f in enumerate(chunk_features):
                all_states[f] = states[i]
                all_posteriors[f] = posteriors[i].astype(np.float16)
        except Exception as e:
            print(f"Prediction failed for chunk {chunk_start}: {str(e)}")
            continue

    if not model_fitted:
        raise RuntimeError("Could not fit HMM to any data chunks")

    #Store in adata.var
    adata.var[f"{output_prefix}_state"] = -1
    for i in range(n_components):
        adata.var[f"{output_prefix}_prob_state_{i}"] = np.nan

    for idx, gene in enumerate(adata.var_names):
        adata.var.loc[gene, f"{output_prefix}_state"] = int(all_states[idx])
        for i in range(n_components):
            adata.var.loc[gene, f"{output_prefix}_prob_state_{i}"] = all_posteriors[idx, i]

    print("Stored HMM states and posteriors in adata.var")
    return adata

#Using argparse to set argument for GitHub use
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Detect CNVs using HMM on z-score data.")
    parser.add_argument("h5ad_file", type=str, help="Input AnnData .h5ad file with z-scores in obsm")
    parser.add_argument("--matrix_name", type=str, default="filtered_z", help="Key in obsm with z-score matrix")
    parser.add_argument("--n_components", type=int, default=3, help="Number of hidden states for HMM")
    parser.add_argument("--n_iter", type=int, default=50, help="Number of iterations for HMM fitting")
    parser.add_argument("--output_prefix", type=str, default="hmm_cnv", help="Prefix for output fields in adata.var")
    parser.add_argument("--min_non_nan", type=int, default=20, help="Minimum number of non-NaN cells for a feature")
    parser.add_argument("--chunk_size", type=int, default=100, help="Chunk size for processing features")
    parser.add_argument("--output", type=str, default="hmm_cnv_results.h5ad", help="Output AnnData file")

    args = parser.parse_args()

    #Load the annData object
    adata = sc.read_h5ad(args.h5ad_file)

    #Call and use the function
    adata = detect_cnvs_with_hmm_final(
        adata=adata,
        matrix_name=args.matrix_name,
        n_components=args.n_components,
        n_iter=args.n_iter,
        output_prefix=args.output_prefix,
        min_non_nan=args.min_non_nan,
        chunk_size=args.chunk_size
    )
