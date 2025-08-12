
# Original Author: Abeer Abdulhakeem Mansour Alhasbary
# Source: https://github.com/Alhasbary/CTAPred
# Original Paper: https://www.sciencedirect.com/science/article/pii/S0010482524014367
# License: MIT License
# Date of Original File: 10/10/2024
# Date of this File: 8/1/2025
# This file has been included with modifications for parallelized processing for phytotherapy-ontology.



import argparse
import sys
import os
import pandas as pd
import numpy as np
import gc
from functools import partial
from multiprocessing import Pool
from SharedFunc import * 


def parse_args():
    parser = argparse.ArgumentParser(
        description="Predict potential protein target(s) by processing a large chemical file in parallel chunks.")

    parser.add_argument("--input", action="store", dest='inputFile', type=str, default='input', required=True,
                        help='Full path to the single input CSV file containing query compounds.')
    
    parser.add_argument("--output", action="store", dest='outputPath', type=str, default='output',
                        help='Full path to the directory where results are saved. Default is "output".')

    parser.add_argument("--fingerprint", action="store", default='ecfp', dest="fingerprint",
                        choices=['avalon', 'ecfp', 'fcfp', 'maccs'],
                        help="Desired fingerprint type. Default='ecfp'.")
    
    parser.add_argument("--nBits", action="store", type=int_range(8, 2048), default=2048, dest="nBits",
                        help="Length of the generated fingerprint (for avalon, ecfp, fcfp). Default=2048.")
    
    parser.add_argument("--radius", action="store", type=int, default=2, choices=[2, 3], dest="radius",
                        help="Radius for Morgan (ECFP/FCFP) fingerprints. Default=2.")

    parser.add_argument("--k", action="store", type=int, default=[1], dest="top_k", nargs='*',
                        help="Value(s) for 'top-k' reference compounds. Default=1.")
    
    parser.add_argument("--n_jobs", action="store", type=int, default=-1, dest="n_jobs",
                        help="Number of CPU cores to use. Default=-1 (use all available cores).")
    
    parser.add_argument("--verbose", action="store_true", default=False,
                        help="Enable this flag to print error messages during SMILES preprocessing.")

    args = parser.parse_args()
    return args


def process_chunk(chunk_df, args, CTA):
    """
    Processes a single DataFrame chunk: cleans SMILES, generates fingerprints, and conducts a similarity search.

    Args:
        chunk_df (DataFrame): A chunk of the input data with 'smiles' and 'np_id'.
        args (Namespace): The script's command-line arguments.
        CTA (DataFrame): The pre-loaded and pre-fingerprinted reference target dataset.

    Returns:
        DataFrame: A DataFrame with similarity search results for the chunk.
    """

    # Preprocess SMILES
    clean_query_structures = preprocess_Query_SMILES(None, args, chunk_df, n_jobs=1)
    if clean_query_structures.empty:
        return pd.DataFrame()
    clean_query_structures.rename(columns={'np_id': 'molregno'}, inplace=True)

    # Generate Fingerprints
    query_fps = gen_fps(None, args, clean_query_structures, 'Query Chunk', n_jobs=1)
    query_fps.rename(columns={'molregno': 'np_id'}, inplace=True)
    if query_fps.empty:
        return pd.DataFrame()

    # Similarity Search
    ChEMBL_similarity = Query_conductSS(None, args, query_fps[['np_id', 'fp']], CTA[['molregno', 'target_chembl_id', 'uniprot', 'fp']])

    return ChEMBL_similarity


def main():
    # Setup
    args = parse_args()

    os.makedirs(args.outputPath, exist_ok=True)
    logPath = os.path.join(args.outputPath, "predict_targets_log.txt")
    _logFile = open(logPath, 'a')
    log = Tee(sys.stdout, _logFile)

    print(" ".join(sys.argv), file=log)
    print(f"Processing input file: {args.inputFile}", file=log)

    if not os.path.exists(args.inputFile):
        print(f"Error: Input file not found at {args.inputFile}", file=log)
        exit()

    if args.fingerprint in ['avalon', 'maccs']:
        args.radius = None
        if args.fingerprint == 'maccs':
            args.nBits = 166

    # Prepare Reference Data
    print("\nLoading and preparing reference ChEMBL dataset", file=log)
    CTA, CTA_id = select_CTA_reference_dataset(log, args)
    CTA = CTA[['molregno', 'clean_smiles', 'uniprot', 'target_chembl_id']].drop_duplicates(ignore_index=True)
    CTA.dropna(inplace=True)

    # Generate Reference Fingerprints
    num_fp_jobs = os.cpu_count() if args.n_jobs == -1 else args.n_jobs
    print("num_fp_jobs", num_fp_jobs)
    
    CTA_fps = gen_fps(log, args, CTA[['molregno', 'clean_smiles']], 'Reference', n_jobs=num_fp_jobs)
    CTA = pd.merge(CTA, CTA_fps, on='molregno', how='inner')
    CTA.dropna(inplace=True)
    print(f"Reference dataset ready with {len(CTA)} compound-target pairs", file=log)

    # Split Input
    query_list = pd.read_csv(args.inputFile)
    print(f"\nLoaded {len(query_list)} compounds from input file", file=log)

    df_chunks = np.array_split(query_list, num_fp_jobs)
    print(f"Input data split into {len(df_chunks)} chunks for parallel processing.", file=log)

    worker_func = partial(process_chunk, args=args, CTA=CTA)

    with Pool(processes=num_fp_jobs) as pool:
        list_of_results = pool.map(worker_func, df_chunks)

    # Combine Results
    print("\nCombining results", file=log)
    ChEMBL_similarity = pd.concat(list_of_results, ignore_index=True)
    del list_of_results
    gc.collect()

    # Rank
    if not ChEMBL_similarity.empty:
        print(f"Found {len(ChEMBL_similarity)} total similarity hits", file=log)
        final_output_dir = os.path.join(args.outputPath, "potential_targets/")
        os.makedirs(final_output_dir, exist_ok=True)


        base_name = os.path.splitext(os.path.basename(args.inputFile))[0]
        potential_targets_path = os.path.join(final_output_dir, f"{base_name}_potential_targets_k_{args.top_k}")

        rankTargets(log, args, ChEMBL_similarity, potential_targets_path)
    else:
        print("No potential targets identified", file=log)

    _logFile.close()



if __name__ == '__main__':
    main()
 