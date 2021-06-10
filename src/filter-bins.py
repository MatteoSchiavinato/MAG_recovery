#!/usr/bin/env python3

# modules
import sys
import pandas as pd
import argparse as ap
import os
import glob
from shutil import copyfile

# parser
p = ap.ArgumentParser()
p.add_argument("--bins-directory", help="Folder containing the original bins that were run through CheckM merge", required=True)
p.add_argument("--extension", default="fna", help="File extension to search for in the --bins-directory")
p.add_argument("--lineage-file", help="Output of CheckM merge, usually called \"merger.tsv\"", required=True)
p.add_argument("--max-cont", help="Maxium contamination detected in bin", default=5.0, type=float)
p.add_argument("--min-comp", help="Minimum completeness detected in bin", default=70.0, type=float)
p.add_argument("--max-het", help="Maximum heterogeneity of bin to retrieve", default=20.0, type=float)
p.add_argument("--output-directory", help="Directory where to place the merged bins, if non-existent will be created", required=True)
args = p.parse_args()

# create workspace
os.makedirs(args.output_directory, exist_ok=True)
out_bin_dir = args.output_directory + "/" + "bins"
os.makedirs(out_bin_dir, exist_ok=True)

# read lineage file
df = pd.read_csv(args.lineage_file, header="infer", index_col=False, sep="\t")

# generate report for unfiltered file
groups = df.loc[:,["Marker_lineage", "Completeness", "Contamination", "Strain_heterogeneity"]].groupby("Marker_lineage").agg(["count", "mean", "std"]).round(2)
groups = groups.reset_index()
groups.to_csv(f"{args.output_directory}/marker_lineage_counts.raw.tsv", sep="\t", index=False, na_rep="NaN")

# filter df based on criteria
mask_1 = df["Contamination"] <= args.max_cont
mask_2 = df["Completeness"] >= args.min_comp
mask_3 = df["Strain_heterogeneity"] <= args.max_het
sys.stderr.write(f"Bins with contamination <= {args.max_cont}: {mask_1.sum()}\n")
sys.stderr.write(f"Bins with completeness >= {args.min_comp}: {mask_2.sum()}\n")
sys.stderr.write(f"Bins with heterogeneity <= {args.max_het}: {mask_3.sum()}\n")
df = df.loc[mask_1 & mask_2 & mask_3, : ]
sys.stderr.write(f"Retained bins: {df.shape[0]}\n")

# generate report for filtered file
groups = df.loc[:,["Marker_lineage", "Completeness", "Contamination", "Strain_heterogeneity"]].groupby("Marker_lineage").agg(["count", "mean", "std"]).round(2)
groups = groups.reset_index()
groups.to_csv(f"{args.output_directory}/marker_lineage_counts.filtered.tsv", sep="\t", index=False, na_rep="NaN")

# select only bins that are in filtered df
Files = [ (x.split("/")[-1].rstrip(".fa"), x) for x in glob.glob(args.bins_directory + "/" + "*" + args.extension) ]
Valid_bins = df["Bin_Id"].to_list()
Files = [ (x[0], x[1]) for x in Files if x[0] in Valid_bins ]

counter = 0
total = len(Files)
for bin_name, bin_path in Files:
    src = bin_path
    dst = out_bin_dir + "/" + bin_name + ".fa"
    copyfile(src, dst)
    counter += 1
    sys.stderr.write(f"Copying in progress: {counter} / {total}\r")

sys.stderr.write("\n")
