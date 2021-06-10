#!/usr/bin/env python3.8

# import modules
import pandas as pd
import sys
import argparse as ap

# read arguments
p = ap.ArgumentParser()
p.add_argument("--bac-summary", help="Summary file produced by gtdb-tk for the 120 bacterial markers", required=True)
p.add_argument("--bin-sizes", help="File containing MAG names and their total sizes in two columns", required=True)
p.add_argument("--output-file", help="Output table containing taxonomical annotations", required=True)
args = p.parse_args()

# read bac summary
Bac_summary = pd.read_csv(args.bac_summary, sep="\t", header="infer", index_col=0)
Bac_summary = Bac_summary.loc[ : , ["fastani_reference", "msa_percent", "classification"] ]
Bac_summary.classification = Bac_summary.classification.str.replace(" ", "_")
Bac_summary.index = Bac_summary.index.str.split("/").str[-1]
Bac_summary.index = Bac_summary.index.rename("MAG")

# read bin sizes
Bin_sizes = pd.read_csv(args.bin_sizes, sep="\t", header="infer", index_col=0)
Bin_sizes.index = Bin_sizes.index.str.split("/").str[-1]
Bin_sizes.index = Bin_sizes.index.str.replace(".fa", "")
Bin_sizes.index = Bin_sizes.index.rename("MAG")

# merge dataframes
df = Bac_summary.merge(right=Bin_sizes, how="outer", on="MAG")
df = df.loc[ : , ["Size(Mbp)", "msa_percent", "fastani_reference", "classification"] ]
df = df.sort_values(by=["Size(Mbp)", "msa_percent"], ascending=False)

# create new columns for each classification
df["species"] = df["classification"].str.split(";").str[-1].str.replace("s__", "")
df["genus"] = df["classification"].str.split(";").str[-2].str.replace("g__", "")
df["family"] = df["classification"].str.split(";").str[-3].str.replace("f__", "")
df["order"] = df["classification"].str.split(";").str[-4].str.replace("o__", "")
df["class"] = df["classification"].str.split(";").str[-5].str.replace("c__", "")
df["phylum"] = df["classification"].str.split(";").str[-6].str.replace("p__", "")
df["domain"] = df["classification"].str.split(";").str[-7].str.replace("d__", "")

# fill empty
df = df.fillna("---")
df = df.replace("", "---")

# reorder
df = df.loc[ : , ["Size(Mbp)", "msa_percent", "fastani_reference", "species", "genus", "family", "order", "class", "phylum", "domain"] ]
df = df.reset_index()

# write to output
df.to_csv(args.output_file, sep="\t", index=False, header=True, na_rep="NA")
