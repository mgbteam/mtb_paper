#!/usr/bin/env python3

import argparse
import os

import pandas as pd

parser = argparse.ArgumentParser(
        description="Filter proteins based on FDRBench FDP estimates"
)
parser.add_argument(
        "-i", "--input", metavar=["TSV", "CSV"], required=True, nargs=2, action="append",
        help="Sample level protein results (tsv) and matching FDRBench output (csv)"
)
parser.add_argument(
        "-p", "--proteins", metavar="FILE", required=True,
        help="List of identified proteins on strain level in tsv format"
)
parser.add_argument(
        "-o", "--output", metavar="FILE", required=True,
        help="Output tsv file of proteins with added FDP filter column"
)
args = parser.parse_args()

filt_prots = None

for prot_file, fdp_file in args.input:
    df = pd.read_csv(fdp_file, sep=",")
    df.sort_values("combined_fdp", inplace=True)
    prev_row = {"q_value": 0, "Protein Probability": 1, "Top Peptide Probability": 1}

    for index, row in df.iterrows():
        if row["combined_fdp"] > 0.01:
            max_qval = prev_row["q_value"]
            break

        prev_row = row

    df = pd.read_csv(prot_file, sep="\t")
    filt_prots_sample = df[df["Protein Qvalue"] <= max_qval]

    if filt_prots is None:
        filt_prots = filt_prots_sample
    else:
        filt_prots = pd.concat([filt_prots, filt_prots_sample], axis=0)

df = pd.read_csv(args.proteins, sep="\t")
df["FDP Filter"] = df["protein"].isin(filt_prots["Protein"])
df = df[["protein", "FDP Filter"]]
df.to_csv(args.output, sep="\t", index=False)