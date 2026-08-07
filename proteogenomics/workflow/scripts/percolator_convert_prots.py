#!/usr/bin/env python3

import argparse
import csv
import os

parser = argparse.ArgumentParser(
        description="Deduplicate and split percolator output"
)
parser.add_argument(
        "percolator",
        help="Percolator output file (target proteins)"
)
parser.add_argument(
        "outdir",
        help="Output folder where the subsets will be written"
)
parser.add_argument(
        "--refseq-prefix", default="REF",
        help="RefSeq ID prefix (default: REF)"
)
parser.add_argument(
        "--decoy-prefix", default="rev_",
        help="Decoy ID prefix (default: rev_)"
)
parser.add_argument(
        "--contam-prefix", default="sp|",
        help="Contaminant ID prefix (default: sp|)"
)
args = parser.parse_args()

new_header = [
    "protein",
    "score",
    "q_value",
    "posterior_error_prob",
]


def add_if_better(subset, protein, row):
    if protein not in subset or int(row[1]) < int(subset[protein][1]):
        subset[protein] = [protein] + row[1:4]


# split peptides into subsets and only keep the best scoring PSM
subsets = {"all": {}, "refseq": {}, "novels": {}, "contams": {}}

with open(args.percolator, "r") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)

    for row in reader:
        # only add first protein of each protein group (like FDRBench)
        prot = row[0].rstrip(";").split(";")[0]
        add_if_better(subsets["all"], prot, row)

        if prot.startswith(args.contam_prefix):
            add_if_better(subsets["contams"], prot, row)
        else:
            base = prot.split("|")[0]

            if prot.startswith(args.refseq_prefix):
                if not any([s in base for s in ["_-", "_+", "_p_", "_p|"]]):
                    add_if_better(subsets["refseq"], prot, row)
                else:
                    add_if_better(subsets["novels"], prot, row)
            else:
                add_if_better(subsets["novels"], prot, row)

# write subsets to files
os.makedirs(args.outdir, exist_ok=True)

for subset_name, subset_dict in subsets.items():
    with open(f"{args.outdir}/{subset_name}.tsv", "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(new_header)
        writer.writerows(subset_dict.values())
