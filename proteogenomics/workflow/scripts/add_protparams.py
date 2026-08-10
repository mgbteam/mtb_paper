#!/usr/bin/env python3

import argparse
import csv
import os

parser = argparse.ArgumentParser(description="Join TSV files by common column")
parser.add_argument(
        "-t", "--table", required=True,
        help="Protein table in tsv format"
)
parser.add_argument(
        "-p", "--protparams", required=True,
        help="Protein parameters in tsv format"
)
parser.add_argument(
        "-o", "--outfile", required=True,
        help="Merged output table in tsv format"
)
parser.add_argument(
        "-T", "--table-column", default="protein",
        help="Table column to combine by (default: 'protein')"
)
parser.add_argument(
        "-P", "--protparam-column", default="Protein",
        help="Protparam column to combine by (default: 'Protein')"
)
parser.add_argument(
        "-d", "--decoy-prefix", default="rev_",
        help="Decoy prefix (default: rev_)"
)
parser.add_argument(
        "-e", "--entrapment_suffix", default="_p_target",
        help="Entrapment suffix (default: _p_target)"
)
args = parser.parse_args()

os.makedirs(os.path.dirname(args.outfile), exist_ok=True)

with open(args.protparams, "r") as fi:
    reader = csv.DictReader(fi, delimiter="\t")
    pp_cols = [col for col in reader.fieldnames if col != args.protparam_column]
    protparams = {}

    for row in reader:
        params = {k: v for k, v in row.items() if k != args.protparam_column}
        protparams[row[args.protparam_column]] = params

outrows = []

with open(args.table, "r") as fi:
    reader = csv.DictReader(fi, delimiter="\t")
    prot_cols = reader.fieldnames

    for row in reader:
        prot_id = row[args.table_column]
        prot_id = prot_id.removeprefix(args.decoy_prefix)
        prot_id = prot_id.removesuffix(args.entrapment_suffix)
        row.update(protparams[prot_id])
        outrows.append(row)

with open(args.outfile, "w") as fo:
    out_cols = prot_cols + pp_cols
    writer = csv.DictWriter(fo, delimiter="\t", fieldnames=out_cols)
    writer.writeheader()
    writer.writerows(outrows)
