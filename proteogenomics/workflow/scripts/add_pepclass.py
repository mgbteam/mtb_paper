#!/usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser(
        description="Add PeptideClassifier information"
)
parser.add_argument(
        "-t", "--table", required=True,
        help="Protein table in TSV format"
)
parser.add_argument(
        "-p", "--pepclass", required=True,
        help="PeptideClassifier results"
)
parser.add_argument(
        "-o", "--outfile", required=True,
        help="Output file in TSV format"
)
args = parser.parse_args()

pepclasses = ["unique", "ambiguous", "contam"]
proteins = {}

with open(args.pepclass, "r") as fi:
    reader = csv.reader(fi, delimiter="\t")

    for row in reader:
        for protein in row[2].split(";"):
            if protein not in proteins:
                proteins[protein] = {c: 0 for c in pepclasses}

            proteins[protein][row[1]] += 1

with open(args.table, "r") as fi, open(args.outfile, "w") as fo:
    reader = csv.DictReader(fi, delimiter="\t")
    cols = reader.fieldnames + pepclasses

    writer = csv.DictWriter(fo, delimiter="\t", fieldnames=cols)
    writer.writeheader()

    for row in reader:
        if row["protein"] not in proteins:
            print(f"WARNING: no classified peptide for {row['protein']}")
            row.update({pc: 0 for pc in pepclasses})
        else:
            row.update(proteins[row["protein"]])

        writer.writerow(row)
