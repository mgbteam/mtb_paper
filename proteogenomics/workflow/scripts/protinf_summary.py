#!/usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser(
        description="Combine protein inference summary files"
)
parser.add_argument(
        "-i", "--input", required=True,
        nargs=2, metavar=("STRAIN", "FILE"), action="append",
        help="Strain and filter stats tsv (flag can be used multiple times)"
)
parser.add_argument(
        "-o", "--output", metavar="FILE", required=True,
        help="Output combined summary across all strains in tsv format"
)
args = parser.parse_args()

cols = ["strain"]
table = []

for strain, file in args.input:
    outrow = {"strain": strain}

    with open(file, "r") as fi:
        reader = csv.DictReader(fi, delimiter="\t")

        for row in reader:
            if row["groups"] not in cols:
                cols.append(row["groups"])

            outrow[row["groups"]] = int(row["proteins"])

    table.append(outrow)

with open(args.output, "w") as fo:
    writer = csv.DictWriter(fo, delimiter="\t", fieldnames=cols, restval=0)
    writer.writeheader()
    writer.writerows(table)
