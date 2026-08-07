#!/usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser(
        description="Filter ambiguous protein identifications"
)
parser.add_argument(
        "-i", "--input", required=True,
        help="Protein table including PeptideClassifier results"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="Output of unambiguous proteins"
)
parser.add_argument(
        "-a", "--ambiguous",
        help="Optional output of filtered out ambiguous proteins"
)
parser.add_argument(
        "-c", "--contams",
        help="Optional output of contaminant proteins"
)
parser.add_argument(
        "-p", "--contam-prefix", default="sp|",
        help="Prefix of contaminant proteins (default: 'sp|')"
)
args = parser.parse_args()

selected = []
discarded = []
contams = []

with open(args.input, "r") as fi:
    reader = csv.DictReader(fi, delimiter="\t")
    cols = reader.fieldnames

    for row in reader:
        if row["protein"].startswith(args.contam_prefix):
            contams.append(row)
        elif int(row["unique"]) > 0:
            selected.append(row)
        else:
            discarded.append(row)

print(f"selected:  {len(selected)}")
print(f"discarded: {len(discarded)}")
print(f"contams: {len(contams)}")

with open(args.output, "w") as fo:
    writer = csv.DictWriter(fo, delimiter="\t", fieldnames=cols)
    writer.writeheader()
    writer.writerows(selected)

if args.ambiguous:
    with open(args.ambiguous, "w") as fo:
        writer = csv.DictWriter(fo, delimiter="\t", fieldnames=cols)
        writer.writeheader()
        writer.writerows(discarded)

if args.contams:
    with open(args.contams, "w") as fo:
        writer = csv.DictWriter(fo, delimiter="\t", fieldnames=cols)
        writer.writeheader()
        writer.writerows(contams)
