#!/usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser(
    description="Remove refound genes from presence/absence table"
)
parser.add_argument(
    "-i", "--input", required=True,
    help="presence absence/table from Panaroo"
)
parser.add_argument(
    "-o", "--output", required=True,
    help="presence/absence table with refound genes removed"
)
args = parser.parse_args()

with open(args.input, "r") as fi, open(args.output, "w") as fo:
    reader = csv.DictReader(fi, delimiter=",")
    writer = csv.DictWriter(fo, delimiter=",", fieldnames=reader.fieldnames)
    writer.writeheader()

    for row in reader:
        for key in row:
            if "_refound_" in row[key]:
                row[key] = ""

        writer.writerow(row)
