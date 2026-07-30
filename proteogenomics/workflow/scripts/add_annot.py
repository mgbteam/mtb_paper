#!/usr/bin/env python3

import argparse
import csv

from Bio import SeqIO

parser = argparse.ArgumentParser(
        description="Add RefSeq annotation to protein table"
)
parser.add_argument(
        "-i", "--input", required=True, help="Input protein table"
)
parser.add_argument(
        "-a", "--annot", required=True, help="Input GenBank annotation"
)
parser.add_argument(
        "-o", "--output", required=True, help="Output protein table"
)
args = parser.parse_args()

annotation = {}

for record in SeqIO.parse(args.annot, "genbank"):
    for feature in record.features:
        if "product" in feature.qualifiers:
            annotation[feature.qualifiers["locus_tag"][0]] = feature.qualifiers["product"][0]

with open(args.input, "r") as fi, open(args.output, "w") as fo:
    reader = csv.DictReader(fi, delimiter="\t")
    cols = reader.fieldnames + ["Annotation"]
    writer = csv.DictWriter(fo, delimiter="\t", fieldnames=cols)
    writer.writeheader()

    for row in reader:
        locus_tag = row["protein"].split("|")[0].split("_+")[0].split("_-")[0].removesuffix("_fCDS")

        if locus_tag in annotation:
            row["Annotation"] = annotation[locus_tag]

        writer.writerow(row)