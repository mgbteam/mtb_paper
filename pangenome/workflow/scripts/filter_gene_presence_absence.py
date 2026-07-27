#!/usr/bin/env python

import argparse
import csv

from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Filter presence/absence matrix by pan-genome fasta"
)
parser.add_argument(
    "-r", "--reference", required=True,
    help="pan_genome_reference.fa from Roary or Panaroo"
)
parser.add_argument(
    "-m", "--matrix", required=True,
    help="gene_presence_absence.csv from Roary or Panaroo"
)
parser.add_argument(
    "-o", "--output", required=True,
    help="Filtered gene_presence_absence.csv as output"
)
args = parser.parse_args()

groups = set()

for record in SeqIO.parse(args.reference, "fasta"):
    parts = record.description.split(" ")

    if len(parts) > 1:
        groups.add(" ".join(parts[1:]))
    else:
        groups.add(parts[0])

with open(args.matrix, "r") as fi, open(args.output, "w") as fo:
    reader = csv.DictReader(fi, delimiter=",")
    writer = csv.DictWriter(fo, delimiter=",", fieldnames=reader.fieldnames)
    writer.writeheader()

    removed_groups = 0

    for row in reader:
        if row["Gene"] in groups:
            writer.writerow(row)
        else:
            removed_groups += 1

print(f"Removed {removed_groups} groups")
