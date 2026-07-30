#!/usr/bin/env python3

import argparse
import csv

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

parser = argparse.ArgumentParser(
        description="Calculate physicochemical properties of protein"
)
parser.add_argument(
        "-i", "--input", required=True,
        help="FASTA file with protein sequences"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="Output file in TSV format"
)
args = parser.parse_args()

protein_seqs = list(SeqIO.parse(args.input, "fasta"))
results = []

for seq in protein_seqs:
    analysis = ProteinAnalysis(str(seq.seq))
    row = [
        seq.id,
        analysis.molecular_weight(),
        analysis.instability_index(),
        analysis.aromaticity(),
        analysis.isoelectric_point()
    ]

    results.append(row)

with open(args.output, "w") as outfile:
    writer = csv.writer(outfile, delimiter="\t")
    writer.writerow([
        "Protein",
        "MW",
        "Instability",
        "Aromaticity",
        "Isoelectric point"
    ])
    writer.writerows(results)
