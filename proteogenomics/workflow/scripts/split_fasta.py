#!/usr/bin/env python3

import argparse
import os

from Bio import SeqIO

parser = argparse.ArgumentParser(
        description="Split fasta file containing multiple sequences"
)
parser.add_argument(
        "-i", "--input", required=True,
        help="Input fasta file containing multiple sequences"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="Output folder containing individual fasta files"
)
args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)

for record in SeqIO.parse(args.input, "fasta"):
    SeqIO.write([record], f"{args.output}/{record.id}.faa", "fasta")
