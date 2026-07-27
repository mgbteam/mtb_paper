#!/usr/bin/env python3

import argparse

from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Convert GenBank file to Roary compatible GFF"
)
parser.add_argument(
    "-i", "--input", required=True, help="Input GenBank file"
)
parser.add_argument(
    "-o", "--output", required=True, help="Output Fasta file"
)
args = parser.parse_args()

records = list(SeqIO.parse(args.input, "genbank"))
SeqIO.write(records, args.output, "fasta")
