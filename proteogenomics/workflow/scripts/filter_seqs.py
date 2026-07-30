#!/usr/bin/env python3

import argparse

from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Filter sequences in fasta file based on length"
)
parser.add_argument(
    "-s", "--shortest", type=int,
    help="Sequences shorter than this will be filtered out"
)
parser.add_argument(
    "-l", "--longest", type=int,
    help="Sequences longer than this will be filtered out"
)
parser.add_argument(
    "-i", "--input", required=True,
    help="Input file in fasta format"
)
parser.add_argument(
    "-o", "--output", required=True,
    help="Output file in fasta format"
)
args = parser.parse_args()

records = []

for record in SeqIO.parse(args.input, "fasta"):
    if args.shortest and len(record.seq) < args.shortest:
        continue
    if args.longest and len(record.seq) > args.longest:
        continue

    records.append(record)

SeqIO.write(records, args.output, "fasta")