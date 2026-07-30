#!/usr/bin/env python3

import argparse
import csv

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(
        description="Extract full sequences of proteins based on iPtgxDB ID"
)
parser.add_argument(
        "-i", "--input", metavar="FILE", required=True,
        help="Input tsv file of proteins"
)
parser.add_argument(
        "-g", "--genome", metavar="FILE", required=True,
        help="Input fasta file of genome sequence"
)
parser.add_argument(
        "-o", "--output", metavar="FILE", required=True,
        help="Output fasta file with full sequences of proteins"
)
args = parser.parse_args()

chroms = SeqIO.index(args.genome, "fasta")
proteins = []

with open(args.input, "r") as infile, open(args.output, "w") as outfile:
    reader = csv.DictReader(infile, delimiter="\t")

    for row in reader:
        # extract chromosome, start and stop from line
        pos = row["protein"].split("|")[-1].split("_")
        chrom = "_".join(pos[:-5])
        strand = pos[-3][0]
        start = int(pos[-5])
        stop = int(pos[-4])

        # get DNA sequence, take reverse complement if start bigger than stop
        if strand == "-":
            dna_seq = chroms[chrom][stop-1:start].reverse_complement().seq
        else:
            dna_seq = chroms[chrom][start-1:stop].seq

        # translate and make the first amino acid a methionine
        aa_seq = MutableSeq(dna_seq.translate()).replace("I", "L")
        aa_seq[0] = "M"
        aa_seq = aa_seq.rstrip("*")

        seqrecord = SeqRecord(
            aa_seq,
            id=row["protein"],
            description=row["protein"],
        )

        # add translated sequence to output
        proteins.append(seqrecord)

SeqIO.write(proteins, args.output, "fasta")
