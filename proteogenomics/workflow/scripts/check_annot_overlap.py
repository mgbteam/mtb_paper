#!/usr/bin/env python3

import argparse
import csv

from Bio import SeqIO

parser = argparse.ArgumentParser(
        description="Check overlap between novel proteins and annotation"
)
parser.add_argument(
        "-i", "--input", required=True,
        help="Input TSV file containing novel proteins"
)
parser.add_argument(
        "-a", "--annot", required=True,
        help="Input annotation in GenBank format"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="output TSV file containing overlap info"
)
parser.add_argument(
        "-n", "--name", default="annotation",
        help="Name of annotation (default: annotation)"
)
args = parser.parse_args()

refseq_by_stop = {}

for record in SeqIO.parse(args.annot, "genbank"):
    refseq_by_stop[record.id] = {"+": {}, "-": {}}

    for feature in record.features:
        if feature.type != "CDS":
            continue

        strand = "-" if feature.location.strand < 0 else "+"

        if strand == "+":
            stop = int(feature.location.end)
        else:
            stop = int(feature.location.start) + 1

        locus_tag = feature.qualifiers["locus_tag"][0]
        refseq_by_stop[record.id][strand][stop] = locus_tag

with open(args.input, "r") as fi, open(args.output, "w") as fo:
    reader = csv.DictReader(fi, delimiter="\t")
    writer = csv.writer(fo, delimiter="\t")
    writer.writerow(["protein", args.name])

    for row in reader:
        pos = row["protein"].split("|")[-1].split("_")
        chrom = "_".join(pos[:-5])
        start = int(pos[-5])
        stop = int(pos[-4])
        strand = pos[-3][0]

        tag = refseq_by_stop.get(chrom, {}).get(strand, {}).get(stop, None)
        tag = "" if tag is None else tag
        writer.writerow([row["protein"], tag])
