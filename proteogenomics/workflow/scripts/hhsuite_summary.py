#!/usr/bin/env python3

import argparse
import csv
import os

parser = argparse.ArgumentParser(
        description="Summarize hhsuite results"
)
parser.add_argument(
        "-i", "--input", nargs="+", required=True,
        help="hhsuite results in txt format, one file per protein"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="output in TSV format"
)
parser.add_argument(
        "-e", "--exclude", nargs="+",
        help="list of keywords to exclude matches"
)
args = parser.parse_args()

outrows = []

for file in args.input:
    with open(file, "r") as fi:
        protein = os.path.splitext(os.path.basename(file))[0]

        for line in fi:
            if line.startswith(">"):
                description = line.split(" ", 1)[1].split(" n=")[0]

                if args.exclude:
                    if any([kw.lower() in description.lower() for kw in args.exclude]):
                        continue

                probability = next(fi).removeprefix("Probab=").split(" ", 1)[0]
                outrows.append([protein, probability, description])
                break

with open(args.output, "w") as fo:
    writer = csv.writer(fo, delimiter="\t")
    writer.writerow(["protein", "hhsuite probability", "hhsuite description"])
    writer.writerows(outrows)
