#!/usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser(
        description="Split identified proteins into categories"
)
parser.add_argument(
        "-i", "--input", required=True,
        help="Input protein table in TSV format"
)
parser.add_argument(
        "-a", "--all",
        help="Output RefSeq protein"
)
parser.add_argument(
        "-r", "--refseq",
        help="Output RefSeq protein"
)
parser.add_argument(
        "-p", "--pseudo",
        help="Output expressed pseudogenes"
)
parser.add_argument(
        "-s", "--starts",
        help="Output novel start sites"
)
parser.add_argument(
        "-n", "--novels",
        help="Output novel proteins"
)
parser.add_argument(
        "-N", "--novelty",
        help="Output novel protein, start sites and pseudogenes combined"
)
args = parser.parse_args()

outrows = {
    "all": {"file": args.all, "rows": []},
    "refseq": {"file": args.refseq, "rows": []},
    "pseudo": {"file": args.pseudo, "rows": []},
    "starts": {"file": args.starts, "rows": []},
    "novels": {"file": args.novels, "rows": []},
    "novelty": {"file": args.novelty, "rows": []},
}

with open(args.input, "r") as fi:
    reader = csv.DictReader(fi, delimiter="\t")
    cols = reader.fieldnames

    for row in reader:
        outrows["all"]["rows"].append(row)

        if row["protein-type"] == "RefSeq protein":
            outrows["refseq"]["rows"].append(row)
        elif row["protein-type"] == "RefSeq pseudogene":
            outrows["pseudo"]["rows"].append(row)
            outrows["novelty"]["rows"].append(row)
        elif row["protein-type"] == "extension/reduction to RefSeq protein":
            outrows["starts"]["rows"].append(row)
            outrows["novelty"]["rows"].append(row)
        else:
            outrows["novels"]["rows"].append(row)
            outrows["novelty"]["rows"].append(row)

print("groups", "proteins", sep="\t")

for name, category in outrows.items():
    print(name, len(category["rows"]), sep="\t")

    if category["file"]:
        with open(category["file"], "w") as fo:
            writer = csv.DictWriter(fo, delimiter="\t", fieldnames=cols)
            writer.writeheader()
            writer.writerows(category["rows"])
