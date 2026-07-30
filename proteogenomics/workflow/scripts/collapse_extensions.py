#!/usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser(
        description="Collapse iPtgxDB clusters of novel proteins and RefSeq extensions"
)
parser.add_argument(
        "-i", "--input", metavar="FILE", required=True,
        help="Input tsv file of proteins"
)
parser.add_argument(
        "-o", "--output", metavar="FILE", required=True,
        help="Output protein tsv file with collapsed clusters"
)
args = parser.parse_args()

cols_to_sum = [
    "total - psms",
    "total - peptides",
    "total - spectra",
    "1a", "1b", "2a", "2b", "3a", "3b",
]

clusters = {}

with open(args.input, "r") as fi:
    reader = csv.DictReader(fi, delimiter="\t")
    cols = reader.fieldnames

    for row in reader:
        # convert number columns from string to int or float
        for col in cols_to_sum:
            if col in row:
                row[col] = int(row[col])

        row["count proteins"] = int(row["count proteins"])
        row["sequence length"] = int(row["sequence length"])
        row["total - abundance"] = float(row["total - abundance"])

        # group novels and RefSeq extensions (but not anchors) by cluster
        if row["protein-type"] == "RefSeq protein":
            clusters[row["protein"] + "_anchor"] = [row]
        else:
            if row["protein-reference-id"] not in clusters:
                clusters[row["protein-reference-id"]] = []
            
            clusters[row["protein-reference-id"]].append(row)

outrows = []
merged_seqs = []

for refid, entries in clusters.items():
    if len(entries) == 1:
        # if only one entry is present for this cluster just add it
        outrows.append(entries[0])
    else:
        # if multiple entries are present merge values of shorter cluster entries to longest one
        entries_by_size = sorted(entries, key=lambda row: int(row["protein"].split("_")[-1]))
        longest_entry = entries_by_size[-1]
        base = longest_entry["protein"].split("|")[0]

        # update columns of longest entry
        for i, row in enumerate(entries_by_size[:-1]):
            base = longest_entry["protein"].split("|")[0]

            # if this is a reduction and the anchor is not a pseudogene, add it as separate entry
            if "_-" in base and "pseudogene" not in row["protein-type"]:
                outrows.append(row)
                continue

            # most counts can be summed (peptides, psms, etc.)
            for col in cols_to_sum:
                if col in longest_entry:
                    longest_entry[col] += row[col]

            # if this extension has a more specific peptide, update the protein group
            if row["count proteins"] < longest_entry["count proteins"]:
                longest_entry["protein-group"] = row["protein-group"]
                longest_entry["count proteins"] = row["count proteins"]

            # use the larger abundance of the two entries
            if row["total - abundance"] > longest_entry["total - abundance"]:
                longest_entry["total - abundance"] = row["total - abundance"]

        outrows.append(longest_entry)

# update lengths to be the total length of the protein, not that of the iPtgxDB entry
for outrow in outrows:
    outrow["sequence length"] = int(outrow["protein"].split("_")[-1])

with open(args.output, "w") as fo:
    writer = csv.DictWriter(fo, delimiter="\t", fieldnames=cols)
    writer.writeheader()
    writer.writerows(outrows)