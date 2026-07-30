#!/usr/bin/env python3

import argparse
import csv

from statistics import median

parser = argparse.ArgumentParser(
        description="Add information to pangenome table"
)
parser.add_argument(
        "-p", "--pangenome", required=True,
        help="gene_presence_absence.csv from Panaroo or Roary"
)
parser.add_argument(
        "-i", "--iptgxdbs", nargs="+", required=True,
        help="iPtgxDB fasta file per strain"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="Output table in tsv format",
)
parser.add_argument(
        "-r", "--refseq-prefix", required=True,
        help="Prefix of RefSeq gene identifiers"
)
args = parser.parse_args()

prot_lengths = {}

for fasta in args.iptgxdbs:
    with open(fasta, "r") as fi:
        for line in fi:
            if not line.startswith(">"):
                continue

            length = int(line.strip().split("_")[-1])
            identifier = line.strip().lstrip(">")
            anchor = identifier.split("|")[0]

            if identifier.startswith(args.refseq_prefix):
                if not any([s in anchor for s in ["_p_", "_+", "_-"]]):
                    identifier = anchor

            prot_lengths[identifier] = length

cluster_type_counts = {"RefSeq": 0, "Pseudo": 0, "New Start": 0, "Novel": 0}

outcols = [
        "Gene",
        "Annotation",
        "No. isolates",
        "No. sequences",
        "Group type",
        "Min. len. aa",
        "Max. len. aa",
        "Mdn. len. aa",
]

with open(args.pangenome, "r") as fi, open(args.output, "w") as fo:
    reader = csv.DictReader(fi, delimiter=",")
    cols = reader.fieldnames
    outcols.extend(cols[14:])
    writer = csv.DictWriter(fo, delimiter="\t", fieldnames=outcols, extrasaction="ignore")
    writer.writeheader()

    for row in reader:
        lengths = []
        types = []

        for col in cols[14:]:
            if row[col] == "":
                continue

            for full_id in row[col].split(";"):
                identifier = ".".join(full_id.split(".")[:-1])
                lengths.append(prot_lengths[identifier])

                if "novel" in row["Annotation"]:
                    types.append("Novel")
                elif "pseudogene" in row["Annotation"]:
                    types.append("Pseudo")
                elif "_ext" in full_id or "_red" in full_id:
                    types.append("New Start")
                else:
                    types.append("RefSeq")

        if "Novel" in types:
            prot_type = "Novel"
        elif "Pseudo" in types:
            prot_type = "Pseudo"
        elif "New Start" in types:
            prot_type = "New Start"
        else:
            prot_type = "RefSeq"

        row["Group type"] = prot_type
        row["Min. len. aa"] = min(lengths)
        row["Max. len. aa"] = max(lengths)
        row["Mdn. len. aa"] = median(lengths)
        writer.writerow(row)

        cluster_type_counts[prot_type] += 1

for cluster_type, count in cluster_type_counts.items():
    print(cluster_type, count)
