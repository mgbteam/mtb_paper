#!/usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser(
    description="UpSetPlot of gene presence/absence"
)
parser.add_argument(
    "-p", "--panaroo", required=True,
    help="gene_presence_absence.csv from Panaroo"
)
parser.add_argument(
    "-r", "--roary", required=True,
    help="gene_presence_absence.csv from Roary"
)
parser.add_argument(
    "-C", "--common",
    help="Output common groups"
)
parser.add_argument(
    "-P", "--panaroo-unique",
    help="Output Panaroo unique groups"
)
parser.add_argument(
    "-R", "--roary-unique",
    help="Output Roary unique groups"
)
parser.add_argument(
    "-U", "--unique-matched",
    help="Output unique groups matched by genes"
)
args = parser.parse_args()


def read_presence_absence(path, startcol):
    groups = {}

    with open(path, "r") as f:
        reader = csv.reader(f)
        header = next(reader)

        if "Avg group size nuc" in header:
            startcol = header.index("Avg group size nuc") + 1
        else:
            startcol = header.index("Annotation") + 1

        for row in reader:
            group_genes = []

            for i in range(startcol, len(row)):
                row[i] = row[i].replace("\t", ";")

                for gene in row[i].split(";"):
                    if gene == "":
                        continue
                    if "refound" in gene:
                        continue

                    group_genes.append(gene.split(".")[0])

            groups[";".join(sorted(group_genes))] = row

    return header, groups


startcol = 14
panaroo_header, panaroo_groups = read_presence_absence(args.panaroo, startcol)
roary_header, roary_groups = read_presence_absence(args.roary, startcol)

common = {k: v for k, v in panaroo_groups.items() if k in roary_groups}
panaroo = {k: v for k, v in panaroo_groups.items() if k not in roary_groups}
roary = {k: v for k, v in roary_groups.items() if k not in panaroo_groups}

if args.common:
    with open(args.common, "w") as fo:
        writer = csv.writer(fo, delimiter="\t")
        writer.writerow(panaroo_header)
        writer.writerows(common.values())

if args.panaroo_unique:
    with open(args.panaroo_unique, "w") as fo:
        writer = csv.writer(fo, delimiter="\t")
        writer.writerow(panaroo_header)
        writer.writerows(panaroo.values())

if args.roary_unique:
    with open(args.roary_unique, "w") as fo:
        writer = csv.writer(fo, delimiter="\t")
        writer.writerow(roary_header)
        writer.writerows(roary.values())

if args.unique_matched:
    matched_panaroo_groups = set()
    matched_roary_groups = set()

    with open(args.unique_matched, "w") as fo:
        writer = csv.writer(fo, delimiter="\t")
        header2 = ["Gene", "Gene name", "Annotation"]
        header2.extend([f"Panaroo {n}" for n in panaroo_header[startcol:]])
        header2.extend([f"Roary {n}" for n in roary_header[startcol:]])
        writer.writerow(header2)

        for identifier, group in panaroo.items():
            row = group[:3] + group[startcol:]
            matching_groups = set()

            for gene in identifier.split(";"):
                for identifier2, group2 in roary.items():
                    if gene in identifier2.split(";"):
                        matched_panaroo_groups.add(identifier)
                        matched_roary_groups.add(identifier2)
                        fields = group2[startcol:]
                        group_id = ";".join(fields)

                        if group_id not in matching_groups:
                            matching_groups.add(group_id)
                            writer.writerow(row + fields)

        for identifier, group in panaroo.items():
            if identifier not in matched_panaroo_groups:
                writer.writerow(row + [""] * (len(group2) - startcol))

        for identifier2, group2 in roary.items():
            if identifier2 not in matched_roary_groups:
                row = group2[:3] + [""] * (len(group) - startcol)
                row += group2[startcol:]
                writer.writerow(row)

