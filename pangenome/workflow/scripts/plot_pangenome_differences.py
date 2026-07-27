#!/usr/bin/env python3

import argparse
import csv

import matplotlib.pyplot as plt

from matplotlib_venn import venn2

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
    "-v", "--venn", required=True,
    help="Output venn diagramm"
)
parser.add_argument(
    "-b", "--boxplot", required=True,
    help="Output venn diagramm"
)
parser.add_argument(
    "-i", "--include-refound", action="store_true",
    help="Include refound genes from Panaroo"
)
args = parser.parse_args()


def read_presence_absence(path):
    groups = set()

    with open(path, "r") as f:
        reader = csv.DictReader(f)

        if "Avg group size nuc" in reader.fieldnames:
            start = reader.fieldnames.index("Avg group size nuc") + 1
        else:
            start = reader.fieldnames.index("Annotation") + 1

        strains = reader.fieldnames[start:]

        for row in reader:
            group_genes = []

            for strain in strains:
                for gene in row[strain].replace("\t", ";").split(";"):
                    if gene == "":
                        continue
                    if "refound" in gene and not args.include_refound:
                        continue

                    group_genes.append(gene.split(".")[0])

            groups.add(";".join(sorted(group_genes)))

    return groups


panaroo_groups = read_presence_absence(args.panaroo)
roary_groups = read_presence_absence(args.roary)

common = [len(g.split(";")) for g in panaroo_groups if g in roary_groups]
panaroo = [len(g.split(";")) for g in panaroo_groups if g not in roary_groups]
roary = [len(g.split(";")) for g in roary_groups if g not in panaroo_groups]

print("\nPanaroo", len(panaroo_groups), sep="\t")
print("Roary", len(roary_groups), sep="\t")
print("Common", len(common), sep="\t")

venn2((panaroo_groups, roary_groups), set_labels=("Panaroo", "Roary"))
plt.savefig(args.venn)
plt.close()

vals = [panaroo, common, roary]
labels = ["Panaroo Only", "Common", "Roary Only"]

plt.figure()
plt.boxplot(vals, labels=labels)
plt.ylabel("Genes in Group")
plt.savefig(args.boxplot)
