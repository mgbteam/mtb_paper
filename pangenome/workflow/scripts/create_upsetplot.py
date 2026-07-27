#!/usr/bin/env python3

import argparse
import csv
import re

import matplotlib.pyplot as plt

import upsetplot

# Parse command line arguments
parser = argparse.ArgumentParser(
    description="UpSetPlot of gene presence/absence"
)
parser.add_argument(
    "-i", "--input", required=True,
    help="gene_presence_absence.csv from Panaroo"
)
parser.add_argument(
    "-o", "--output", required=True,
    help="Output UpSetPlot to this file"
)
parser.add_argument(
    "-r", "--no-refound", action="store_true",
    help="Exclude refound genes"
)
parser.add_argument(
    "-p", "--no-pseudo", action="store_true",
    help="Exclude pseudogenes"
)
parser.add_argument(
    "-R", "--refound-only", action="store_true",
    help="Only plot refound genes"
)
parser.add_argument(
    "-P", "--pseudo-only", action="store_true",
    help="Only plot pseudogenes"
)
parser.add_argument(
    "-s", "--subset-regex", type=str,
    help="Subset to gene products which match regex"
)
parser.add_argument(
    "-l", "--subset-list", type=str,
    help="List of locus tags to subset genes by"
)
parser.add_argument(
    "-c", "--color", nargs="+", action="append",
    help="Color and strains to color (space separated)"
)
args = parser.parse_args()

if args.subset_regex:
    regex = re.compile(args.subset_regex)

if args.subset_list:
    with open(args.subset_list, "r") as fi:
        subset_locus_tags = [line.strip() for line in fi]

# Read in the gene presence/absence table
with open(args.input, "r") as f:
    reader = csv.DictReader(f)

    if "Avg group size nuc" in reader.fieldnames:
        start = reader.fieldnames.index("Avg group size nuc") + 1
    else:
        start = reader.fieldnames.index("Annotation") + 1

    strains = reader.fieldnames[start:][::-1]
    gene_dict = {strain: [] for strain in strains}

    for row in reader:
        for strain in strains:
            if row[strain] == "":
                continue

            if all(["refound" in pid for pid in row[strain].split(";")]):
                if args.no_refound:
                    continue
            else:
                if args.refound_only:
                    continue

            if all(["pseudo" in pid for pid in row[strain].split(";")]):
                if args.no_pseudo:
                    continue
            else:
                if args.pseudo_only:
                    continue

            if args.subset_regex:
                if not regex.search(row["Annotation"]):
                    continue

            if args.subset_list:
                has_subset_gene = False

                for strain2 in strains:
                    for locus_tag in row[strain2].split(";"):
                        for subset_locus_tag in subset_locus_tags:
                            if locus_tag.startswith(subset_locus_tag):
                                has_subset_gene = True
                                break
                        if has_subset_gene:
                            break
                    if has_subset_gene:
                        break

                if not has_subset_gene:
                    continue

            gene_dict[strain].append(row["Gene"])

# Create an upset plot
try:
    upset = upsetplot.UpSet(
            upsetplot.from_contents(gene_dict),
            sort_by="cardinality",
            sort_categories_by=None,
            show_counts=True
    )

    if args.color:
        for colorstrains in args.color:
            color = colorstrains[0]
            strains = colorstrains[1:]
            upset.style_categories(strains, bar_facecolor=color)

    upset.plot()
except IndexError:
    plt.close()
    plt.figure()
    plt.plot([], [])

plt.savefig(args.output)
plt.close()
