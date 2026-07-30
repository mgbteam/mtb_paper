#!/usr/bin/env python3

import argparse
import csv
import os

from ete3 import NCBITaxa

parser = argparse.ArgumentParser(
        description="Summarize conservation BLAST results"
)
parser.add_argument(
        "-b", "--blast-results", nargs="+", required=True,
        help="BLAST results as TSV file per protein"
)
parser.add_argument(
        "-t", "--tax-levels", nargs="+", required=True,
        help="Taxonomy levels to include in the summary"
)
parser.add_argument(
        "-c", "--min-coverage", default=60, type=int,
        help="Minimum coverage percentage (default: 60)"
)
parser.add_argument(
        "-i", "--min-identity", default=40, type=int,
        help="Minimum identity percentage (default: 40)"
)
parser.add_argument(
        "-e", "--max-evalue", default=0.01, type=float,
        help="Maximum e-value (default: 0.01)"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV file"
)
parser.add_argument(
        "-s", "--subtaxa-counts",
        help="Optionally output table of subtaxa counts"
)
parser.add_argument(
        "-d", "--taxdb",
        help="Optionally provide path to taxonomy4blast.sqlite3"
)
parser.add_argument(
        "-D", "--taxdump",
        help="Optionally provide path to taxonomy4blast.sqlite3"
)

args = parser.parse_args()

if args.taxdb and args.taxdump:
    ncbi = NCBITaxa(dbfile=args.taxdb, taxdump_file=args.taxdump)
elif args.taxdb:
    ncbi = NCBITaxa(dbfile=args.taxdb)
else:
    ncbi = NCBITaxa()

header = ["Protein"] + args.tax_levels
outrows = []

for infile in args.blast_results:
    with open(infile, "r") as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter="\t")
        protein = os.path.basename(infile.removesuffix(".tsv"))
        unique_lineages = set()
        tax_counts = {}

        for row in reader:
            if 100 * int(row["length"])/int(row["qlen"]) < args.min_coverage:
                continue

            if 100 * int(row["nident"])/int(row["length"]) < args.min_identity:
                continue

            if float(row["evalue"]) > args.max_evalue:
                continue

            for taxid in row["staxids"].split(";"):
                try:
                    lineage = ncbi.get_lineage(int(taxid))
                except ValueError:
                    continue

                names = ncbi.get_taxid_translator(lineage).values()
                lineage_path = ";".join(sorted(names))

                if lineage_path in unique_lineages:
                    continue

                unique_lineages.add(lineage_path)

                for name in names:
                    if name not in tax_counts:
                        tax_counts[name] = 0

                    tax_counts[name] += 1

        outrow = {"Protein": protein}
        prev_tax_count = 0

        for tax_level in args.tax_levels:
            if tax_level in tax_counts:
                outrow[tax_level] = tax_counts[tax_level] - prev_tax_count
                prev_tax_count = tax_counts[tax_level]
            else:
                outrow[tax_level] = 0
                prev_tax_count = 0

        outrows.append(outrow)

with open(args.output, "w") as outfile:
    writer = csv.DictWriter(outfile, delimiter="\t", fieldnames=header)
    writer.writeheader()
    writer.writerows(outrows)

if args.subtaxa_counts:
    counts = {t: len(ncbi.get_descendant_taxa(t)) for t in args.tax_levels}

    with open(args.subtaxa_counts, "w") as fo:
        writer = csv.writer(fo, delimiter="\t")
        writer.writerow(["Level", "Subtaxa"])

        for taxon, count in counts.items():
            writer.writerow([taxon, count])
