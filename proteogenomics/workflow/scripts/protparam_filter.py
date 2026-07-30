#!/usr/bin/env python3

import argparse
import csv


def parse_args():
    parser = argparse.ArgumentParser(
           description="Filter proteins based on PSMs, peptides and length"
    )
    parser.add_argument(
            "-p", "--proteins", metavar="FILE", required=True,
            help="Input tsv file of proteins"
    )
    parser.add_argument(
            "-m", "--min-psms", metavar=("NAME", "INT"),
            nargs=2, action="append",
            help="Annotation name and minimum number of PSMs"
    )
    parser.add_argument(
            "-M", "--min-peps", metavar="INT", type=int, default=2,
            help="Minimum number of peptides for (default: 2)"
    )
    parser.add_argument(
            "-w", "--min-weight", metavar="INT", type=int, default=15000,
            help="Minimum weight for peptide filter in Dalton (default: 15000)"
    )
    parser.add_argument(
            "-o", "--output", metavar="FILE", required=True,
            help="Output tsv file containing filtered proteins"
    )
    parser.add_argument(
            "-O", "--output-discarded", metavar="FILE",
            help="Output tsv file containing discarded proteins"
    )

    return parser.parse_args()


def read_file(file, id_col, pep_col, psm_col, min_weight, min_peps, min_psms):
    proteins = {"selected": [], "discarded": []}

    with open(file, "r") as fi:
        reader = csv.DictReader(fi, delimiter="\t")

        for row in reader:
            identified = False
            discarded = False

            # discard the protein if it does not pass the peptide filter
            above_min_weight = float(row["MW"]) > min_weight
            too_few_peps = int(row[pep_col]) < min_peps

            if above_min_weight and too_few_peps:
                discarded = True

            # discard the protein if it does not pass the psm filter
            for annot, min_psms_annot in min_psms:
                if annot in row[id_col]:
                    identified = True

                    if int(row[psm_col]) < int(min_psms_annot):
                        discarded = True
                        break

            if identified:
                if discarded:
                    proteins["discarded"].append(row)
                else:
                    proteins["selected"].append(row)
            else:
                print("Warning, no psm threshold for " + row["protein"])

    return proteins


def main():
    args = parse_args()

    proteins = read_file(
            args.proteins,
            "protein-type",
            "total - peptides",
            "total - psms",
            args.min_weight,
            args.min_peps,
            args.min_psms
    )

    for category, entries in proteins.items():
        print(category, len(entries), sep="\t")

    with open(args.output, "w") as fo:
        cols = proteins["selected"][0].keys()
        writer = csv.DictWriter(fo, delimiter="\t", fieldnames=cols)
        writer.writeheader()
        writer.writerows(proteins["selected"])

    if args.output_discarded:
        with open(args.output_discarded, "w") as fo:
            cols = proteins["discarded"][0].keys()
            writer = csv.DictWriter(fo, delimiter="\t", fieldnames=cols)
            writer.writeheader()
            writer.writerows(proteins["discarded"])


if __name__ == "__main__":
    main()
