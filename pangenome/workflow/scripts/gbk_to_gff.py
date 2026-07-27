#!/usr/bin/env python3

import argparse
import csv

from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Convert GenBank file to Roary compatible GFF"
)
parser.add_argument(
    "-i", "--input", required=True, help="Input GenBank file"
)
parser.add_argument(
    "-o", "--output", required=True, help="Output GFF file"
)
parser.add_argument(
    "-n", "--no-pseudo", action="store_true", help="Exclude pseudogenes"
)
args = parser.parse_args()

pseudo_keys = ["pseudo", "pseudogene"]

with open(args.output, "w") as fo:
    for record in SeqIO.parse(args.input, "genbank"):
        writer = csv.writer(fo, delimiter="\t")
        prev_locus_tag = None
        subfeature_count = 0

        for feature in record.features:
            if "locus_tag" not in feature.qualifiers:
                continue

            locus_tag = feature.qualifiers["locus_tag"][0]
            pid = locus_tag
            is_pseudo = False
            attrs = []

            if locus_tag != prev_locus_tag:
                subfeature_count = 0
                prev_locus_tag = locus_tag
            else:
                subfeature_count += 1
                attrs.append(f"Parent={locus_tag}")

            for key, val in feature.qualifiers.items():
                if key == "translation":
                    continue

                if key in pseudo_keys:
                    key = "pseudo"
                    val = ["true"]
                    is_pseudo = True

                attrs.append(f"{key}={val[0]}")

            if args.no_pseudo and is_pseudo:
                continue

            for i, part in enumerate(feature.location.parts):
                if i > 0:
                    subfeature_count += 1

                if subfeature_count == 0:
                    id_attrs = [f"ID={pid}"] + attrs
                else:
                    id_attrs = [f"ID={pid}.{subfeature_count}"] + attrs

                if feature.type == "CDS":
                    if i == 0:
                        phase = 0
                    elif part.strand == 1:
                        phase = (3-((part.start-feature.location.start+1) % 3)) % 3
                    else:
                        phase = (3-((part.end-feature.location.end+1) % 3)) % 3
                else:
                    phase = "."

                writer.writerow([
                        record.id,
                        "GenBank",
                        feature.type,
                        part.start + 1,
                        part.end,
                        ".",
                        "-" if part.strand < 0 else "+",
                        phase,
                        ";".join(id_attrs)
                ])

        writer.writerow(["##FASTA"])
        SeqIO.write(record, fo, "fasta")
