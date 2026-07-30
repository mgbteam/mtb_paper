#!/usr/bin/env python3

import argparse
import csv

from Bio import SeqIO

parser = argparse.ArgumentParser(
        description="Reformat Abacus results"
)
parser.add_argument(
        "-i", "--input", required=True,
        help="Protein table from Abacus"
)
parser.add_argument(
        "-f", "--fasta", required=True,
        help="Database fasta file"
)
parser.add_argument(
        "-m", "--mapping", required=True,
        help="Database mapping in TSV format"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="Reformatted output"
)
parser.add_argument(
        "-r", "--refseq-prefix", default="REF",
        help="RefSeq identifier prefix (default: REF)"
)
parser.add_argument(
        "-c", "--contam-prefix", default="sp|",
        help="Contaminant prefix (default: sp|)"
)
args = parser.parse_args()

database = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

with open(args.mapping, "r") as fi:
    reader = csv.reader(fi, delimiter="\t")
    header = next(reader)
    id_cols = {}

    for col_index, col in enumerate(header):
        if col.endswith(" ID"):
            id_cols[col.removesuffix(" ID")] = col_index

    mapped_ids_index = header.index("mapped IDs")
    id_info = {}

    for row in reader:
        mapped_ids = row[mapped_ids_index:]
        major_id_found = False
        offset = 1
        major_id = mapped_ids[0]

        if major_id.split("|")[0].endswith("_p"):
            prot_type = "pseudogene"
        else:
            prot_type = "protein"

        if major_id in database:
            major_id_found = True
        else:
            while offset < len(mapped_ids):
                major_id = mapped_ids[offset]
                offset += 1

                if major_id in database:
                    major_id_found = True
                    break

        if not major_id_found:
            continue

        additional_ids = []

        for i in range(offset, len(mapped_ids)):
            additional_ids.append(mapped_ids[i])

        for source, col_index in id_cols.items():
            if row[col_index] != "":
                id_info[major_id] = {
                    "type": f"{source} {prot_type}",
                    "reference-id": major_id,
                }

                for additional_id in additional_ids:
                    id_info[additional_id] = {
                        "type": f"extension/reduction to {source} {prot_type}",
                        "reference-id": major_id,
                    }

                break

with open(args.input, "r") as fi, open(args.output, "w") as fo:
    reader = csv.DictReader(fi, delimiter="\t")
    intensity_cols = [n for n in reader.fieldnames if n.endswith("Intensity")]

    outcols = [
        "protein",
        "protein-group",
        "count proteins",
        "sequence length",
        "protein-type",
        "protein-reference-id",
        "total - psms",
        "total - peptides",
        "total - spectra",
        "total - abundance",
    ]

    writer = csv.DictWriter(fo, delimiter="\t", fieldnames=outcols)
    writer.writeheader()

    for row in reader:
        if row["Indistinguishable Proteins"] == "":
            protein_group = row["Protein"]
            n_proteins = 1
        else:
            indistinguishable = row["Indistinguishable Proteins"].split(", ")
            protein_group = ";".join([row["Protein"]] + indistinguishable)
            n_proteins = 1 + len(indistinguishable)

        intensities = [float(row[col]) for col in intensity_cols]
        seq_len = len(database[row["Protein"]].seq) if row["Protein"] in database else -1

        if row["Protein"].startswith(args.contam_prefix):
            prot_info = {"type": "Contaminant protein", "reference-id": row["Protein"]}
        else:
            prot_info = id_info.get(row["Protein"], {"type": "", "reference-id": ""})

        writer.writerow({
            "protein": row["Protein"],
            "protein-group": protein_group,
            "count proteins": n_proteins,
            "sequence length": seq_len,
            "protein-type": prot_info["type"],
            "protein-reference-id": prot_info["reference-id"],
            "total - psms": row["Combined Total Spectral Count"],
            "total - peptides": row["Combined Total Peptides"],
            "total - spectra": row["Combined Total Spectral Count"],
            "total - abundance": sum(intensities) / len(intensity_cols),
        })