#!/usr/bin/env python3

import argparse
import csv

from Bio import SeqIO

parser = argparse.ArgumentParser(
        description="Reformat filter and report results"
)
parser.add_argument(
        "-i", "--input", required=True,
        help="Protein table from report command"
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
parser.add_argument(
        "-d", "--decoy-prefix", default="rev_",
        help="Decoy prefix (default: rev_)"
)
parser.add_argument(
        "-e", "--entrapment_suffix", default="_p_target",
        help="Entrapment suffix (default: _p_target)"
)
parser.add_argument(
        "-C", "--filter-contams", action="store_true",
        help="Filter out contaminants"
)
parser.add_argument(
        "-D", "--filter-decoys", action="store_true",
        help="Filter out decoys"
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
                    "type": source,
                    "reference-id": major_id,
                }

                for additional_id in additional_ids:
                    id_info[additional_id] = {
                        "type": f"extension/reduction to {source}",
                        "reference-id": major_id,
                    }

                break

with open(args.input, "r") as fi:
    reader = csv.DictReader(fi, delimiter="\t")
    intensity_cols = [n for n in reader.fieldnames if n.endswith("Intensity")]
    outrows = []

    for row in reader:
        if row["Indistinguishable Proteins"] == "":
            protein_group = row["Protein"]
            n_proteins = 1
        else:
            indistinguishable = row["Indistinguishable Proteins"].split(", ")
            protein_group = ";".join([row["Protein"]] + indistinguishable)
            n_proteins = 1 + len(indistinguishable)

        intensities = [float(row[col]) for col in intensity_cols]

        prot_id = row["Protein"].removeprefix("rev_").removesuffix("_p_target")
        seq_len = len(database[prot_id].seq) if prot_id in database else -1
        prot_info = id_info.get(prot_id, {"type": "", "reference-id": ""})

        if row["Protein"].removeprefix("rev_").startswith(args.contam_prefix):
            prot_info["type"] = "Contaminant"

        anchor = row["Protein"].split("|")[0]

        if "_p_-" in anchor or "_p_+" in anchor or anchor.endswith("_p"):
            suffix = " pseudogene"
        else:
            suffix = " protein"

        if row["Protein"].endswith(args.entrapment_suffix):
            suffix += " entrapment"
        if row["Protein"].startswith(args.decoy_prefix):
            suffix += " decoy"

        outrows.append({
            "protein": row["Protein"],
            "protein-group": protein_group,
            "count proteins": n_proteins,
            "sequence length": seq_len,
            "protein-type": prot_info["type"] + suffix,
            "protein-reference-id": prot_info["reference-id"],
            "total - psms": row["Total Spectral Count"],
            "total - peptides": row["Total Peptides"],
            "total - spectra": row["Total Spectral Count"],
            "total - abundance": sum(intensities) / len(intensity_cols),
            "q_value": float(row["Protein Qvalue"]),
            "probability": float(row["Protein Probability"]),
        })


if args.filter_contams or args.filter_decoys:
    outrows_filtered = []

    for row in outrows:
        if args.filter_contams and "Contaminant" in row["protein-type"]:
            continue

        if args.filter_decoys and "decoy" in row["protein-type"]:
            continue

        outrows_filtered.append(row)
else:
    outrows_filtered = outrows

outrows_sorted = sorted(outrows_filtered, key=lambda d: (d["q_value"], -d["probability"]))

for i in range(len(outrows_sorted)):
    outrows_sorted[i]["score"] = i+1

with open(args.output, "w") as fo:
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
        "q_value",
        "probability",
        "score",
    ]

    writer = csv.DictWriter(fo, delimiter="\t", fieldnames=outcols)
    writer.writeheader()
    writer.writerows(outrows_sorted)
