#!/usr/bin/env python3

import argparse
import csv

from Bio import SeqIO

parser = argparse.ArgumentParser(
        description="Create a GFF file containing annotated and novel proteins"
)
parser.add_argument(
        "-a", "--annot", required=True,
        help="Input GenBank file of annotation"
)
parser.add_argument(
        "-p", "--proteins", required=True,
        help="Input table of detected proteins"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="Output GFF file"
)
parser.add_argument(
        "-n", "--no-gene-attr", action="store_true",
        help="Do not add gene name to attribute columns"
)
parser.add_argument(
        "-t", "--tag-ptx", action="store_true",
        help="Add tag to RefSeq identifier if detected in proteomics"
)
parser.add_argument(
        "-f", "--fasta-append", action="store_true",
        help="Optionally append fasta sequence to GFF file"
)
args = parser.parse_args()

pseudo_keys = ["pseudo", "pseudogene", "ribosomal_slippage"]

identified_refseq_prots = set()
identified_refseq_exts = {}
identified_refseq_reds = {}

gff_entries = {}
seqs = {}

# read in the list of detected proteins
with open(args.proteins, "r") as fi:
    reader = csv.DictReader(fi, delimiter="\t")

    for row in reader:
        if row["protein-type"].endswith("RefSeq protein"):
            if row["protein-type"] == "RefSeq protein":
                protid = row["protein"].split("|")[0]
                identified_refseq_prots.add(protid)
            if row["protein-type"] == "extension/reduction to RefSeq protein":
                protid = row["protein"].split("|")[0]

                if "_+" in protid:
                    protid_split = protid.split("_+")
                    protid = protid_split[0]
                    length = int(protid_split[1].split("aa")[0])
                    identified_refseq_exts[protid] = length
                elif "_-" in protid:
                    protid_split = protid.split("_-")
                    protid = protid_split[0]
                    length = int(protid_split[1].split("aa")[0])
                    identified_refseq_reds[protid] = length
                else:
                    print(f"Unable to parse extensions/reduction {protid}")

            continue
        elif "pseudogene" in row["protein-type"]:
            prottype = "pseudogene"
        else:
            prottype = "novel"

        protid = row["protein"]
        pos = protid.split("|")[-1].split("_")
        chrom = "_".join(pos[:-5])
        start = int(pos[-5])
        end = int(pos[-4])
        strand = pos[-3][0]

        if chrom not in gff_entries:
            gff_entries[chrom] = []

        if strand == "-":
            start, end = end, start

        gene_attrs = ";".join([
                "ID=" + protid,
                "locus_tag=" + protid,
        ])

        cds_attrs = ";".join([
                "ID=" + protid + ".1",
                "Parent=" + protid,
                "locus_tag=" + protid,
                "codon_start=1",
                "transl_table=11",
                "product=" + prottype,
        ])

        src = "GenBank"
        gene = [chrom, src, "gene", start, end, ".", strand, ".", gene_attrs]
        cds = [chrom, src, "CDS", start, end, ".", strand, 0, cds_attrs]

        gff_entries[chrom].append(gene)
        gff_entries[chrom].append(cds)


# read in the annotation
for record in SeqIO.parse(args.annot, "genbank"):
    prev_locus_tag = None
    subfeature_count = 0

    if record.id not in gff_entries:
        gff_entries[record.id] = []

    seqs[record.id] = record

    for feature in record.features:
        if "locus_tag" not in feature.qualifiers:
            continue

        locus_tag = feature.qualifiers["locus_tag"][0]
        is_pseudo = False
        attrs = []
        pid = locus_tag

        if locus_tag != prev_locus_tag:
            subfeature_count = 0
            prev_locus_tag = locus_tag
        else:
            subfeature_count += 1

            if args.tag_ptx:
                parentid = pid

                if locus_tag in identified_refseq_prots:
                    parentid += "_ptx"

                if locus_tag in identified_refseq_exts:
                    parentid += f"_ext_{identified_refseq_exts[locus_tag]}aa"

                if locus_tag in identified_refseq_reds:
                    parentid += f"_red_{identified_refseq_reds[locus_tag]}aa"
            else:
                parentid = pid

            attrs.append(f"Parent={parentid}")

        for key, val in feature.qualifiers.items():
            if key == "translation":
                continue

            if args.no_gene_attr and key == "gene":
                continue

            if key in pseudo_keys:
                key = "pseudo"
                val = ["true"]
                is_pseudo = True

            attrs.append(f"{key}={val[0]}")

        if is_pseudo:
            continue

        for i, part in enumerate(feature.location.parts):
            if i > 0:
                subfeature_count += 1

            if subfeature_count > 0:
                pid = f"{pid}.{subfeature_count}"

            if args.tag_ptx:
                if locus_tag in identified_refseq_prots:
                    pid += "_ptx"

                if locus_tag in identified_refseq_exts:
                    pid += f"_ext_{identified_refseq_exts[locus_tag]}aa"

                if locus_tag in identified_refseq_reds:
                    pid += f"_red_{identified_refseq_reds[locus_tag]}aa"

            id_attrs = [f"ID={pid}"] + attrs

            if feature.type == "CDS":
                if i == 0:
                    phase = 0
                elif part.strand == 1:
                    phase = (3-((part.start-feature.location.start+1) % 3)) % 3
                else:
                    phase = (3-((part.end-feature.location.end+1) % 3)) % 3
            else:
                phase = "."

            strand = "-" if part.strand < 0 else "+"

            gff_entries[record.id].append([
                    record.id,
                    "GenBank",
                    feature.type,
                    part.start + 1,
                    part.end,
                    ".",
                    strand,
                    phase,
                    ";".join(id_attrs)
            ])

with open(args.output, "w") as fo:
    writer = csv.writer(fo, delimiter="\t")

    for chrom, rows in gff_entries.items():
        writer.writerows(sorted(rows, key=lambda row: row[4]))

        if args.fasta_append:
            writer.writerow(["##FASTA"])
            SeqIO.write(seqs[chrom], fo, "fasta")
