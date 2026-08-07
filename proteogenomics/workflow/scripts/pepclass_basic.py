#!/usr/bin/env python3

import argparse
import csv

from Bio import SeqIO

parser = argparse.ArgumentParser(
        description="Identify ambiguous 3b peptides"
)
parser.add_argument(
        "-i", "--iptgxdb", required=True,
        help="iPtgxDB fasta file"
)
parser.add_argument(
        "-p", "--peptides", nargs="+", required=True,
        help="Peptide input files"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="Output file of peptide classes"
)
parser.add_argument(
        "-c", "--contam-prefix", default="sp|",
        help="Contaminant ID prefix (default: 'sp|')"
)
parser.add_argument(
        "-d", "--decoy-prefix", default="rev_",
        help="Decoy ID prefix (default: 'rev_')"
)
parser.add_argument(
        "-e", "--entrapment-suffix", default="_p_target",
        help="Entrapment ID suffix (default: '_p_target')"
)
args = parser.parse_args()

proteins = SeqIO.to_dict(SeqIO.parse(args.iptgxdb, "fasta"))
peptides = {}


def cluster_from_iptgxdb_id(identifier):
    baseid = identifier.removesuffix(args.entrapment_suffix)
    pos = baseid.split("|")[-1].split("_")
    chrom = "_".join(pos[:-5])
    stop = int(pos[-4])
    strand = pos[-3][0]

    return f"{chrom}:{stop}:{strand}"


for file in args.peptides:
    with open(file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            peptide = row["Peptide"]

            if peptide not in peptides:
                peptides[peptide] = []

            if row["Protein"] not in peptides[peptide]:
                peptides[peptide].append(row["Protein"])

            if row["Mapped Proteins"] != "":
                for protein in row["Mapped Proteins"].split(", "):
                    if protein not in peptides[peptide]:
                        peptides[peptide].append(protein)

with open(args.output, "w") as fo:
    writer = csv.writer(fo, delimiter="\t")

    for pep, prots in peptides.items():
        if any([p.startswith(args.contam_prefix) for p in prots]):
            writer.writerow([pep, "contam", ";".join(prots)])
            continue
        if any([p.startswith(args.decoy_prefix + args.contam_prefix) for p in prots]):
            writer.writerow([pep, "contam", ";".join(prots)])
            continue

        pepclass = "unique"

        if len(prots) > 1:
            seqs = [str(proteins[prot].seq) for prot in prots]

            for i in range(0, len(prots)):
                for j in range(i+1, len(prots)):
                    if seqs[i] != seqs[j]:
                        cluster_i = cluster_from_iptgxdb_id(prots[i])
                        cluster_j = cluster_from_iptgxdb_id(prots[j])

                        if cluster_i != cluster_j:
                            pepclass = "ambiguous"
                            break

        writer.writerow([pep, pepclass, ";".join(prots)])