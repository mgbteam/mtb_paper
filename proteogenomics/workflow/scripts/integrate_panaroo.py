#!/usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser(
        description="Merge post-processing results by pan-genome clusters"
)
parser.add_argument(
        "-g", "--gene-presence-absence", required=True,
        help="Gene presence/absence table in Roary format"
)
parser.add_argument(
        "-p", "--proteins", nargs=2, action="append", required=True,
        help="Column name in roary output and corresponding protein table"
)
parser.add_argument(
        "-t", "--taxa-levels", nargs="+", required=True,
        help="Taxonomy levels used in conservation BLAST"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="Output file in TSV format"
)
args = parser.parse_args()

unmerged_cols = {
        "total - psms": "PSMs",
        "total - peptides": "Peps",
        "FDP Filter": "FDP",
}

unique_val_cols = {
        "sequence length": "Length",
        "Phyre2 Confidence (%)": "Phyre2 Conf",
        "Phyre2 Hit info 2": "Phyre2 Description",
        "hhsuite probability": "hhsuite probability",
        "hhsuite description": "hhsuite description",
        "eggNOG evalue": "eggNOG evalue",
        "eggNOG score": "eggNOG score",
        "eggNOG Description": "eggNOG Description",
        "AMP2 Class": "AMP2 Class",
        "AMP2 Probability": "AMP2 Probability",
        "TMHMM": "TMHMM",
        "SignalP": "SignalP",
        "LipoP Score": "LipoP Score",
        "LipoP": "LipoP",
        "PSORTb Score": "PSORTb Score",
        "PSORTb Localization": "PSORTb Localization",
        "OperonMapper COG": "OperonMapper COG",
        "OperonMapper BLASTP": "OperonMapper BLASTP",
        "Selection": "Selection",
        "MW": "MW",
        "Instability": "Instability",
        "Aromaticity": "Aromaticity",
        "Isoelectric point": "Isoelectric point",
        "Annotation": "Annotation",
}

unique_member_cols = {
        "Phobius": "Phobius",
        "Operon Functions": "Operon Gene Functions",
        "Overlap": "Overlap",
        "Overlap Type": "Overlap Type",
        "Overlap Annots": "Overlap Annots",
}

# read post-processing results of novel proteins for each strain
prot_info = {}

for strain, file in args.proteins:
    prot_info[strain] = {}

    with open(file, "r") as fi:
        reader = csv.DictReader(fi, delimiter="\t")

        for row in reader:
            refseq_cols = [col for col in row if col.startswith("refseq")]

            for col in refseq_cols:
                if row[col] != "":
                    row[col] = "Yes"
                if col not in unique_val_cols:
                    unique_val_cols[col] = col

            prot_info[strain][row["protein"]] = row

# read gene presence/absence pan-genome grouping from Roary or Panaroo
groups = {}
outcols = ["Novel", "Strains", "RefSeq in other strain"]

with open(args.gene_presence_absence, "r") as fi:
    reader = csv.DictReader(fi, delimiter=",")

    for row in reader:
        annot = row["Annotation"]

        if "novel" not in annot:
            continue

        ids = {f"ID {s}": ".".join(row[s].split(".")[:-1]) for s in prot_info}
        groups[row["Gene"]] = ids

        # Remove mentions of novel from annotation column
        if annot != "novel":
            annot = annot.replace(";novel", "").replace("novel;", "")
            groups[row["Gene"]]["RefSeq in other strain"] = annot
        else:
            groups[row["Gene"]]["RefSeq in other strain"] = ""

# Consolidate post-processing results for each pan-genome group
for i, group in enumerate(groups):
    # get the protein info per strain
    info = {}

    for strain in prot_info:
        prot_id = groups[group][f"ID {strain}"]
        info[strain] = prot_info.get(strain, {}).get(prot_id, {})

    # process unmerged columns
    for oldcol, newcol in unmerged_cols.items():
        for strain in prot_info:
            groups[group][f"{newcol} {strain}"] = info[strain].get(oldcol, "")

            if i == 0:
                outcols.append(f"{newcol} {strain}")

    # process unique value columns
    for oldcol, newcol in unique_val_cols.items():
        vals = set()

        for strain in prot_info:
            vals.add(info[strain].get(oldcol, ""))

        if i == 0:
            outcols.append(newcol)

        if "" in vals:
            vals.remove("")

        groups[group][newcol] = ";".join(vals)

    # process unique member columns (members separated by ';')
    for oldcol, newcol in unique_member_cols.items():
        vals = set()

        for strain in prot_info:
            vals.update(info[strain].get(oldcol, "").split(";"))

        if i == 0:
            outcols.append(newcol)

        if "" in vals:
            vals.remove("")

        groups[group][newcol] = ", ".join(vals)

    # process BLAST subtaxa matches as min value columns
    for col in args.taxa_levels:
        vals = []

        for strain in prot_info:
            val_str = info[strain].get(col, "")

            if val_str != "":
                vals.append(int(val_str))

        if i == 0:
            outcols.append(col)

        if vals:
            groups[group][col] = max(vals)
        else:
            groups[group][col] = ""

outcols.extend([f"ID {s}" for s in prot_info])

# sort by number of strains the novel protein was seen in
for group in groups:
    strains = 0

    for strain in prot_info:
        if groups[group][f"PSMs {strain}"] != "":
            strains += 1

    groups[group]["Strains"] = strains

outrows = sorted(groups.values(), key=lambda row: row["Strains"], reverse=True)

# number the novels
for i in range(len(outrows)):
    outrows[i]["Novel"] = f"novel_{i+1}"

# remove all empty columns from output
outcols_final = []

for col in outcols:
    if any([row[col] != "" for row in outrows]):
        outcols_final.append(col)

with open(args.output, "w") as fo:
    writer = csv.DictWriter(fo, delimiter="\t", fieldnames=outcols_final, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(outrows)
