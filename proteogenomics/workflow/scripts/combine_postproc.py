#!/usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser(
        description="Combine tables from all post-processing analyses"
)
parser.add_argument(
        "-t", "--table", required=True,
        help="TSV table of identified proteins"
)
parser.add_argument(
        "-f", "--fdpfilter",
        help="FDP filtering results based on FDRBench"
)
parser.add_argument(
        "-a", "--annotoverlap", nargs="+",
        help="Annotation overlap TSV results"
)
parser.add_argument(
        "-p", "--phyre2",
        help="Phyre2 results (summaryinfo.txt)"
)
parser.add_argument(
        "-e", "--eggnog",
        help="eggNOG results (*.emapper.annotations)"
)
parser.add_argument(
        "-i", "--interpro",
        help="InterProScan results in TSV format"
)
parser.add_argument(
        "-s", "--psortb",
        help="PsortB results in txt format"
)
parser.add_argument(
        "-l", "--lipop",
        help="LipoP results in GFF format"
)
parser.add_argument(
        "-g", "--codongcfreq",
        help="Codon GC frequency analysis results in TSV format"
)
parser.add_argument(
        "-H", "--hhsuite",
        help="hhsuite results summary"
)
parser.add_argument(
        "-S", "--ampscanner2",
        help="ampscanner2 results in csv format"
)
parser.add_argument(
        "-c", "--conservationblast",
        help="Conservation BLAST results summary"
)
parser.add_argument(
        "-m", "--operonmapper",
        help="Folder with OperonMapper results (remove suffix from files)"
)
parser.add_argument(
        "-A", "--alphafold",
        help="Confidence ranking from alphafold"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="Output combined table"
)
args = parser.parse_args()
combined = {}

# get novels from protein table
with open(args.table, "r") as fi:
    reader = csv.DictReader(fi, delimiter="\t")
    colnames = reader.fieldnames

    for row in reader:
        if any([s in row["protein-type"] for s in ["pseudogene", "RefSeq"]]):
            continue

        identifier = row["protein"].replace("|", "_")
        combined[identifier] = row

novels_present = len(combined) > 0

# add FDP filter
if args.fdpfilter and novels_present:
    with open(args.fdpfilter, "r") as fi:
        reader = csv.DictReader(fi, delimiter="\t")
        colnames.append("FDP Filter")

        for row in reader:
            identifier = row["protein"].replace("|", "_")
            combined[identifier]["FDP Filter"] = row["FDP Filter"]

# Add annotation overlap
if args.annotoverlap and novels_present:
    for file in args.annotoverlap:
        with open(file, "r") as fi:
            reader = csv.DictReader(fi, delimiter="\t")

            for row in reader:
                identifier = row["protein"].replace("|", "_")
                combined[identifier].update(row)

            colnames.extend([c for c in row if c != "protein"])

# add conservation BLAST
if args.conservationblast and novels_present:
    with open(args.conservationblast, "r") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        colnames.extend(reader.fieldnames[1:])

        for row in reader:
            identifier = row.pop("Protein").replace("|", "_")
            combined[identifier].update(row)

# add Phyre2
if args.phyre2 and novels_present:
    colnames.extend(["Phyre2 Confidence (%)", "Phyre2 Hit info 2"])

    with open(args.phyre2, "r") as fi:
        header = next(fi).strip().removeprefix("# ").split(" | ")

        for line in fi:
            row = {k: v for (k, v) in zip(header, line.strip().removesuffix(" |").split(" | "))}
            identifier = row["Description"].replace("|", "_")

            if "Confidence (%)" in row and row["Confidence (%)"] is not None:
                if float(row["Confidence (%)"]) < 30:
                    continue

                combined[identifier]["Phyre2 Confidence (%)"] = row["Confidence (%)"]

            if "Hit info 2" in row and row["Hit info 2"] is not None:
                combined[identifier]["Phyre2 Hit info 2"] = row["Hit info 2"].split(":")[-1]

# add eggNOG
if args.eggnog and novels_present:
    eggnog_cols = ["evalue", "score", "Description", "Preferred_name"]
    colnames.extend(["eggNOG " + col for col in eggnog_cols])

    with open(args.eggnog, "r") as fi:
        reader = csv.DictReader((l for l in fi if not l.startswith("##")), delimiter="\t")

        for row in reader:
            identifier = row["#query"].replace("|", "_")

            if identifier in combined:
                combined[identifier].update({"eggNOG " + col: row[col] for col in eggnog_cols})

# add hhsuite
if args.hhsuite and novels_present:
    with open(args.hhsuite, "r") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        colnames.extend(reader.fieldnames[1:])

        for row in reader:
            identifier = row.pop("protein").replace("|", "_")
            combined[identifier].update(row)


# add AMP Scanner Version 2
if args.ampscanner2 and novels_present:
    colnames.extend(["AMP2 Class", "AMP2 Probability", "Phobius"])

    with open(args.ampscanner2, "r") as infile:
        reader = csv.DictReader(infile, delimiter=",")

        for row in reader:
            identifier = row["SeqID"].replace("|", "_")

            if identifier.endswith("*"):
                continue

            combined[identifier]["AMP2 Class"] = row["Prediction_Class"]
            combined[identifier]["AMP2 Probability"] = row["Prediction_Probability"]


# add InterProScan
if args.interpro and novels_present:
    colnames.extend(["TMHMM", "SignalP", "Phobius"])

    with open(args.interpro, "r") as fi:
        reader = csv.reader(fi, delimiter="\t")

        for row in reader:
            identifier = row[0].replace("|", "_")

            if row[3] == "TMHMM":
                combined[identifier]["TMHMM"] = "true"
            elif row[3] == "SignalP_GRAM_POSITIVE":
                combined[identifier]["SignalP"] = "gram+"
            elif row[3] == "SignalP_GRAM_NEGATIVE":
                combined[identifier]["SignalP"] = "gram-"
            elif row[3] == "Phobius":
                if "Phobius" not in combined[identifier]:
                    combined[identifier]["Phobius"] = set()

                if row[4].startswith("SIGNAL"):
                    combined[identifier]["Phobius"].add("signal")
                elif row[4].startswith("CYTOPLASMIC"):
                    combined[identifier]["Phobius"].add("cytoplasmic")
                elif row[4].startswith("NON_CYTOPLASMIC"):
                    combined[identifier]["Phobius"].add("non-cytoplasmic")
                elif row[4].startswith("TRANSMEMBRANE"):
                    combined[identifier]["Phobius"].add("transmembrane")

        for row in combined.values():
            if "Phobius" in row:
                row["Phobius"] = ";".join(row["Phobius"])

# add LipoP
if args.lipop and novels_present:
    colnames.extend(["LipoP", "LipoP Score"])

    with open(args.lipop, "r") as fi:
        reader = csv.reader(fi, delimiter="\t")

        for row in reader:
            if row[0].startswith("#") or not row[1].endswith("Best"):
                continue

            if float(row[5]) > 0:
                identifier = row[0].replace("|", "_")
                combined[identifier]["LipoP"] = row[2]
                combined[identifier]["LipoP Score"] = row[5]

# add PsortB
if args.psortb and novels_present:
    colnames.extend(["PSORTb Localization", "PSORTb Score"])

    with open(args.psortb, "r") as infile:
        line = next(infile, None)

        while line is not None:
            if line.startswith("SeqID"):
                identifier = line.rstrip().split(": ", 1)[1].replace("|", "_")
            elif line.startswith("  Final Prediction:"):
                line = next(infile).lstrip().rstrip()

                if line.startswith("Unknown"):
                    (location, score) = ("", "")
                else:
                    (location, score) = line.rsplit(None, 1)

                combined[identifier]["PSORTb Localization"] = location
                combined[identifier]["PSORTb Score"] = score

            line = next(infile, None)

# add codon GC frequency analysis results
if args.codongcfreq and novels_present:
    colnames.extend(["Overlap", "Overlap Type", "Overlap Genes", "Overlap Annots", "Selection"])

    with open(args.codongcfreq, "r") as infile:
        reader = csv.DictReader(infile, delimiter="\t")

        for row in reader:
            identifier = row["protein"].replace("|", "_")
            selection = "Yes" if row["Selection Signal"] == "True" else ""
            combined[identifier]["Overlap"] = row["Overlap"]
            combined[identifier]["Overlap Type"] = row["Overlap Type"]
            combined[identifier]["Overlap Genes"] = row["Overlap Genes"]
            combined[identifier]["Overlap Annots"] = row["Overlap Annots"]
            combined[identifier]["Selection"] = selection

# add OperonMapper results
if args.operonmapper and novels_present:
    colnames.extend(["OperonMapper COG", "OperonMapper BLASTP", "Operon Genes", "Operon Functions"])

    with open(f"{args.operonmapper}/list_of_operons", "r") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        current_operon_id = None
        current_operon = []
        operons = {}

        for row in reader:
            if row["Operon"] != "":
                if len(current_operon) > 1 and current_operon_id is not None:
                    operons[current_operon_id] = current_operon

                current_operon_id = int(row["Operon"])
                current_operon = {}
            else:
                function = row["Function"].split("] ")[-1].strip()
                current_operon[row["IdGene"]] = function
                identifier = ".".join(row["IdGene"].split(".")[:-1]).replace("|", "_")

                if identifier in combined:
                    if function != "NA":
                        combined[identifier]["OperonMapper COG"] = function

        for operon_id, operon in operons.items():
            for gene1 in operon:
                identifier = ".".join(gene1.split(".")[:-1]).replace("|", "_")

                if identifier in combined:
                    genes = []
                    functions = []

                    for gene2, function2 in operon.items():
                        if gene2 != gene1:
                            genes.append(gene2)
                            functions.append(function2)

                    combined[identifier]["Operon Genes"] = ";".join(genes)
                    combined[identifier]["Operon Functions"] = ";".join(functions)

    with open(f"{args.operonmapper}/functional_descriptions", "r") as infile:
        reader = csv.DictReader(infile, delimiter="\t")

        for row in reader:
            identifier = ".".join(row["#gene_id"].split(".")[:-1]).replace("|", "_")

            if identifier in combined:
                combined[identifier]["OperonMapper BLASTP"] = row["product"]

# add alphafold confidence
if args.alphafold and novels_present:
    with open(args.alphafold, "r") as fi:
        reader = csv.reader(fi, delimiter="\t")
        colnames.append("alphafold confidence")

        for row in reader:
            identifier = row[1].replace("|", "_")
            combined[identifier]["alphafold confidence"] = row[0]

# write combined table
with open(args.output, "w") as fo:
    at_end = ["MW", "Instability", "Aromaticity", "Isoelectric point"]
    remove = at_end + ["protein-group", "count proteins", "protein-type", "protein-reference-id", "total - spectra", "total - abundance"]
    colnames = [e for e in colnames if e not in remove]
    colnames.extend(at_end)
    writer = csv.DictWriter(fo, delimiter="\t", fieldnames=colnames, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(combined.values())
