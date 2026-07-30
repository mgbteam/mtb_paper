#!/usr/bin/env python3

import argparse
import csv
import sys


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="check overlap with RefSeq pseudogenes"
    )
    parser.add_argument(
        "-p", "--proteins", metavar="TSV", required=True,
        help="list of proteins in tsv format"
    )
    parser.add_argument(
        "-a", "--annot", metavar="GFF", required=True,
        help="annotation in GFF format"
    )
    parser.add_argument(
        "-o", "--output", metavar="TSV", required=True,
        help="output list of reassigned proteins in tsv format")
    parser.add_argument(
        "-d", "--details", metavar="TSV",
        help="output details about pseudogene assignments (optional)"
    )
    return parser.parse_args()


def read_pseudogene_gff(path):
    pseudogenes = []

    with open(path, "r") as infile:
        reader = csv.reader(infile, delimiter="\t")

        for line in reader:
            if line[0].startswith("#"):
                continue

            attrs = [attr.split("=") for attr in line[8].split(";")]
            attr_dict = {attr[0]: attr[1] for attr in attrs}

            if "pseudo" in attr_dict:
                if line[6] == "+":
                    start = int(line[3])
                    end = int(line[4]) + 1
                elif line[6] == "-":
                    start = int(line[4]) + 1
                    end = int(line[3])
                else:
                    print(f"Error: strand in GFF must be '+' or '-', cannot be '{line[6]}'")
                    sys.exit(1)

                pseudogenes.append({"sequence": line[0],
                                    "start": start,
                                    "end": end,
                                    "strand": line[6],
                                    "attributes": attr_dict})

    return pseudogenes


def check_pseudo_overlap(pseudogenes, proteins, output, details):
    detail_rows = []
    reassigned = 0

    with open(proteins, "r") as infile, open(output, "w") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        colnames = reader.fieldnames
        writer = csv.DictWriter(outfile, delimiter="\t", fieldnames=colnames)
        writer.writeheader()

        for row in reader:
            if row["protein-type"] == "RefSeq protein":
                writer.writerow(row)
                continue

            pos_str = row["protein"].split("|")[-1].split("_")
            sequence = "_".join(pos_str[:-5])
            start = int(pos_str[-5])
            end = int(pos_str[-4])
            strand = pos_str[-3][0]
            match = ""
            is_pseudo = False

            if strand == "+":
                end += 1
            else:
                start += 1

            minpos = min(start, end)
            maxpos = max(start, end)

            for pseudo in pseudogenes:
                if pseudo["sequence"] == sequence and pseudo["strand"] == strand:
                    pstart = pseudo["start"]
                    pend = pseudo["end"]

                    pminpos = min(pstart, pend)
                    pmaxpos = max(pstart, pend)

                    if pminpos <= minpos <= pmaxpos or pminpos <= maxpos <= pmaxpos:
                        if pstart == start:
                            if pend == end:
                                is_pseudo = True
                                match = "identical to"
                            else:
                                is_pseudo = True
                                match = "matching start with"
                        elif pend == end:
                            is_pseudo = True
                            match = "matching end with"
                        elif (pstart - start) % 3 == 0 or (pend - end) % 3 == 0:
                            is_pseudo = True
                            match = "in-frame overlap with"
                        else:
                            match = "out-of-frame overlap with"

                        overlap = (min(maxpos, pmaxpos) - max(minpos, pminpos))
                        overlap_rel_pseudo = round(100 * overlap / (pmaxpos - pminpos), 1)
                        overlap_rel_novel = round(100 * overlap / (maxpos - minpos), 1)

                        if overlap == 0:
                            is_pseudo = False
                            match = ""
                            continue

                        if max(overlap_rel_pseudo, overlap_rel_novel) < 50:
                            is_pseudo = False

                        if row["protein-type"] != "RefSeq pseudogene":
                            reassigned += 1

                        row["protein-type"] = "RefSeq pseudogene"

                        match += " pseudogene " + pseudo["attributes"]["ID"]
                        detail_rows.append([
                            row["protein"],
                            minpos,
                            maxpos,
                            pminpos,
                            pmaxpos,
                            strand,
                            match,
                            overlap,
                            overlap_rel_novel,
                            overlap_rel_pseudo,
                            is_pseudo
                        ])

                        break

            writer.writerow(row)

    if details:
        with open(details, "w") as fo:
            writer = csv.writer(fo, delimiter="\t")
            writer.writerow([
                "protein",
                "novel from",
                "novel to",
                "pseudo from",
                "pseudo to",
                "strand",
                "match",
                "overlap (bp)",
                "overlap novel %",
                "overlap pseudo %",
                "pseudo"
            ])
            writer.writerows(detail_rows)

    print(f"Reassigned {reassigned} novel proteins as pseudogenes")


def main():
    args = parse_arguments()
    pseudogenes = read_pseudogene_gff(args.annot)
    check_pseudo_overlap(pseudogenes, args.proteins, args.output, args.details)


if __name__ == "__main__":
    main()
