#!/usr/bin/env python3

import argparse
import csv

from Bio import SeqIO

import matplotlib.pyplot as plt

from scipy.stats import fisher_exact


def parse_arguments():
    # parse command line arguments
    parser = argparse.ArgumentParser(
            description="Signal of selection prediction based on G/C skew"
    )
    parser.add_argument(
            "-i", "--input", required=True,
            help="Table of protiens in TSV format"
    )
    parser.add_argument(
            "-g", "--genome", required=True,
            help="Genome in GenBank format"
    )
    parser.add_argument(
            "-o", "--output", required=True,
            help="Output file"
    )
    parser.add_argument(
            "-p", "--plot",
            help="Optionally output a plot"
    )
    return parser.parse_args()


def get_non_overlapping_novel_seq(genome, identifier):
    # extract the location of the novel ORF
    loc = identifier.split("|")[-1].split("_")
    chrom = "_".join(loc[:-5])
    start = int(loc[-5])
    end = int(loc[-4])
    strand = loc[-3][0]

    # get the DNA sequence of the novel ORF without the start and stop codon
    if strand == "+":
        seq = genome[chrom].seq[start-1+3:end-3]
    else:
        seq = genome[chrom].seq[end-1+3:start-3]
        start, end = end, start

    untrimmed_len = len(seq)

    # check if the location overlaps with a known ORF
    overlap_type = "none"
    overlap_genes = []
    overlap_annots = []

    for feature in genome[chrom].features:
        if feature.type == "CDS":
            if "pseudo" in feature.qualifiers:
                continue

            start_known = feature.location.start
            end_known = feature.location.end
            strand_known = "-" if feature.location.strand < 0 else "+"

            # if the novel ORF starts within a known ORF
            if start_known <= start <= end_known:
                overlap_genes.append(feature.qualifiers["locus_tag"][0])
                overlap_annots.append(feature.qualifiers["product"][0])

                # if the novel ORF also ends within the known ORF
                if start_known <= end <= end_known:
                    # discard the novel ORF
                    if strand != strand_known:
                        overlap_type = "antisense"
                    elif (end_known - end) % 3 == 0:
                        overlap_type = "in-frame"
                    else:
                        overlap_type = "out-of-frame"

                    seq = ""
                    break

                border = end_known - start

                # determine overlap type
                if strand != strand_known:
                    overlap_type = "antisense"
                elif border % 3 != 0:
                    overlap_type = "out-of-frame"
                else:
                    overlap_type = "in-frame"

                # round to complete codons
                if border % 3 != 0:
                    border = border - border % 3 + 3

                # discard the overlapping region
                seq = seq[border:]
                start = start + border

            # if the novel ORF only ends within a known ORF
            elif start_known <= end <= end_known:
                overlap_genes.append(feature.qualifiers["locus_tag"][0])
                overlap_annots.append(feature.qualifiers["product"][0])
                border = start_known - start

                # determine overlap type
                if strand != strand_known:
                    overlap_type = "antisense"
                elif border % 3 != 0:
                    overlap_type = "out-of-frame"
                else:
                    overlap_type = "in-frame"

                # round to complete codons
                if border % 3 != 0:
                    border = border - border % 3

                # discard the overlapping region
                seq = seq[:border]
                end = start + border

            # if a known ORF is completely within the novel ORF
            elif start <= start_known and end >= end_known:
                overlap_genes.append(feature.qualifiers["locus_tag"][0])
                overlap_annots.append(feature.qualifiers["product"][0])
                border_left = start_known - start
                border_right = end_known - start

                # determine overlap type
                if strand != strand_known:
                    overlap_type = "antisense"
                elif border_left % 3 != 0 or border_right % 3 != 0:
                    overlap_type = "out-of-frame"
                else:
                    overlap_type = "in-frame"

                # round to complete codons
                if border_left % 3 != 0:
                    border_left = border_left - border_left % 3

                if border_right % 3 != 0:
                    border_right = border_right - border_right % 3 + 3

                # extract the non-overlapping regions left and right
                seq = seq[:border_left] + seq[border_right:]

                # stop the analyis of this split novel ORF as it
                # can lead to issues in rare cases with more overlaps
                break

    overlap_percent = 1 - len(seq) / untrimmed_len

    if strand == "-" and seq != "":
        seq = seq.reverse_complement()

    return overlap_type, overlap_percent, overlap_genes, overlap_annots, seq


def count_nucs_per_codon_pos(seq):
    # loop over codons and count the occurrence of each nucleotide per position
    nuc_counts = {i: {"GC": 0, "AT": 0} for i in range(3)}
    number_of_codons = 0

    for i in range(0, len(seq), 3):
        number_of_codons += 1
        codon = seq[i:i+3]

        for j in range(3):
            if codon[j] in ["G", "C"]:
                nuc_counts[j]["GC"] += 1
            elif codon[j] in ["A", "T"]:
                nuc_counts[j]["AT"] += 1
            else:
                print(f"Unknown nucleotide: {codon[j]}")

    return number_of_codons, nuc_counts


def calc_gc_skew_and_significance(nuc_counts):
    # calculate the G/C skew
    if nuc_counts[1]["GC"] == 0:
        gc_skew = 0
    else:
        gc_skew = float(nuc_counts[2]["GC"]) / nuc_counts[1]["GC"]

    # calculate Fisher's exact test for the G/C frequency
    counts = [[nuc_counts[1]["GC"], nuc_counts[1]["AT"]],
              [nuc_counts[2]["GC"], nuc_counts[2]["AT"]]]

    _, p_value = fisher_exact(counts, alternative="two-sided")

    return gc_skew, p_value


def calc_average_gc_frequencies(nuc_counts_total):
    # calculate the average G/C frequency for each position in the codon
    gc_freqs = []

    for i in range(3):
        total_gc = nuc_counts_total[i]["GC"]
        total_nt = nuc_counts_total[i]["GC"] + nuc_counts_total[i]["AT"]
        gc_freqs.append(total_gc / total_nt)

    return gc_freqs


def count_gc_freqs_novels(genome, inputfile, outputfile):
    # load novel ORFs and extract GC skew for each ORF
    outcols = [
        "protein",
        "Overlap",
        "Overlap Type",
        "Overlap Genes",
        "Overlap Annots",
        "Non-Overlapping nt-Sequence",
        "G/C Position 1",
        "G/C Position 2",
        "G/C Position 3",
        "# Codons",
        "G/C Skew",
        "G/C Skew p-value",
        "Selection Signal",
    ]

    nuc_counts_total = {i: {"GC": 0, "AT": 0} for i in range(3)}

    with open(inputfile, "r") as infile, open(outputfile, "w") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        writer = csv.DictWriter(outfile, delimiter="\t", fieldnames=outcols)
        writer.writeheader()
        novels_total = 0
        novels_included = 0

        for row in reader:
            # get the sequence of the novel that doesn't overlap a RefSeq gene
            novels_total += 1
            oltype, olperc, olgenes, olannots, seq = get_non_overlapping_novel_seq(genome, row["protein"])

            if len(seq) > 0:
                novels_included += 1

            # count the nucleotides in non-overlapping sequence, update totals
            codon_count, nuc_counts = count_nucs_per_codon_pos(seq)

            for i in range(3):
                for key in ["GC", "AT"]:
                    nuc_counts_total[i][key] += nuc_counts[i][key]

            # calculate GC skew and significance, write results to output table
            gc_skew, p_value = calc_gc_skew_and_significance(nuc_counts)

            writer.writerow({
                "protein": row["protein"],
                "Overlap": olperc,
                "Overlap Type": oltype,
                "Overlap Genes": ";".join(olgenes),
                "Overlap Annots": ";".join(olannots),
                "Non-Overlapping nt-Sequence": seq,
                "G/C Position 1": nuc_counts[0]["GC"],
                "G/C Position 2": nuc_counts[1]["GC"],
                "G/C Position 3": nuc_counts[2]["GC"],
                "# Codons": codon_count,
                "G/C Skew": gc_skew,
                "G/C Skew p-value": p_value,
                "Selection Signal": p_value < 0.1 and gc_skew > 1,
            })

    # print total number of novel ORFs processed
    print(f"Total novel ORFs: {novels_total}")
    print(f"Novel ORFs included: {novels_included}")

    # calculate average G/C frequency for each position in the novel codons
    return calc_average_gc_frequencies(nuc_counts_total)


def count_gc_freqs_refseq(genome):
    # calculate the average G/C frequeny for all refseq codons
    nuc_counts_total = {i: {"GC": 0, "AT": 0} for i in range(3)}

    for record_id, record in genome.items():
        for feature in record.features:
            if feature.type != "CDS" or "pseudo" in feature.qualifiers:
                continue

            # extract the sequence without the start and stop codon
            seq = feature.location.extract(record).seq[3:-3]

            # count the nucleotides in sequence, update totals
            codon_count, nuc_counts = count_nucs_per_codon_pos(seq)

            for i in range(3):
                for key in ["GC", "AT"]:
                    nuc_counts_total[i][key] += nuc_counts[i][key]

    # calculate average G/C frequency for each position in the refseq codons
    return calc_average_gc_frequencies(nuc_counts_total)


def main():
    args = parse_arguments()
    genome = SeqIO.to_dict(SeqIO.parse(args.genome, "genbank"))
    avg_gc_freq_novels = count_gc_freqs_novels(genome, args.input, args.output)

    if args.plot:
        colors_novel = ["#ff8e1c", "#c76301", "#713900"]
        colors_refseq = ["#7071ff", "#1d1cfe", "#000070"]
        avg_gc_freq_refseq = count_gc_freqs_refseq(genome)

        fig, (ax1, ax2) = plt.subplots(1, 2)

        ax1.bar(range(3), avg_gc_freq_refseq, color=colors_refseq)
        ax1.set_xticks(range(3), ["1st", "2nd", "3rd"])
        ax1.set_xlabel("Position in codon")
        ax1.set_ylabel("G/C frequency")
        ax1.set_ylim(0, 1.0)
        ax1.set_title("RefSeq")

        ax2.bar(range(3), avg_gc_freq_novels, color=colors_novel)
        ax2.set_xticks(range(3), ["1st", "2nd", "3rd"])
        ax2.set_xlabel("Position in codon")
        ax2.set_ylim(0, 1.0)
        ax2.set_title("Novel")

        fig.savefig(args.plot)


if __name__ == "__main__":
    main()
