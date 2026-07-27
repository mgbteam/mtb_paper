#!/usr/bin/env python3

import argparse
import csv
import os

import matplotlib.pyplot as plt

import numpy as np

parser = argparse.ArgumentParser(
        description="Create a barplot from refound vs. pseudogene results"
)
parser.add_argument(
        "-i", "--input", nargs="+", required=True,
        help="TSV files containing refound vs. pseudogene statistics"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="Barplot output"
)
args = parser.parse_args()

counts = {
        "Pseudo: Same Start/Stop": [],
        "Pseudo: In Frame": [],
        "Novel": [],
}

colors = ["silver", "tab:blue", "tab:red"]
labels = []

for file in args.input:
    labels.append(os.path.splitext(os.path.basename(file))[0])

    with open(file, "r") as fi:
        reader = csv.reader(fi, delimiter="\t")

        for row in reader:
            if row[0] in counts:
                counts[row[0]].append(int(row[1]))

plt.figure()
bottom = np.zeros(len(labels))

for (category, values), color in zip(counts.items(), colors):
    plt.bar(labels, values, label=category, bottom=bottom, color=color)
    bottom += np.array(values)

plt.title("Panaroo refound genes vs. annotated pseudogenes")
plt.legend()
plt.savefig(args.output)
