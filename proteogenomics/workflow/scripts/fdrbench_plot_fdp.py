#!/usr/bin/env python3

import os
import sys

import matplotlib.pyplot as plt

import pandas as pd

if len(sys.argv) != 3:
    print("Usage: fdrbench_plot_fdp.py <input_file> <output_file>")
    sys.exit(1)

plt.figure()

if os.path.getsize(sys.argv[1]) > 0:
    data = pd.read_csv(sys.argv[1])
    data = data[data["q_value"] <= 0.1]

    data["q_value"] = 100 * data["q_value"]
    data["combined_fdp"] = 100 * data["combined_fdp"]
    data["paired_fdp"] = 100 * data["paired_fdp"]
    data["lower_bound_fdp"] = 100 * data["lower_bound_fdp"]

    data.plot(x="q_value", y=["combined_fdp", "paired_fdp", "lower_bound_fdp"], kind="line")
    plt.plot(data["q_value"], data["q_value"], color="lightgrey")
    plt.axvline(x=1, color="blue", linestyle="--")
    plt.xlabel("FDR threshold (%)")
    plt.ylabel("Estimated FDP (%)")
    plt.title("Estimated FDP vs FDR threshold")

plt.savefig(sys.argv[2])
