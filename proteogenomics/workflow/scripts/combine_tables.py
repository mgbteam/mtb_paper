#!/usr/bin/env python3

import argparse

import pandas as pd

parser = argparse.ArgumentParser(description="Join TSV files by common column")
parser.add_argument(
        "-l", "--left-column", required=True,
        help="Left column to combine by"
)
parser.add_argument(
        "-r", "--right-column", required=True,
        help="Right column to combine by"
)
parser.add_argument(
        "-s", "--separator", default="\t",
        help="Separator"
)
parser.add_argument(
        "-m", "--method", default="outer",
        help="left, right, [outer], inner, cross"
)
parser.add_argument(
        "-o", "--output", required=True,
        help="Output file"
)
parser.add_argument(
        "input", nargs="+",
        help="Input files"
)
args = parser.parse_args()

merged = pd.read_csv(args.input[0], sep=args.separator)

for f in args.input[1:]:
    df = pd.read_csv(f, sep=args.separator)
    merged = pd.merge(merged, df, left_on=args.left_column, right_on=args.right_column, how=args.method)

    if args.left_column != args.right_column:
        merged.drop(args.right_column, axis=1, inplace=True)

merged.to_csv(args.output, sep=args.separator, index=False)
