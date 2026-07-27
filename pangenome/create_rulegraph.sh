#!/bin/sh

snakemake --rulegraph | dot -Tpng > rulegraph.png
snakemake --rulegraph | dot -Tsvg > rulegraph.svg
