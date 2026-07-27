#!/bin/sh

snakemake --cores $(nproc) --use-conda --printshellcmds "$@"
