#!/bin/sh

snakemake --cores $(nproc) --use-conda --conda-frontend conda --printshellcmds "$@"
