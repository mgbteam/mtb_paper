rule protinf_pseudo_overlap:
    input:
        proteins = "results/searches/{search}{subsearch}/protinf/collapse_extensions/{strain}/proteins.collapsed.tsv",
        annot = "data/annotations/refseq/{strain}.gff"
    output:
        proteins = "results/searches/{search}{subsearch}/protinf/pseudo_overlap/{strain}/proteins.tsv",
        details = "results/searches/{search}{subsearch}/protinf/pseudo_overlap/{strain}/details.tsv"
    log:
        "results/searches/{search}{subsearch}/protinf/pseudo_overlap/{strain}/log.txt"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/pseudo_overlap.py \
-p '{input.proteins}' \
-a '{input.annot}' \
-o '{output.proteins}' \
-d '{output.details}' \
2>&1 | tee '{log}'
        """
