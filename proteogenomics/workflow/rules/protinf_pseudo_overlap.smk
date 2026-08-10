rule protinf_pseudo_overlap:
    input:
        proteins = "results/searches/{search}{subsearch}/protinf/collapse_extensions/{strain}{subsample}/proteins.collapsed.tsv",
        annot = "data/annotations/refseq/{strain}.gff"
    output:
        proteins = "results/searches/{search}{subsearch}/protinf/pseudo_overlap/{strain}{subsample}/proteins.tsv",
        details = "results/searches/{search}{subsearch}/protinf/pseudo_overlap/{strain}{subsample}/details.tsv"
    log:
        "results/searches/{search}{subsearch}/protinf/pseudo_overlap/{strain}{subsample}/log.txt"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/pseudo_overlap.py \
-p '{input.proteins}' \
-a '{input.annot}' \
-o '{output.proteins}' \
-d '{output.details}' \
-e '_p_target' \
2>&1 | tee '{log}'
        """
