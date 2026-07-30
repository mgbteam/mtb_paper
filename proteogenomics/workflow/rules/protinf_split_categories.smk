checkpoint protinf_split_categories:
    input:
        "results/searches/{search}{subsearch}/protinf/pseudo_overlap/{strain}/proteins.tsv"
    output:
        refseq = "results/searches/{search}{subsearch}/protinf/split_categories/{strain}/refseq.tsv",
        pseudo = "results/searches/{search}{subsearch}/protinf/split_categories/{strain}/pseudo.tsv",
        starts = "results/searches/{search}{subsearch}/protinf/split_categories/{strain}/starts.tsv",
        novels = "results/searches/{search}{subsearch}/protinf/split_categories/{strain}/novels.tsv",
        summary = "results/searches/{search}{subsearch}/protinf/split_categories/{strain}/summary.tsv"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/split_categories.py \
-i '{input}' \
-r '{output.refseq}' \
-p '{output.pseudo}' \
-s '{output.starts}' \
-n '{output.novels}' \
> '{output.summary}'
        """
