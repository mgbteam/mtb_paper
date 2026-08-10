rule protinf_split_categories:
    input:
        "results/searches/{search}{subsearch}/protinf/pseudo_overlap/{strain}{subsample}/proteins.tsv"
    output:
        all = "results/searches/{search}{subsearch}/protinf/split_categories/{strain}{subsample}/all.tsv",
        refseq = "results/searches/{search}{subsearch}/protinf/split_categories/{strain}{subsample}/refseq.tsv",
        pseudo = "results/searches/{search}{subsearch}/protinf/split_categories/{strain}{subsample}/pseudo.tsv",
        starts = "results/searches/{search}{subsearch}/protinf/split_categories/{strain}{subsample}/starts.tsv",
        novels = "results/searches/{search}{subsearch}/protinf/split_categories/{strain}{subsample}/novels.tsv",
        novelty = "results/searches/{search}{subsearch}/protinf/split_categories/{strain}{subsample}/novelty.tsv",
        summary = "results/searches/{search}{subsearch}/protinf/split_categories/{strain}{subsample}/summary.tsv"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/split_categories.py \
-i '{input}' \
-a '{output.all}' \
-r '{output.refseq}' \
-p '{output.pseudo}' \
-s '{output.starts}' \
-n '{output.novels}' \
-N '{output.novelty}' \
> '{output.summary}'
        """
