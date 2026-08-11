def get_split_cat_input(wildcards):
    basedir = f"results/searches/{wildcards.search}{wildcards.subsearch}/protinf"

    if wildcards.splitcatssfx == "":
        return f"{basedir}/pseudo_overlap/{wildcards.strain}{wildcards.subsample}/proteins.tsv"
    else:
        return f"{basedir}/pepclass_filter/{wildcards.strain}{wildcards.subsample}/proteins.uniq.tsv"


rule protinf_split_categories:
    input:
        get_split_cat_input
    output:
        all = "results/searches/{search}{subsearch}/protinf/split_cats{splitcatssfx}/{strain}{subsample}/all.tsv",
        refseq = "results/searches/{search}{subsearch}/protinf/split_cats{splitcatssfx}/{strain}{subsample}/refseq.tsv",
        pseudo = "results/searches/{search}{subsearch}/protinf/split_cats{splitcatssfx}/{strain}{subsample}/pseudo.tsv",
        starts = "results/searches/{search}{subsearch}/protinf/split_cats{splitcatssfx}/{strain}{subsample}/starts.tsv",
        novels = "results/searches/{search}{subsearch}/protinf/split_cats{splitcatssfx}/{strain}{subsample}/novels.tsv",
        novelty = "results/searches/{search}{subsearch}/protinf/split_cats{splitcatssfx}/{strain}{subsample}/novelty.tsv",
        summary = "results/searches/{search}{subsearch}/protinf/split_cats{splitcatssfx}/{strain}{subsample}/summary.tsv"
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
