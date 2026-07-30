rule postproc_add_annot:
    input:
        proteins = "results/searches/{search}{subsearch}/protinf/split_categories/{strain}/{file}.tsv",
        annot = "data/annotations/refseq/{strain}.gb"
    output:
        "results/searches/{search}{subsearch}/postproc/add_annot/{strain}/{file}.tsv"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/add_annot.py \
-i '{input.proteins}' \
-a '{input.annot}' \
-o '{output}' \
        """
