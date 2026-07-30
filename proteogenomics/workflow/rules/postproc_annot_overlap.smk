rule postproc_annot_overlap:
    input:
        novels = "results/searches/{search}/protinf/split_categories/{strain}/novels.tsv",
        annot = lambda wildcards: "data/annotations/{annot}/" + config["postproc"]["annot_overlap"][wildcards.annot]["filename"]
    output:
        "results/searches/{search}/postproc/novels/annot_overlap/{annot}/{strain}.tsv"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/check_annot_overlap.py -i '{input.novels}' -a '{input.annot}' -n '{wildcards.annot}' -o '{output}'
        """
