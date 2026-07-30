rule postproc_operon_mapper:
    input:
        annot = "data/annotations/refseq/{strain}.gb",
        proteins = "results/searches/{search}/protinf/pseudo_overlap/{strain}/proteins.tsv"
    output:
        "results/searches/{search}/postproc/operon_mapper/input/{strain}.gff"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/create_gff_with_novels.py \
-a {input.annot} \
-p {input.proteins} \
-o {output} \
--no-gene-attr \
--tag-ptx
        """
