rule postproc_integrate_panaroo:
    input:
        panaroo = "results/searches/{search}/postproc/panaroo/output/gene_presence_absence_roary.csv",
        novels = expand("results/searches/{{search}}/postproc/novels/combined/{strain}.tsv", strain=config["strains"])
    output:
        "results/searches/{search}/postproc/novels/novels_panaroo_postproc.tsv"
    params:
        input_flags = lambda wildcards, input: [f"-p '{os.path.splitext(os.path.basename(f))[0]}' '{f}'" for f in input.novels],
        tax_levels = config["postproc"]["conservation_blast"]["tax_levels"]
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/integrate_panaroo.py \
-g {input.panaroo} \
{params.input_flags} \
-t {params.tax_levels} \
-o {output}
        """
