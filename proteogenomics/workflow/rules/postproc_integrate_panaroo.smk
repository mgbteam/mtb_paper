def collect_combined_postproc(wildcards):
    folder = f"results/searches/{wildcards.search}"

    for strain in config["strains"]:
        with open(f"{folder}/protinf/split_categories/{strain}/novels.tsv", "r") as fi:
            if len([True for line in fi]) > 1:
                yield f"{folder}/postproc/novels/combined/{strain}.tsv"


def collect_combined_postproc_flags(wildcards):
    for file in collect_combined_postproc(wildcards):
        strain = os.path.splitext(os.path.basename(file))[0]
        yield f"-p {strain} {file}"


rule postproc_integrate_panaroo:
    input:
        panaroo = "results/searches/{search}/postproc/panaroo/output/gene_presence_absence_roary.csv",
        novels = collect_combined_postproc
    output:
        "results/searches/{search}/postproc/novels/novels_panaroo_postproc.tsv"
    params:
        input_flags = collect_combined_postproc_flags,
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
