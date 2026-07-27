rule upsetplots_regex:
    input:
        "results/{strains}/{tool}/run/gene_presence_absence.csv"
    output:
        report("results/{strains}/{tool}/upsetplots/{subset}_regex/{plot}.png", category="{strains} {tool}", subcategory="UpSetPlots {subset}")
    params:
        color_flags = create_color_flags,
        strains_regex = lambda wildcards: config["upsetplots"]["subsets_regex"][f"{wildcards.subset}"],
        plot_flags = lambda wildcards: config["upsetplots"]["plots"][f"{wildcards.plot}"]
    conda:
        "../envs/upsetplot.yml"
    shell:
        """
        workflow/scripts/create_upsetplot.py -i '{input}' -o '{output}' -s '{params.strains_regex}' {params.color_flags} {params.plot_flags}
        """


rule upsetplots_list:
    input:
        "results/{strains}/{tool}/run/gene_presence_absence.csv"
    output:
        report("results/{strains}/{tool}/upsetplots/{subset}_list/{plot}.png", category="{strains} {tool}", subcategory="UpSetPlots {subset}")
    params:
        color_flags = create_color_flags,
        strains_list = lambda wildcards: config["upsetplots"]["subsets_list"][f"{wildcards.subset}"],
        plot_flags = lambda wildcards: config["upsetplots"]["plots"][f"{wildcards.plot}"]
    conda:
        "../envs/upsetplot.yml"
    shell:
        """
        workflow/scripts/create_upsetplot.py -i '{input}' -o '{output}' -l '{params.strains_list}' {params.color_flags} {params.plot_flags}
        """
