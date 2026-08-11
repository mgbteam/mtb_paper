rule entrapment_plot_fdp:
    input:
        "results/searches/{search}/entrapment/fdrbench/{stage}/calc_fdp/{strain}/{sample}/{subset}_fdp.csv"
    output:
        report("results/searches/{search}/entrapment/fdrbench/{stage}/plot_fdp/{strain}/{sample}/{subset}_fdp.png", category="{search}", subcategory="FDRBench", labels={"Stage": "{stage}", "Strain": "{strain}", "Subset": "{subset}", "Sample": "{sample}"})
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/fdrbench_plot_fdp.py '{input}' '{output}'
        """
