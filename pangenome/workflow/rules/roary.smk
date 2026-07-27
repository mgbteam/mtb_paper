rule roary_run:
    input:
        lambda wildcards: expand("results/{{strains}}/convert/gff/{strain}.gff", strain=config["strains"][wildcards.strains])
    output:
        folder = directory("results/{strains}/roary/run"),
        groups = report("results/{strains}/roary/run/gene_presence_absence.csv", category="{strains} roary", subcategory="Tables"),
        seqs = "results/{strains}/roary/run/pan_genome_reference.fa"
    params:
        flags = config["roary"]["flags"]
    threads:
        config["roary"]["threads"]
    conda:
        "../envs/roary.yml"
    shell:
        """
        tmpdir="/tmp/$(uuidgen)"
        roary -p {threads} {params.flags} -f "$tmpdir" {input}

        mkdir -p "{output.folder}"
        mv "$tmpdir"*/* "{output.folder}"
        rm -r "$tmpdir"*
        """
