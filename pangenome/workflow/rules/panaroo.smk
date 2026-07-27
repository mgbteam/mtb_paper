rule panaroo_run:
    input:
        lambda wildcards: expand("results/{{strains}}/convert/gff/{strain}.gff", strain=config["strains"][wildcards.strains])
    output:
        folder = directory("results/{strains}/panaroo/run"),
        groups = report("results/{strains}/panaroo/run/gene_presence_absence.csv", category="{strains} panaroo", subcategory="Tables"),
        seqs = "results/{strains}/panaroo/run/pan_genome_reference.fa",
        entropy = "results/{strains}/panaroo/run/alignment_entropy.csv",
        alignments = directory("results/{strains}/panaroo/run/aligned_gene_sequences")
    params:
        flags = config["panaroo"]["flags"]
    threads:
        config["panaroo"]["threads"]
    conda:
        "../envs/panaroo.yml"
    shell:
        """
        panaroo -t {threads} -i {input} -o {output.folder} {params.flags}
        """


rule panaroo_filter_refound:
    input:
        "results/{strains}/panaroo/run/gene_presence_absence.csv"
    output:
        report("results/{strains}/panaroo/run/gene_presence_absence_no_refound.csv", category="{strains} panaroo", subcategory="Tables")
    conda:
        "../envs/panaroo.yml"
    shell:
        """
        workflow/scripts/filter_refound.py -i '{input}' -o '{output}'
        """


rule panaroo_generate_gffs:
    input:
        gffs = lambda wildcards: expand("results/{{strains}}/convert/gff/{strain}.gff", strain=config["strains"][wildcards.strains]),
        panaroo_dir = "results/{strains}/panaroo/run"
    output:
        directory("results/{strains}/panaroo/run/postpanaroo_gffs")
    conda:
        "../envs/panaroo.yml"
    shell:
        """
        panaroo-generate-gffs -i {input.gffs} -o '{input.panaroo_dir}'
        """


rule panaroo_refound_vs_pseudo:
    input:
        gbk = "data/annotation/{strain}.gbff",
        gffdir = "results/{strains}/panaroo/run/postpanaroo_gffs"
    output:
        "results/{strains}/panaroo/refound_vs_pseudo/{strain}.tsv"
    conda:
        "../envs/panaroo.yml"
    shell:
        """
        workflow/scripts/refound_vs_pseudo.py -g '{input.gbk}' -p '{input.gffdir}/{wildcards.strain}_panaroo.gff' -o '{output}'
        """


rule panaroo_plot_refound_vs_pseudo:
    input:
        lambda wildcards: expand("results/{{strains}}/panaroo/refound_vs_pseudo/{strain}.tsv", strain=config["strains"][wildcards.strains])
    output:
        report("results/{strains}/panaroo/refound_vs_pseudo/refound_vs_pseudo.png", category="{strains} panaroo", subcategory="Refinding")
    conda:
        "../envs/panaroo.yml"
    shell:
        """
        workflow/scripts/plot_refound_vs_pseudo.py -i {input} -o '{output}'
        """
