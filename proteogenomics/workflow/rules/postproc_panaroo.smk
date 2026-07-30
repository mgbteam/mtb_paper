rule panaroo_refseq_gff_with_novels:
    input:
        annot = "data/annotations/refseq/{strain}.gb",
        proteins = "results/searches/{search}/protinf/pseudo_overlap/{strain}/proteins.tsv"
    output:
        "results/searches/{search}/postproc/panaroo/input/{strain}.gff"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/create_gff_with_novels.py \
-a {input.annot} \
-p {input.proteins} \
-o {output} \
--tag-ptx \
--fasta-append
        """


rule panaroo_run:
    input:
        expand("results/searches/{{search}}/postproc/panaroo/input/{strain}.gff", strain=config["strains"])
    output:
        folder = directory("results/searches/{search}/postproc/panaroo/output"),
        groups_roary = "results/searches/{search}/postproc/panaroo/output/gene_presence_absence_roary.csv",
        groups = report("results/searches/{search}/postproc/panaroo/output/gene_presence_absence.csv", category="{search}", subcategory="Panaroo")
    params:
        flags = config["postproc"]["panaroo"]["flags"]
    threads:
        config["postproc"]["panaroo"]["threads"]
    conda:
        "../envs/panaroo.yml"
    shell:
        """
panaroo -t {threads} -i {input} -o {output.folder} {params.flags}
        """


rule panaroo_postproc:
    input:
        pangenome = "results/searches/{search}/postproc/panaroo/output/gene_presence_absence_roary.csv",
        iptgxdbs = lambda wildcards: [f"data/iptgxdbs/{get_iptgxdb(wildcards.search)}/{s}/iptgxdb.fasta" for s in config["strains"]]
    output:
        "results/searches/{search}/postproc/panaroo/gene_presence_absence.tsv"
    log:
        "results/searches/{search}/postproc/panaroo/gene_presence_absence.log"
    params:
        refseq_prefix = config["prefixes"]["refseq"]
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/postproc_panaroo.py -p {input.pangenome} -i {input.iptgxdbs} -r {params.refseq_prefix} -o {output} > {log}
        """
