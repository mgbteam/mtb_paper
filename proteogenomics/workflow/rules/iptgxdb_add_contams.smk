rule iptgxdb_add_contams:
    input:
        db = "data/iptgxdbs/{iptgxdb}/{strain}/iptgxdb.fasta",
        contams = "data/contams.fasta"
    output:
        "results/iptgxdbs/{iptgxdb}/{strain}/iptgxdb_with_contams.fasta"
    shell:
        """
cat {input.db} > {output}
cat {input.contams} >> {output}
        """


rule iptgxdb_protparams:
    input:
        "data/iptgxdbs/{iptgxdb}/{strain}/iptgxdb.fasta"
    output:
        "results/iptgxdbs/{iptgxdb}/{strain}/iptgxdb_protparams.tsv"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/calc_prot_params.py -i '{input}' -o '{output}'
        """