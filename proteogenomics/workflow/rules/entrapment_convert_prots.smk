rule entrapment_convert_prots:
    input:
        proteins = "results/searches/{search}/entrapment/protinf/filter_and_report/{strain}/{sample}/protein.tsv",
        iptgxdb = lambda wildcards: "data/iptgxdbs/" + get_iptgxdb(wildcards.search) + "/{strain}/iptgxdb.fasta",
        mapping = lambda wildcards: "data/iptgxdbs/" + get_iptgxdb(wildcards.search) + "/{strain}/iptgxdb.tsv"
    output:
        "results/searches/{search}/entrapment/protinf/convert_prots/{strain}/{sample}/protein.reformat.tsv"
    params:
        refseq_prefix = config["prefixes"]["refseq"],
        contam_prefix = config["prefixes"]["contam"]
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/reformat_report.py \
-i '{input.proteins}' \
-f '{input.iptgxdb}' \
-m '{input.mapping}' \
-o '{output}' \
-r '{params.refseq_prefix}' \
-c '{params.contam_prefix}' \
-d 'rev_' \
-e '_p_target' \
--filter-contams \
--filter-decoys
        """
