rule entrapment_convert_prots:
    input:
        "results/searches/{search}/entrapment/search/percolator/{strain}/{sample}_target_prots.tsv"
    output:
        all = "results/searches/{search}/entrapment/fdrbench/convert/{strain}/{sample}/all.tsv",
        refseq = "results/searches/{search}/entrapment/fdrbench/convert/{strain}/{sample}/refseq.tsv",
        novels = "results/searches/{search}/entrapment/fdrbench/convert/{strain}/{sample}/novels.tsv",
        contams = "results/searches/{search}/entrapment/fdrbench/convert/{strain}/{sample}/contams.tsv"
    params:
        refseq_prefix = config["prefixes"]["refseq"],
        contam_prefix = config["prefixes"]["contam"]
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/percolator_convert_prots.py \
'{input}' \
"$(dirname '{output.all}')" \
--refseq-prefix '{params.refseq_prefix}' \
--contam-prefix '{params.contam_prefix}'
        """
