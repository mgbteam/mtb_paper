rule postproc_conservation_blast:
    input:
        "results/searches/{search}/postproc/novels/extract_seqs/{strain}/individual/{novel}.faa"
    output:
        "results/searches/{search}/postproc/novels/conservation_blast/{strain}/results/{novel}.tsv"
    threads:
        config["postproc"]["conservation_blast"]["threads"]
    params:
        blastdb = config["postproc"]["conservation_blast"]["blastdb"],
        flags = config["postproc"]["conservation_blast"]["flags"],
        columns = " ".join(config["postproc"]["conservation_blast"]["columns"])
    conda:
        "../envs/blast.yml"
    shell:
        """
echo "$(echo '{params.columns}' | sed 's/ /\t/g')" > '{output}'

tblastn \
-num_threads {threads} \
{params.flags} \
-outfmt '6 {params.columns}' \
-query '{input}' \
-db '{params.blastdb}' \
>> '{output}'
        """


def get_blast_results_per_novel(wildcards):
    folder = f"results/searches/{wildcards.search}/postproc/novels/conservation_blast/{wildcards.strain}/results"
    checkpoint_output = checkpoints.postproc_extract_individual_seqs.get(**wildcards).output[0]
    novel, = glob_wildcards(os.path.join(checkpoint_output, "{novel}.faa"))
    return expand(os.path.join(folder, "{novel}.tsv"), novel=novel)


rule postproc_conservation_blast_summary:
    input:
        get_blast_results_per_novel
    output:
        summary = "results/searches/{search}/postproc/novels/conservation_blast/{strain}/summary.tsv",
        subtaxa_counts = "results/searches/{search}/postproc/novels/conservation_blast/{strain}/subtaxa_counts.tsv"
    params:
        quoted_input = lambda wildcards, input: [f"'{f}'" for f in input],
        tax_levels = [f"'{t}'" for t in config["postproc"]["conservation_blast"]["tax_levels"]],
        min_coverage = config["postproc"]["conservation_blast"]["min_coverage"],
        min_identity = config["postproc"]["conservation_blast"]["min_identity"],
        max_evalue = config["postproc"]["conservation_blast"]["max_evalue"],
        taxdbflag = "-d " + config["postproc"]["conservation_blast"]["taxdb"] if config["postproc"]["conservation_blast"]["taxdb"] else "",
    conda:
        "../envs/ete3.yml"
    shell:
        """
workflow/scripts/conservation_blast_summary.py \
-b {params.quoted_input} \
-t {params.tax_levels} \
-c {params.min_coverage} \
-i {params.min_identity} \
-e {params.max_evalue} \
-o {output.summary} \
-s {output.subtaxa_counts} \
{params.taxdbflag}
        """
