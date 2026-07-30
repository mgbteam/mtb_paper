rule protinf_abacus_run:
    input:
        protxmldir = "results/searches/{search}{subsearch}/protinf/proteinprophet/{strain}",
        dirs = get_filter_folders_for_strain
    output:
        "results/searches/{search}{subsearch}/protinf/abacus/{strain}/combined_protein.tsv"
    log:
        "results/searches/{search}{subsearch}/protinf/abacus/{strain}/log.txt"
    params:
        input_relative = lambda wildcards: [f"{workflow.basedir}/../{f}" for f in get_filter_folders_for_strain(wildcards)],
        settings = lambda wildcards: get_search_config_section(wildcards.search, "philosopher")
    shell:
        """
outdir="$(dirname '{output}')"
mkdir -p "$outdir"
cd "$outdir"
cp {workflow.basedir}/../{input.protxmldir}/interact.prot.xml combined.prot.xml
{params.settings[executable]} workspace --init --nocheck --analytics false 2>&1 | tee {workflow.basedir}/../{log}
{params.settings[executable]} abacus {params.settings[abacus_flags]} --tag rev_ --protein {params.input_relative} 2>&1 | tee -a {workflow.basedir}/../{log}
        """


rule protinf_abacus_reformat:
    input:
        proteins = "results/searches/{search}{subsearch}/protinf/abacus/{strain}/combined_protein.tsv",
        database = lambda wildcards: "data/iptgxdbs/" + get_iptgxdb(wildcards.search) + "/{strain}/iptgxdb.fasta",
        mapping = lambda wildcards: "data/iptgxdbs/" + get_iptgxdb(wildcards.search) + "/{strain}/iptgxdb.tsv"
    output:
        "results/searches/{search}{subsearch}/protinf/abacus/{strain}/proteins.reformat.tsv"
    params:
        refseq_prefix = config["prefixes"]["refseq"],
        contam_prefix = config["prefixes"]["contam"]
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/reformat_abacus.py \
-i '{input.proteins}' \
-f '{input.database}' \
-m '{input.mapping}' \
-o '{output}' \
-r '{params.refseq_prefix}' \
-c '{params.contam_prefix}'
        """