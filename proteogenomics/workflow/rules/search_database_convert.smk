rule create_mainsearch_prot_db:
    input:
        lambda wildcards: "results/iptgxdbs/" + get_iptgxdb(wildcards.search) + "/{strain}/iptgxdb_with_contams.fasta"
    output:
        directory("results/searches/{search}/search/database/{strain}")
    log:
        "results/searches/{search}/search/database/{strain}/log.txt"
    params:
        settings = lambda wildcards: get_search_config_section(wildcards.search, "philosopher")
    shell:
        """
cd {output}
{params.settings[executable]} workspace --init --nocheck --analytics false 2>&1 | tee {workflow.basedir}/../{log}
{params.settings[executable]} database --custom {workflow.basedir}/../{input} 2>&1 | tee -a {workflow.basedir}/../{log}
        """