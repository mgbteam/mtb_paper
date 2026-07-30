rule percolator:
    input:
        dbdir = "results/searches/{search}{subsearch}/search/database/{strain}",
        folder = get_search_folder_for_strain
    output:
        target_psms = "results/searches/{search}{subsearch}/search/percolator/{strain}/{sample}_target_psms.tsv",
        target_peps = "results/searches/{search}{subsearch}/search/percolator/{strain}/{sample}_target_peps.tsv",
        target_prots = "results/searches/{search}{subsearch}/search/percolator/{strain}/{sample}_target_prots.tsv",
        decoy_psms = "results/searches/{search}{subsearch}/search/percolator/{strain}/{sample}_decoy_psms.tsv",
        decoy_peps = "results/searches/{search}{subsearch}/search/percolator/{strain}/{sample}_decoy_peps.tsv",
        decoy_prots = "results/searches/{search}{subsearch}/search/percolator/{strain}/{sample}_decoy_prots.tsv"
    threads:
        lambda wildcards: get_search_config_value(wildcards.search, "percolator/threads")
    params:
        settings = lambda wildcards: get_search_config_section(wildcards.search, "percolator")
    shell:
        """
db=$(find '{input.dbdir}'/ -name *.fas)

{params.settings[executable]} {params.settings[flags]} \
--num-threads {threads} \
--results-psms '{output.target_psms}' \
--decoy-results-psms '{output.decoy_psms}' \
--results-peptides '{output.target_peps}' \
--decoy-results-peptides '{output.decoy_peps}' \
--picked-protein "$db" \
--protein-decoy-pattern 'rev_' \
--results-proteins '{output.target_prots}' \
--decoy-results-proteins '{output.decoy_prots}' \
'{input.folder}/{wildcards.sample}.pin'
        """


rule percolator_convert:
    input:
        pindir = get_search_folder_for_strain,
        fragdir = "results/searches/{search}{subsearch}/search/msfragger/{strain}",
        target = "results/searches/{search}{subsearch}/search/percolator/{strain}/{sample}_target_psms.tsv",
        decoy = "results/searches/{search}{subsearch}/search/percolator/{strain}/{sample}_decoy_psms.tsv"
    output:
        "results/searches/{search}{subsearch}/search/percolator/{strain}/interact-{sample}.pep.xml"
    params:
        settings = lambda wildcards: get_search_config_section(wildcards.search, "percolator")
    conda:
        "../envs/openjdk.yml"
    shell:
        """
java -cp "{params.settings[fragpipe_classpath]}/*" \
com.dmtavt.fragpipe.tools.percolator.PercolatorOutputToPepXML \
"{input.pindir}/{wildcards.sample}.pin" \
"{input.fragdir}/{wildcards.sample}" \
"{input.target}" \
"{input.decoy}" \
"$(echo '{output}' | sed 's/.pep.xml$//')" \
DDA \
0 \
""
        """