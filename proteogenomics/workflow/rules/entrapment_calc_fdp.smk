rule entrapment_calc_prot_fdp_percolator:
    input:
        "results/searches/{search}/entrapment/fdrbench/convert/{strain}/{sample}/{subset}.tsv"
    output:
        "results/searches/{search}/entrapment/fdrbench/calc_fdp/{strain}/{sample}/{subset}_fdp.csv"
    params:
        settings = lambda wildcards: get_search_config_section(wildcards.search, "entrapment")
    conda:
        "../envs/openjdk.yml"
    shell:
        """
java -jar {params.settings[jar]} \
{params.settings[fdp_flags]} \
-level protein \
-i '{input}' \
-o '{output}' \
-score 'score:0' \
-fold 1 \
-pick first
        """