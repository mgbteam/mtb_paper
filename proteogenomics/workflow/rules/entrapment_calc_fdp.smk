def get_entrapment_input(wildcards):
    basedir = f"results/searches/{wildcards.search}/entrapment/protinf"

    if wildcards.stage == "pre_protparam":
        folder = "split_cats_pre_protparam"
    else:
        folder = "split_cats"

    return f"{basedir}/{folder}/{wildcards.strain}/{wildcards.sample}/{wildcards.subset}.tsv"


rule entrapment_calc_prot_fdp:
    input:
        get_entrapment_input
    output:
        "results/searches/{search}/entrapment/fdrbench/{stage}/calc_fdp/{strain}/{sample}/{subset}_fdp.csv"
    params:
        settings = lambda wildcards: get_search_config_section(wildcards.search, "entrapment")
    conda:
        "../envs/openjdk.yml"
    shell:
        """
if [ "$(cat '{input}' | wc -l)" -gt "1" ]; then
    java -jar {params.settings[jar]} \
    {params.settings[fdp_flags]} \
    -level protein \
    -i '{input}' \
    -o '{output}' \
    -score 'score:0' \
    -fold 1 \
    -pick first
else
    touch '{output}'
fi
        """