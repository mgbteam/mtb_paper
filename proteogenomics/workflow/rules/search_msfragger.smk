rule msfragger:
    input:
        dbdir = "results/searches/{search}{subsearch}/search/database/{strain}",
        rawdir = ancient("data/raw/{strain}")
    output:
        directory("results/searches/{search}{subsearch}/search/msfragger/{strain}"),
    log:
        "results/searches/{search}{subsearch}/search/msfragger/{strain}/log.txt"
    threads:
        lambda wildcards: get_search_config_section(wildcards.search, "msfragger")["threads"]
    params:
        settings = lambda wildcards: get_search_config_section(wildcards.search, "msfragger")
    conda:
        "../envs/openjdk.yml"
    shell:
        """
db=$(find '{input.dbdir}'/ -name *.fas)
config="{output}/$(basename '{params.settings[config]}')"

mkdir -p '{output}'
cp '{params.settings[config]}' "$config"
sed -i 's;^num_threads.*;num_threads = {threads};' "$config"
sed -i "s;^database_name.*;database_name = $db;" "$config"

cd '{output}'
ln -s '{workflow.basedir}/../{input.rawdir}'/*.d .
cd -

java -Xmx{params.settings[memory]} -jar '{params.settings[executable]}' "$config" '{output}'/*.d 2>&1 | tee {log}
        """