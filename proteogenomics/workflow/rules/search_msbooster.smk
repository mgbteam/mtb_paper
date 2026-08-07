rule msbooster:
    input:
        "results/searches/{search}{subsearch}/search/msfragger/{strain}"
    output:
        directory("results/searches/{search}{subsearch}/search/msbooster/{strain}")
    threads:
        lambda wildcards: get_search_config_value(wildcards.search, "msbooster/threads")
    params:
        msfragger_config = lambda wildcards: get_search_config_value(wildcards.search, "msfragger/config"),
        settings = lambda wildcards: get_search_config_section(wildcards.search, "msbooster")
    conda:
        "../envs/openjdk.yml"
    shell:
        """
mkdir -p "{output}"
msbooster_config="$(basename '{params.settings[config]}')"
msfragger_config="$(basename '{params.msfragger_config}')"

cp "{params.settings[config]}" "{output}/$msbooster_config"
sed -i "s;^numThreads.*;numThreads = {threads};" "{output}/$msbooster_config"
sed -i "s;^fragger.*;fragger = {input}/$msfragger_config;" "{output}/$msbooster_config"
sed -i "s;^pinPepXMLDirectory.*;pinPepXMLDirectory = {input};" "{output}/$msbooster_config"
sed -i "s;^mzmlDirectory.*;mzmlDirectory = {input};" "{output}/$msbooster_config"

java -Xmx{params.settings[memory]} \
-cp '{params.settings[msbooster_classpath]}:{params.settings[batmass_classpath]}' \
Features.MainClass \
--paramsList "{output}/$msbooster_config"

shopt -s nullglob
for file in "{input}"/MSBooster*; do mv "$file" "{output}"/; done
for file in "{input}"/spectraRT*; do mv "$file" "{output}"/; done
for file in "{input}"/*_edited.pin; do mv "$file" "{output}/$(basename "$file" | sed 's/_edited//')"; done
        """