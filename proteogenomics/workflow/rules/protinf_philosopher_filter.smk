def get_philosopher_filter_flags(wildcards):
    if wildcards.subsearch.startswith("/entrapment"):
        key = "filter_flags_entrapment"
    else:
        key = "filter_flags"

    return get_search_config_value(wildcards.search, f"philosopher/{key}")


rule protinf_philosopher_filter:
    input:
        dbdir = "results/searches/{search}{subsearch}/search/database/{strain}",
        pepxmldir = get_rescored_pepxml_for_sample,
        protxmldir = "results/searches/{search}{subsearch}/protinf/proteinprophet/{strain}"
    output:
        folder = directory("results/searches/{search}{subsearch}/protinf/filter_and_report/{strain}/{sample}"),
        proteins = "results/searches/{search}{subsearch}/protinf/filter_and_report/{strain}/{sample}/protein.tsv"
    log:
        "results/searches/{search}{subsearch}/protinf/filter_and_report/{strain}/{sample}/log.txt"
    params:
        pepxml = lambda wildcards: f"{workflow.basedir}/../" + get_rescored_pepxml_for_sample(wildcards),
        executable = lambda wildcards: get_search_config_value(wildcards.search, "philosopher/executable"),
        filter_flags = get_philosopher_filter_flags,
        report_flags = lambda wildcards: get_search_config_value(wildcards.search, "philosopher/report_flags")
    shell:
        """
mkdir -p {output.folder}
cd {output.folder}
{params.executable} workspace --init --nocheck --analytics false 2>&1 | tee {workflow.basedir}/../{log}
{params.executable} database --annotate {workflow.basedir}/../{input.dbdir}/*.fas --prefix rev_
{params.executable} filter {params.filter_flags} --tag rev_ --pepxml {params.pepxml} --protxml {workflow.basedir}/../{input.protxmldir}/interact.prot.xml 2>&1 | tee {workflow.basedir}/../{log}
{params.executable} report {params.report_flags}
        """
