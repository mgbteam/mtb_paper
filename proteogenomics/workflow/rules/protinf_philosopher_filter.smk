rule protinf_philosopher_filter:
    input:
        dbdir = "results/searches/{search}{subsearch}/search/database/{strain}",
        pepxmldir = get_rescored_pepxml_for_sample,
        protxmldir = "results/searches/{search}{subsearch}/protinf/proteinprophet/{strain}"
    output:
        directory("results/searches/{search}{subsearch}/protinf/filter_and_report/{strain}/{sample}")
    log:
        "results/searches/{search}{subsearch}/protinf/filter_and_report/{strain}/{sample}/log.txt"
    params:
        pepxml = lambda wildcards: f"{workflow.basedir}/../" + get_rescored_pepxml_for_sample(wildcards),
        settings = lambda wildcards: get_search_config_section(wildcards.search, "philosopher")
    shell:
        """
mkdir -p {output}
cd {output}
{params.settings[executable]} workspace --init --nocheck --analytics false 2>&1 | tee {workflow.basedir}/../{log}
{params.settings[executable]} database --annotate {workflow.basedir}/../{input.dbdir}/*.fas --prefix rev_
{params.settings[executable]} filter {params.settings[filter_flags]} --tag rev_ --pepxml {params.pepxml} --protxml {workflow.basedir}/../{input.protxmldir}/interact.prot.xml 2>&1 | tee {workflow.basedir}/../{log}
{params.settings[executable]} report {params.settings[report_flags]}
        """
