rule proteinprophet:
    input:
        get_rescored_pepxmls_for_strain
    output:
        directory("results/searches/{search}{subsearch}/protinf/proteinprophet/{strain}")
    log:
        "results/searches/{search}{subsearch}/protinf/proteinprophet/{strain}/log.txt"
    params:
        input_files = lambda wildcards: [f"{workflow.basedir}/../{p}" for p in get_rescored_pepxmls_for_strain(wildcards)],
        settings = lambda wildcards: get_search_config_section(wildcards.search, "philosopher")
    shell:
        """
mkdir -p {output}
cd {output}
{params.settings[executable]} workspace --init --nocheck --analytics false 2>&1 | tee {workflow.basedir}/../{log}
{params.settings[executable]} proteinprophet {params.settings[proteinprophet_flags]} {params.input_files} 2>&1 | tee -a {workflow.basedir}/../{log}
        """