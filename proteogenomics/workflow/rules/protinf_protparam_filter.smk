rule protinf_protparam_add:
    input:
        proteins = "results/searches/{search}{subsearch}/protinf/pepclass_filter/{strain}/proteins.uniq.tsv",
        protparams = lambda wildcards: "results/iptgxdbs/" + get_iptgxdb(wildcards.search) + "/{strain}/iptgxdb_protparams.tsv"
    output:
        "results/searches/{search}{subsearch}/protinf/protparam_filter/{strain}/proteins.uniq.tsv"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/combine_tables.py \
-l 'protein' \
-r 'Protein' \
-m 'left' \
-o '{output}' \
'{input.proteins}' \
'{input.protparams}'
        """


def create_protparam_threshold_flags(search):
    filter_settings = get_search_config_section(search, "filter")

    flags = [f"--min-psms '{cat}' {val}" for cat, val in filter_settings["psms"].items()]
    flags.append(f"--min-peps " + str(filter_settings["peps"]["min_peps"]))
    flags.append(f"--min-weight " + str(filter_settings["peps"]["min_weight"]))
    return " ".join(flags)


rule protinf_protparam_filter:
    input:
        proteins = "results/searches/{search}{subsearch}/protinf/protparam_filter/{strain}/proteins.uniq.tsv"
    output:
        selected = "results/searches/{search}{subsearch}/protinf/protparam_filter/{strain}/proteins.selected.tsv",
        discarded = "results/searches/{search}{subsearch}/protinf/protparam_filter/{strain}/proteins.discarded.tsv"
    log:
        "results/searches/{search}{subsearch}/protinf/protparam_filter/{strain}/log.txt"
    params:
        threshold_flags = lambda wildcards: create_protparam_threshold_flags(wildcards.search)
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/protparam_filter.py \
{params.threshold_flags} \
-p "{input.proteins}" \
-o "{output.selected}" \
-O "{output.discarded}" \
> "{log}"
        """