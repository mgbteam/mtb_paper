def get_fdp_filter_input(wildcards):
    basedir = f"results/searches/{wildcards.search}/entrapment"
    protdir = f"{basedir}/protinf/filter_and_report/{wildcards.strain}"
    fdpdir = f"{basedir}/fdrbench/pre_protparam/calc_fdp/{wildcards.strain}"
    
    for sample in get_strain_samples(wildcards.strain):
        yield f"{protdir}/{sample}/protein.tsv"
        yield f"{fdpdir}/{sample}/novels_fdp.csv"


def get_fdp_filter_input_flags(wildcards):
    basedir = f"results/searches/{wildcards.search}/entrapment"
    protdir = f"{basedir}/protinf/filter_and_report/{wildcards.strain}"
    fdpdir = f"{basedir}/fdrbench/pre_protparam/calc_fdp/{wildcards.strain}"
    
    for sample in get_strain_samples(wildcards.strain):
        yield f"-i '{protdir}/{sample}/protein.tsv' '{fdpdir}/{sample}/novels_fdp.csv'"


rule postproc_fdp_filter:
    input:
        fdp_files = get_fdp_filter_input,
        proteins = "results/searches/{search}/protinf/split_cats/{strain}/novels.tsv"
    output:
        "results/searches/{search}/postproc/novels/fdp_filter/{strain}.tsv"
    params:
        flags = lambda wildcards: get_fdp_filter_input_flags(wildcards)
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/fdp_filter.py \
{params.flags} \
-p '{input.proteins}' \
-o '{output}'
        """