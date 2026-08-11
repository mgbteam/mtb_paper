def create_overview_flags(wildcards):
    flags = []
    protinf_dir = f"results/searches/{wildcards.search}{wildcards.subsearch}/protinf"

    for strain in config["strains"]:
        summary_file = f"split_cats/{strain}/summary.tsv"
        flags.append(f"-i '{strain}' '{protinf_dir}/{summary_file}'")

    return " ".join(flags)


rule protinf_summary:
    input:
        expand("results/searches/{{search}}{{subsearch}}/protinf/split_cats/{strain}/summary.tsv", strain=config["strains"])
    output:
        "results/searches/{search}{subsearch}/protinf/protinf_summary.tsv"
    params:
        flags = create_overview_flags
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/protinf_summary.py {params.flags} -o "{output}"
        """
