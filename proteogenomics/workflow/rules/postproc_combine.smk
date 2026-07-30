# Analyses that can be run automatically
postproc_flags_and_files = {
    "eggnog":             ("-e", "results/searches/{search}/postproc/novels/eggnog/{strain}/eggnog.emapper.annotations"),
    "interpro":           ("-i", "results/searches/{search}/postproc/novels/interpro/{strain}/novels.faa.tsv"),
    "psortb":             ("-s", "results/searches/{search}/postproc/novels/psortb/{strain}/psortb.txt"),
    "lipop":              ("-l", "results/searches/{search}/postproc/novels/lipop/{strain}/lipop.gff"),
    "hhsuite":            ("-H", "results/searches/{search}/postproc/novels/hhsuite/{strain}/summary.tsv"),
    "ampscanner2":        ("-S", "results/searches/{search}/postproc/novels/ampscanner2/{strain}/ampscanner2.csv"),
    "conservation_blast": ("-c", "results/searches/{search}/postproc/novels/conservation_blast/{strain}/summary.tsv"),
    "codon_gc_freq":      ("-g", "results/searches/{search}/postproc/novels/codon_gc_freq/{strain}.tsv"),
}

# Analyses that must be run manually
postproc_manual_flags_and_files = {
    "phyre2":        ("-p", "results/searches/{search}/postproc/novels/phyre2/output/{strain}/summaryinfo.txt"),
    "operon_mapper": ("-m", "results/searches/{search}/postproc/operon_mapper/output/{strain}"),
}


def get_postproc_files(wildcards):
    folder = f"results/searches/{wildcards.search}"
    files = [f"{folder}/protinf/split_categories/{wildcards.strain}/novels.tsv"]

    # Annotation comparison results
    for annot in config["postproc"]["annot_overlap"]:
        filename = config["postproc"]["annot_overlap"][annot]["filename"]
        files.append(f"{folder}/postproc/novels/annot_overlap/{annot}/{wildcards.strain}.tsv")

    # Result files from all analyses that can be run automatically and are enabled
    for analysis, (flag, file) in postproc_flags_and_files.items():
        if config["postproc"][analysis]["enable"]:
            files.append(file.replace("{search}", wildcards.search).replace("{strain}", wildcards.strain))

    return files


def get_postproc_flags(wildcards):
    folder = f"results/searches/{wildcards.search}"
    flags = ["-t", f"{folder}/protinf/split_categories/{wildcards.strain}/novels.tsv"]

    # Annotation comparison results if any annotation given
    if config["postproc"]["annot_overlap"]:
        flags.append("-a")

    for annot in config["postproc"]["annot_overlap"]:
        filename = config["postproc"]["annot_overlap"][annot]["filename"]
        flags.append(f"{folder}/postproc/novels/annot_overlap/{annot}/{wildcards.strain}.tsv")

    # Result files from all analyses that can be run automatically and are enabled
    for analysis, (flag, file) in postproc_flags_and_files.items():
        if config["postproc"][analysis]["enable"]:
            flags.extend([flag, file.replace("{search}", wildcards.search).replace("{strain}", wildcards.strain)])

    # Result files from all analyses that must be run manually, are enabled and results are present
    for analysis, (flag, file) in postproc_manual_flags_and_files.items():
        if config["postproc"][analysis]["enable"]:
            file = file.replace("{search}", wildcards.search).replace("{strain}", wildcards.strain)

            if os.path.exists(file):
                flags.extend([flag, file])

    return flags


rule postproc_combine:
    input:
        get_postproc_files
    output:
        "results/searches/{search}/postproc/novels/combined/{strain}.tsv"
    params:
        flags = get_postproc_flags
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/combine_postproc.py \
{params.flags} \
-o '{output}'
        """
