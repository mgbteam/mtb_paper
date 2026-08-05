rule postproc_hhsuite:
    input:
        "results/searches/{search}/postproc/novels/extract_seqs/{strain}/individual/{novel}.faa"
    output:
        "results/searches/{search}/postproc/novels/hhsuite/{strain}/results/{novel}.txt"
    threads:
        config["postproc"]["hhsuite"]["threads"]
    params:
        database = config["postproc"]["hhsuite"]["database"],
        flags = config["postproc"]["hhsuite"]["flags"]
    conda:
        "../envs/hhsuite.yml"
    shell:
        """
hhblits \
-i '{input}' \
-o '{output}' \
-d '{params.database}' \
{params.flags} \
-cpu {threads} \
        """


def get_hhsuite_results_per_novel(wildcards):
    folder = f"results/searches/{wildcards.search}/postproc/novels/hhsuite/{wildcards.strain}/results"
    checkpoint_output = checkpoints.postproc_extract_individual_seqs.get(**wildcards).output[0]
    novel, = glob_wildcards(os.path.join(checkpoint_output, "{novel}.faa"))
    return expand(os.path.join(folder, "{novel}.txt"), novel=novel)


rule postproc_hhsuite_summary:
    input:
        get_hhsuite_results_per_novel
    output:
        "results/searches/{search}/postproc/novels/hhsuite/{strain}/summary.tsv"
    params:
        quoted_input = lambda wildcards, input: [f"'{f}'" for f in input],
        exclude = config["postproc"]["hhsuite"]["exclude"]
    conda:
        "../envs/python.yml"
    shell:
        """
if [ -n '{input}' ]; then
    workflow/scripts/hhsuite_summary.py \
        -i {params.quoted_input} \
        -e {params.exclude} \
        -o '{output}'
else
    touch '{output}'
fi
        """
