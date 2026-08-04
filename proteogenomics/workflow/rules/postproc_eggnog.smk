rule postproc_eggnog:
    input:
        "results/searches/{search}/postproc/novels/extract_seqs/{strain}/novels.faa"
    output:
        folder = directory("results/searches/{search}/postproc/novels/eggnog/{strain}"),
        annots = "results/searches/{search}/postproc/novels/eggnog/{strain}/eggnog.emapper.annotations"
    threads:
        config["postproc"]["eggnog"]["threads"]
    params:
        folder = config["postproc"]["eggnog"]["folder"],
        flags = config["postproc"]["eggnog"]["flags"]
    conda:
        "../envs/eggnog.yml"
    shell:
        """
mkdir -p '{output.folder}'

emapper.py \
--cpu {threads} \
--data_dir '{params.folder}' \
--scratch_dir /tmp \
--temp_dir /tmp \
{params.flags} \
-i '{input}' \
-o eggnog \
--output_dir '{output.folder}' \
--override
        """
