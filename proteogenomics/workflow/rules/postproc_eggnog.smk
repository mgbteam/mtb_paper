rule postproc_eggnog:
    input:
        "results/searches/{search}/postproc/novels/extract_seqs/{strain}/novels.faa"
    output:
        scratchfolder = temp(directory(config["postproc"]["eggnog"]["tempfolder"] + "/eggnog_scratch/{search}/{strain}")),
        tempfolder = temp(directory(config["postproc"]["eggnog"]["tempfolder"] + "/eggnog_temp/{search}/{strain}")),
        outfolder = directory("results/searches/{search}/postproc/novels/eggnog/{strain}"),
        annots = "results/searches/{search}/postproc/novels/eggnog/{strain}/eggnog.emapper.annotations"
    threads:
        config["postproc"]["eggnog"]["threads"]
    params:
        datafolder = config["postproc"]["eggnog"]["datafolder"],
        flags = config["postproc"]["eggnog"]["flags"]
    conda:
        "../envs/eggnog.yml"
    shell:
        """
mkdir -p '{output.scratchfolder}' '{output.tempfolder}' '{output.outfolder}'

if [ -s '{input}' ]; then
    emapper.py \
        --cpu {threads} \
        --data_dir '{params.datafolder}' \
        --scratch_dir '{output.scratchfolder}' \
        --temp_dir '{output.tempfolder}' \
        {params.flags} \
        -i '{input}' \
        -o eggnog \
        --output_dir '{output.outfolder}' \
        --override
else
    touch '{output.annots}'
fi
        """
