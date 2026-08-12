rule postproc_interpro:
    input:
        "results/searches/{search}/postproc/novels/extract_seqs/{strain}/novels.faa"
    output:
        tmpfolder = temp(directory(config["postproc"]["interpro"]["tempfolder"] + "/interpro/{search}/{strain}")),
        outfolder = directory("results/searches/{search}/postproc/novels/interpro/{strain}"),
        tsv = "results/searches/{search}/postproc/novels/interpro/{strain}/novels.faa.tsv"
    threads:
        config["postproc"]["interpro"]["threads"]
    params:
        datafolder = config["postproc"]["interpro"]["datafolder"],
        flags = config["postproc"]["interpro"]["flags"]
    conda:
        "../envs/interproscan.yml"
    shell:
        """
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib"
mkdir -p '{output.tmpfolder}' '{output.outfolder}'

if [ -s '{input}' ]; then
    "{params.datafolder}/interproscan.sh" \
        -cpu {threads} \
        -i '{input}' \
        -d '{output.outfolder}' \
        -f TSV,JSON,GFF3 \
        -T '{output.tmpfolder}' \
        {params.flags}
else
    touch '{output.tsv}'
fi
        """
