rule postproc_interpro:
    input:
        "results/searches/{search}/postproc/novels/extract_seqs/{strain}/novels.faa"
    output:
        folder = directory("results/searches/{search}/postproc/novels/interpro/{strain}"),
        tsv = "results/searches/{search}/postproc/novels/interpro/{strain}/novels.faa.tsv"
    threads:
        config["postproc"]["interpro"]["threads"]
    params:
        folder = config["postproc"]["interpro"]["folder"],
        flags = config["postproc"]["interpro"]["flags"]
    conda:
        "../envs/interproscan.yml"
    shell:
        """
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib"
mkdir -p {output.folder}

"{params.folder}/interproscan.sh" \
-cpu {threads} \
-i {input} \
-d {output.folder} \
-f TSV,JSON,GFF3 \
{params.flags}
        """
