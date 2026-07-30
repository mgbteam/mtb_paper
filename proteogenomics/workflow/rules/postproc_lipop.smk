rule postproc_lipop:
    input:
        "results/searches/{search}/postproc/novels/extract_seqs/{strain}/novels.faa"
    output:
        folder = directory("results/searches/{search}/postproc/novels/lipop/{strain}"),
        gff = "results/searches/{search}/postproc/novels/lipop/{strain}/lipop.gff"
    params:
        executable = config["postproc"]["lipop"]["executable"]
    shell:
        """
{params.executable} '{input}' -workdir '{output.folder}' > '{output.gff}'
        """
