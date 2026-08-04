rule postproc_psortb:
    input:
        "results/searches/{search}/postproc/novels/extract_seqs/{strain}/novels.faa"
    output:
        folder = directory("results/searches/{search}/postproc/novels/psortb/{strain}"),
        txt = "results/searches/{search}/postproc/novels/psortb/{strain}/psortb.txt"
    params:
        flags = config["postproc"]["psortb"]["flags"]
    shell:
        """
mkdir -p '{output.folder}'
workflow/scripts/psortb {params.flags} -i '{input}' -r '{output.folder}'
mv '{output.folder}'/*.txt '{output.txt}'
        """
