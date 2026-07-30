rule postproc_ampscanner2:
    input:
        "results/searches/{search}/postproc/novels/extract_seqs/{strain}/novels.faa"
    output:
        fasta = "results/searches/{search}/postproc/novels/ampscanner2/{strain}/ampscanner2.fasta",
        csv = "results/searches/{search}/postproc/novels/ampscanner2/{strain}/ampscanner2.csv"
    log:
        "results/searches/{search}/postproc/novels/ampscanner2/{strain}/ampscanner2.log"
    params:
        executable = config["postproc"]["ampscanner2"]["executable"],
        model = config["postproc"]["ampscanner2"]["model"]
    conda:
        "../envs/ampscanner2.yml"
    shell:
        """
{params.executable} -fasta '{input}' -model '{params.model}' -candidates '{output.fasta}' -preds '{output.csv}' 2>&1 | tee '{log}'
        """
