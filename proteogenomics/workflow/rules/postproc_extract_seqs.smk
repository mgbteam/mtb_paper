rule postproc_extract_combined_seqs:
    input:
        novels = "results/searches/{search}/protinf/split_cats/{strain}/novels.tsv",
        genome = "data/genomes/{strain}.fasta"
    output:
        "results/searches/{search}/postproc/novels/extract_seqs/{strain}/novels.faa"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/extract_seqs.py \
-i '{input.novels}' \
-g '{input.genome}' \
-o '{output}'
        """


checkpoint postproc_extract_individual_seqs:
    input:
        "results/searches/{search}/postproc/novels/extract_seqs/{strain}/novels.faa"
    output:
        directory("results/searches/{search}/postproc/novels/extract_seqs/{strain}/individual")
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/split_fasta.py -i '{input}' -o '{output}'
        """
