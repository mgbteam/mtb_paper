rule postproc_phyre2_create_input:
    input:
        "results/searches/{search}/postproc/novels/extract_seqs/{strain}/novels.faa"
    output:
        "results/searches/{search}/postproc/novels/phyre2/input/{strain}.faa"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/filter_seqs.py \
-s 30 \
-i {input} \
-o {output}
        """
