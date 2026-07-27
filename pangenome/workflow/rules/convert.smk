rule convert_to_fna:
    input:
        "data/annotation/{strain}.gbff"
    output:
        "results/{strains}/convert/fna/{strain}.fna"
    conda:
        "../envs/panaroo.yml"
    shell:
        """
        workflow/scripts/gbk_to_fna.py -i '{input}' -o '{output}'
        """


rule convert_to_gff:
    input:
        "data/annotation/{strain}.gbff"
    output:
        "results/{strains}/convert/gff/{strain}.gff"
    params:
        pseudo_flag = "" if config["convert_pseudo"] else "--no-pseudo"
    conda:
        "../envs/panaroo.yml"
    shell:
        """
        workflow/scripts/gbk_to_gff.py -i '{input}' -o '{output}' {params.pseudo_flag}
        """
