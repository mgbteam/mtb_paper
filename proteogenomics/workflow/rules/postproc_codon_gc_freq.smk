rule postproc_codon_gc_freq:
    input:
        novels = "results/searches/{search}/protinf/split_cats/{strain}/novels.tsv",
        genome = "data/annotations/refseq/{strain}.gb"
    output:
        table = "results/searches/{search}/postproc/novels/codon_gc_freq/{strain}.tsv",
        plot = report("results/searches/{search}/postproc/novels/codon_gc_freq/{strain}.png", category="{search}", subcategory="Novels Codon GC Freq.", labels={"Strain": "{strain}"})
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/codon_gc_freq.py -i '{input.novels}' -g '{input.genome}' -o '{output.table}' -p '{output.plot}'
        """
