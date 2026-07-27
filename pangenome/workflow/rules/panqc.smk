rule panqc_asm_paths:
    input:
        lambda wildcards: expand("results/{{strains}}/convert/fna/{strain}.fna", strain=config["strains"][wildcards.strains])
    output:
        "results/{strains}/convert/fna/input_asm_paths.tsv"
    shell:
        """
        echo -e "SampleID\tGenome_ASM_PATH" > '{output}'
        
        for file in {input}; do
            strain="$(basename "$file" .fna)"
            echo -e "$strain\t$file" >> '{output}'
        done
        """


rule panqc_filter_gene_presence_absence:
    input:
        groups = "results/{strains}/{tool}/run/gene_presence_absence.csv",
        seqs = "results/{strains}/{tool}/run/pan_genome_reference.fa"
    output:
        "results/{strains}/{tool}/run/gene_presence_absence_filtered.csv"
    conda:
        "../envs/panaroo.yml"
    shell:
        """
        workflow/scripts/filter_gene_presence_absence.py -r '{input.seqs}' -m '{input.groups}' -o '{output}'
        """


rule panqc_run:
    input:
        asms = "results/{strains}/convert/fna/input_asm_paths.tsv",
        groups = "results/{strains}/{tool}/run/gene_presence_absence_filtered.csv",
        seqs = "results/{strains}/{tool}/run/pan_genome_reference.fa"
    output:
        folder = directory("results/{strains}/{tool}/panqc"),
        file = report("results/{strains}/{tool}/panqc/Step2_SeqClustering/NSC.ClusterInfo.tsv", category="{strains} {tool}", subcategory="PanQC")
    conda:
        "../envs/panqc.yml"
    shell:
        """
        panqc nrc -a '{input.asms}' -r '{input.seqs}' -m '{input.groups}' -o '{output.folder}'
        """


rule panqc_create_panaroo_alignments:
    input:
        new_groups = "results/{strains}/panaroo/panqc/Step2_SeqClustering/NSC.ClusterInfo.tsv",
        alignments = "results/{strains}/panaroo/run/aligned_gene_sequences"
    output:
        directory("results/{strains}/panaroo/panqc/aligned_new_groups")
    conda:
        "../envs/panaroo.yml"
    shell:
        """
        workflow/scripts/panqc_create_alignments.py -n '{input.new_groups}' -f '{input.alignments}' -o '{output}'
        """
