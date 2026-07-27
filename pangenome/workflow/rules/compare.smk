rule compare_plots:
    input:
        roary = "results/{strains}/roary/run/gene_presence_absence.csv",
        panaroo = "results/{strains}/panaroo/run/gene_presence_absence.csv",
    output:
        venn = report("results/{strains}/compare/plots/roary_vs_panaroo_venn.png", category="{strains} comparison", subcategory="Plots"),
        boxplot = report("results/{strains}/compare/plots/roary_vs_panaroo_boxplot.png", category="{strains} comparison", subcategory="Plots")
    conda:
        "../envs/venn.yml"
    shell:
        """
        workflow/scripts/plot_pangenome_differences.py -r '{input.roary}' -p '{input.panaroo}' -v '{output.venn}' -b '{output.boxplot}'
        """


rule compare_groups:
    input:
        roary = "results/{strains}/roary/run/gene_presence_absence.csv",
        panaroo = "results/{strains}/panaroo/run/gene_presence_absence.csv",
    output:
        common = report("results/{strains}/compare/groups/common.tsv", category="{strains} comparison", subcategory="Groups"),
        unique_roary = report("results/{strains}/compare/groups/unique_roary.tsv", category="{strains} comparison", subcategory="Groups"),
        unique_panaroo = report("results/{strains}/compare/groups/unique_panaroo.tsv", category="{strains} comparison", subcategory="Groups"),
        unique_matched = report("results/{strains}/compare/groups/unique_matched.tsv", category="{strains} comparison", subcategory="Groups")
    conda:
        "../envs/venn.yml"
    shell:
        """
        workflow/scripts/compare_pangenome_groups.py \
            -r '{input.roary}' \
            -p '{input.panaroo}' \
            -C '{output.common}' \
            -R '{output.unique_roary}' \
            -P '{output.unique_panaroo}' \
            -U '{output.unique_matched}'
        """


rule compare_alignment_entropy:
    input:
        entropy = "results/{strains}/panaroo/run/alignment_entropy.csv",
        unique_panaroo = "results/{strains}/compare/groups/unique_panaroo.tsv"
    output:
        report("results/{strains}/compare/entropy/alignment_entropy.png", category="{strains} comparison", subcategory="Entropy")
    conda:
        "../envs/panaroo.yml"
    shell:
        """
        workflow/scripts/compare_alignment_entropy.py \
            -e '{input.entropy}' \
            -u '{input.unique_panaroo}' \
            -o '{output}'
        """
