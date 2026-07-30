rule protinf_collapse_extensions:
    input:
        "results/searches/{search}{subsearch}/protinf/protparam_filter/{strain}/proteins.selected.tsv"
    output:
        "results/searches/{search}{subsearch}/protinf/collapse_extensions/{strain}/proteins.collapsed.tsv"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/collapse_extensions.py \
-i '{input}' \
-o '{output}'
        """