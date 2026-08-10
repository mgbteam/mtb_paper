rule protinf_collapse_extensions:
    input:
        "results/searches/{search}{subsearch}/protinf/protparam_filter/{strain}{subsample}/proteins.selected.tsv"
    output:
        "results/searches/{search}{subsearch}/protinf/collapse_extensions/{strain}{subsample}/proteins.collapsed.tsv"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/collapse_extensions.py \
-i '{input}' \
-o '{output}' \
-e '_p_target'
        """