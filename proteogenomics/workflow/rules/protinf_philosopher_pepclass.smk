rule protinf_pepclass_basic:
    input:
        dbdir = "results/searches/{search}{subsearch}/search/database/{strain}",
        folder = get_filter_folders_for_strain
    output:
        "results/searches/{search}{subsearch}/protinf/pepclass_filter/{strain}/pepclass.tsv"
    params:
        files = lambda wildcards: [f'{d}/peptide.tsv' for d in get_filter_folders_for_strain(wildcards)],
        contam_prefix = config["prefixes"]["contam"]
    conda:
        "../envs/python.yml"
    shell:
        """
db=$(find '{input.dbdir}'/ -name *.fas)

workflow/scripts/pepclass_basic.py \
-i "$db" \
-p {params.files} \
-o '{output}' \
-c '{params.contam_prefix}' \
-d 'rev_' \
-e '_p_target'
        """


rule protinf_pepclass_add:
    input:
        proteins = "results/searches/{search}{subsearch}/protinf/abacus/{strain}/proteins.reformat.tsv",
        pepclass = "results/searches/{search}{subsearch}/protinf/pepclass_filter/{strain}/pepclass.tsv"
    output:
        "results/searches/{search}{subsearch}/protinf/pepclass_filter/{strain}/protclass.tsv"
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/add_pepclass.py \
-t '{input.proteins}' \
-p '{input.pepclass}' \
-o '{output}'
        """


rule protinf_pepclass_filter:
    input:
        "results/searches/{search}{subsearch}/protinf/pepclass_filter/{strain}/protclass.tsv"
    output:
        proteins_uniq = "results/searches/{search}{subsearch}/protinf/pepclass_filter/{strain}/proteins.uniq.tsv",
        contams = "results/searches/{search}{subsearch}/protinf/pepclass_filter/{strain}/proteinsContams.tsv",
        prots3b = "results/searches/{search}{subsearch}/protinf/pepclass_filter/{strain}/proteins3b.tsv"
    params:
        contam_prefix = config["prefixes"]["contam"]
    conda:
        "../envs/python.yml"
    shell:
        """
workflow/scripts/pepclass_filter.py \
-i '{input}' \
-o '{output.proteins_uniq}' \
-a '{output.prots3b}' \
-c '{output.contams}' \
-p '{params.contam_prefix}'
        """