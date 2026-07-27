def create_color_flags(wildcards):
    for color, strains in config["upsetplots"]["colors"].items():
        yield f"-c '{color}' " + " ".join([f"'{strain}'" for strain in strains])
