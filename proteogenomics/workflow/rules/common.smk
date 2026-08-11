wildcard_constraints:
    strain="[^/]+",
    search="[^/]+",
    subsearch=".{0}|/.+",
    subsample=".{0}|/.+",
    splitcatssfx=".{0}|[^/]+"


# Get samples for a strain
def get_strain_samples(strain):
    rawdir = config["rawfiles"]["folder"] + f"/{strain}"
    ext = config["rawfiles"]["extension"]
    samples = [os.path.splitext(f)[0] for f in os.listdir(rawdir) if f.endswith(ext)]
    return samples



# Query search settings
def get_search_config_section(search, path):
    final_settings_section = {}
    path_sections = path.split("/")
    settings_to_check = [
        config["search_default_settings"],
        config["searches"][search],
    ]

    for settings in settings_to_check:
        for path_section in path_sections:
            settings = settings.get(path_section, {})

        final_settings_section.update(settings)

    if final_settings_section:
        return final_settings_section

    print(f"Error, no local or global search settings section found for {path}")


def get_search_config_value(search, path):
    final_value = None
    path_sections = path.split("/")
    settings_to_check = [
        config["search_default_settings"],
        config["searches"][search],
    ]

    for settings in settings_to_check:
        for path_section in path_sections[:-1]:
            settings = settings.get(path_section, {})

        if path_sections[-1] in settings:
            final_value = settings[path_sections[-1]]

    if final_value is not None:
        return final_value

    print(f"Error, no local or global search settings value found for {path}")


def is_msbooster_enabled(search):
    return get_search_config_value(search, "msbooster/enable")


def is_entrapment_enabled(search):
    return get_search_config_value(search, "entrapment/enable")


def get_iptgxdb(search):
	iptgxdb = config["searches"][search].get("iptgxdb", config["search_default_settings"]["iptgxdb"])
	return iptgxdb
    

# Determine rule inputs
def get_search_folder_for_strain(wildcards):
    searchdir = f"results/searches/{wildcards.search}{wildcards.subsearch}/search"
    pepdir = "msbooster" if is_msbooster_enabled(wildcards.search) else "msfragger"
    straindir = f"{searchdir}/{pepdir}/{wildcards.strain}"
    return straindir


def get_rescored_pepxml_for_sample(wildcards):
    searchdir = f"results/searches/{wildcards.search}{wildcards.subsearch}/search"
    pepxml = f"{searchdir}/percolator/{wildcards.strain}/interact-{wildcards.sample}.pep.xml"
    return pepxml


def get_rescored_pepxmls_for_strain(wildcards):
    samples = get_strain_samples(wildcards.strain)
    searchdir = f"results/searches/{wildcards.search}{wildcards.subsearch}/search"
    pepxmls = [f"{searchdir}/percolator/{wildcards.strain}/interact-{s}.pep.xml" for s in samples]
    return sorted(pepxmls)


def get_filter_folders_for_strain(wildcards):
    filtdir = f"results/searches/{wildcards.search}{wildcards.subsearch}/protinf/filter_and_report/{wildcards.strain}"
    results = [f"{filtdir}/{sample}" for sample in get_strain_samples(wildcards.strain)]
    return sorted(results)