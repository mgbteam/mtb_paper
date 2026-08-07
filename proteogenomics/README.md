# *M. tuberculosis* proteogenomics Snakemake pipeline
This Snakemake pipeline performs a proteogenomic analysis of six *M. tuberculosis* clinical reference strains from lineages 1 and 2 ([Heiniger et al., 2026](https://doi.org/10.64898/2026.01.27.701740)). The workflow combines genome annotation processing, peptide-spectrum matching, protein inference, stringent filtering and downstream annotation and structural analyses to identify and prioritize novel proteins.

The pipeline converts multiple annotation sources to create 2 integrated proteogenomics databases ([iPtgxDBs](https://iptgxdb.expasy.org/)) for each strain: a standard iPtgxDB ([Omasits et al. 2017](https://doi.org/10.1101/gr.218255.116)) and a much smaller custom iPtgxDB based on Ribo-seq data ([Hadjeras et al. 2023](https://doi.org/10.1093/femsml/uqad012)), here using orthologs of Ribo-seq identifications from strain H37Rv ([Smith et al. 2022](https://doi.org/10.7554/eLife.73980)). Mass spectrometry data from a Bruker timsTOF Pro device ([PRIDE Project PXD081163](https://www.ebi.ac.uk/pride/archive/projects/PXD081163)) is then searched for each strain against the respective iPtgxDBs using [MSFragger](https://msfragger.nesvilab.org/), relying on [MSBooster](https://doi.org/10.1038/s41467-023-40129-9) and [Percolator](https://github.com/percolator/percolator) for accurate peptide identification. Protein inference and FDR filtering is achieved using [Philosopher](https://philosopher.nesvilab.org/), followed by additional stringent filtering as described ([Heiniger et al., 2026](https://doi.org/10.64898/2026.01.27.701740)), including the selection of unambiguous peptides based on the PeptideClassifier software ([Qeli and Anhres, 2010](https://doi.org/10.1038/nbt0710-647)) extended for prokaryotic proteogenomics ([Omasits et al, 2017](https://doi.org/10.1101/gr.218255.116)).

Novel genes not annotated before, corrected start sites of annotated genes and pseudogenes with confirmed expression are reported separately. Novel gene candidates are further analyzed, including the overlap with pseudogenes and newer annotations, signal peptide and subcellular localization prediction, functional annotation, and integration with Panaroo-based pangenome information.

[TOC]: #
## Table of Contents
- [Software Dependencies](#software-dependencies)
  - [Automatically Installed](#automatically-installed)
  - [Manually Installed](#manually-installed)
    - [Mandatory](#mandatory)
    - [Optional](#optional)
- [Input Data and Configuration](#input-data-and-configuration)
  - [Provided Data to Reproduce Results](#provided-data-to-reproduce-results)
  - [Configuration](#configuration)
- [Running the Pipeline](#running-the-pipeline)
  - [Optional Analyses Requiring Manual Intervention](#optional-analyses-requiring-manual-intervention)
- [Output](#output)
- [Citation](#citation)

## Software Dependencies
### Automatically Installed
By default, dependencies are installed automatically using conda when running the pipeline for the first time, except for FragPipe and dependencies used for optional downstream analyses of the identified novel proteins. The path of each manually installed software has to be configured in `config/config.yaml`. To disable automatic dependency management, remove `--use-conda` from `run.sh` and make sure the required tools are installed manually.

If PSORTb is enabled, make sure to [install Docker](https://docs.docker.com/engine/install) and also enable [management of Docker as a non-root user](https://docs.docker.com/engine/install/linux-postinstall).

### Manually Installed
#### Mandatory
FragPipe, MSFragger, Philosopher and have to be installed manually. While MSFragger can be installed from within the FragPipe GUI, Philosopher should be installed manually as a newer version was used here to calculate protein-level q-values. The following versions were used for the pipeline:

|Software|Version|
|--------|------:|
|[FragPipe](https://fragpipe.nesvilab.org)|22.0|
|[MSFragger](https://github.com/Nesvilab/MSFragger)|4.1|
|[Philosopher](https://github.com/Nesvilab/philosopher)|5.1.2|

All other dependencies of FragPipe are either automatically installed with conda or already provided with FragPipe. These include:

|Software|Version|
|--------|------:|
|[Percolator](https://github.com/percolator/percolator)|3.6.5|
|[MSBooster](https://github.com/Nesvilab/MSBooster)|1.2.31|
|[batmass-io](https://batmass.org)|1.33.4|

After installation, the paths to the 6 tools mentioned above have to be adapted in `config/config.yaml`.

#### Optional
These dependencies have to be installed and the path adjusted in `config/config.yaml` or the analysis has to be disabled by setting `enable: False` in the corresponding config section.

Either the software itself has be installed manually:

|Software|Version|Size|Notes|
|--------|------:|---:|-----|
|[FDRBench](https://github.com/Noble-Lab/FDRBench)|0.0.4|74 Mb|Estimation of actual subset FDR|
|[amp-scanner](https://github.com/dan-veltri/amp-scanner-v2)|2|28 Mb|Antimicrobial peptide prediction|
|[InterProScan](https://github.com/ebi-pf-team/interproscan)|5.59-91.0|47 Gb|Identification of functional domains
|[LipoP](https://services.healthtech.dtu.dk/services/LipoP-1.0/)|1.0a|632 Kb|Lipoprotein and signal peptide detection|

Or the software is automatically installed with conda but a dataset has to be downloaded manually:

|Dataset|Size|Software|Notes|
|-------|---:|-------:|-----|
|[core_nt](https://ftp.ncbi.nih.gov/blast/db)|288 Gb|[tblastn](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)|Conservation BLAST
|[eggNOG DB 5.0.2](http://eggnog5.embl.de/download/emapperdb-5.0.2/)|49-90 Gb|[eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)|Functional annotation by orthology
|[UniRef30_2023_02](https://storage.googleapis.com/alphafold-databases/v2.3/UniRef30_2021_03.tar.gz)|328 Gb|[hhsuite](https://github.com/soedinglab/hh-suite)|Search for homologous proteins

## Input Data and Configuration
### Provided Data to Reproduce Results
All data required to reproduce the analysis of six *M. tuberculosis* clinical reference strains from lineages 1 and 2 ([Heiniger et al., 2026](https://doi.org/10.64898/2026.01.27.701740)) is provided in the `data` subfolder except for the MS/MS raw data which has to be downloaded separatebly from PRIDE ([Project PXD081163](https://www.ebi.ac.uk/pride/archive/projects/PXD081163)). The downloaded `.d` files have to be moved into the strain specific folders in `data/raw/`. All other input data is provided in the following folders:

- `data/annotations/`: refseq annotation and annotations to check for presence of novels
- `data/genomes/`: genome sequence files used by the pipeline
- `data/iptgxdbs`: integrated proteogenomics databases used for searches
- `data/contams.fasta`: contaminant protein sequences used during database construction

The default strain set is defined in `config/config.yaml` under `strains` and includes the six analyzed strains: `0052`, `0072`, `0145`, `0153`, `0155`, and `0157`.

### Configuration
Before running the pipeline, check the configuration in `config/config.yaml`. Make sure that all paths to manually installed dependencies are correct and adjust the memory maximum and number of threads of each step according to your device. The provided configuration was used on a machine with 64 CPU cores and 295 Gb of RAM on Ubuntu 22.04.5 LTS.

## Running the Pipeline
Once the input data are prepared and Snakemake is installed, the analysis can be started with:

```sh
./run.sh
```

The pipeline may take a substantial amount of time depending on the number of searches configured and the available CPU resources.

After the workflow finishes, an HTML report can be generated with:

```sh
./create_report.sh
```

### Optional Analyses Requiring Manual Intervention
If enabled in `config/config.yaml`, the Phyre2 and OperonMapper analyses need manual intervention as they can be only used as a website. The pipeline creates input folders in the respective results folders that can be uploaded to the corresponding website. An output folder then has to be created next to the input folder containing the results and the pipeline has to be executed again to include these in the final tables.

|Analysis|Results Folder|Expected Input File(s)|
|--------|--------------|-------------------|
|Phyre2|`results/searches/{search}/postproc/novels/phyre2`|`summaryinfo.txt`|
|OperonMapper|`results/searches/{search}/postproc/operon_mapper`|`list_of_operons`<br>`functional descriptions`|

## Output
The workflow writes its main results under the `results` directory. The converted annotations and generated iPtgxDBs are stored in the `results/annotations` and `results/iptgxdbs` subfolders, respectively. The search results are organized by search and the 4 stages of the analysis:

- `results/searches/{search}/search/`: search against iPtgxDB assigning and rescoring peptide spectral matches
- `results/searches/{search}/entrapment/`: searches against peptide and protein-level entrapment databases to estimate actual FDR
- `results/searches/{search}/protinf/`: protein inference and filtering
- `results/searches/{search}/postproc/`: downstream annotation and post-processing of novel proteins

The generated HTML report (`report.html`) collects the workflow outputs in a more accessible form.

## Citation
**Stringent proteogenomic discovery of novel small proteins in Mycobacterium tuberculosis clinical reference strains**

Benjamin Heiniger, Christian Schori, Mohammad Arefian, Amir Banaei-Esfahani, Martin Schuler, Sonia Borrell, Chloé Loiseau, Daniela Brites, Iñaki Comas, Ruedi Aebersold, Sebastien Gagneux, Ben C. Collins, Christian H. Ahrens

bioRxiv 2026.01.27.701740; doi: https://doi.org/10.64898/2026.01.27.701740
