# *M. tuberculosis* proteogenomics Snakemake pipeline
This Snakemake pipeline performs a proteogenomic analysis of six *M. tuberculosis* clinical reference strains from lineages 1 and 2 ([Heiniger et al., 2026](https://doi.org/10.64898/2026.01.27.701740)). The workflow combines genome annotation processing, peptide-spectrum matching, protein inference, stringent filtering and downstream annotation and structural analyses to identify and prioritize novel proteins.

The pipeline relies on three integrated proteogenomics databases ([iPtgxDBs](https://iptgxdb.expasy.org/)) per strain: a RefSeq iPtgxDB including only the RefSeq annotation, a comprehensive standard iPtgxDB ([Omasits et al., 2017](https://doi.org/10.1101/gr.218255.116)) and a much smaller custom iPtgxDB based on Ribo-seq data ([Hadjeras et al., 2023](https://doi.org/10.1093/femsml/uqad012)), here using orthologs of Ribo-seq identifications from strain H37Rv ([Smith et al. 2022](https://doi.org/10.7554/eLife.73980)). Mass spectrometry data from a Bruker timsTOF Pro device ([PRIDE Project PXD081163](https://www.ebi.ac.uk/pride/archive/projects/PXD081163)) is searched for each strain against the respective iPtgxDBs using [MSFragger](https://msfragger.nesvilab.org/), relying on [MSBooster](https://doi.org/10.1038/s41467-023-40129-9) and [Percolator](https://github.com/percolator/percolator) for accurate peptide identification. Protein inference and false discovery rate (FDR) filtering is achieved using [Philosopher](https://philosopher.nesvilab.org/), followed by additional stringent filtering as described ([Heiniger et al., 2026](https://doi.org/10.64898/2026.01.27.701740)), including the selection of unambiguous peptides based on the PeptideClassifier software ([Qeli and Ahrens, 2010](https://doi.org/10.1038/nbt0710-647)) extended for prokaryotic proteogenomics ([Omasits et al., 2017](https://doi.org/10.1101/gr.218255.116)). The FDR among all proteins and the pre-defined subsets of RefSeq and novel proteins is estimated with an entrapment analysis using FDRBench ([Wen et al., 2025](https://doi.org/10.1038/s41592-025-02719-x)).

Novel genes not annotated before, corrected start sites of annotated genes and pseudogenes with confirmed expression are reported separately. Novel gene candidates are further analyzed, including the overlap with newer annotations and novel CDS detected with ribosome profiling in strain H37Rv ([Smith et al., 2022](https://doi.org/10.7554/eLife.73980)), with analogous prediction of a signal of selection. Further, potential functions based on structural similarity and the inclusion in operons, as well as the conservation at different taxonomic levels are predicted and integrated with [Panaroo](https://github.com/gtonkinhill/panaroo)-based pangenome information. Proteins newly identified with this pipeline include conserved and lineage-specific SEPs, an antitoxin, candidate antimicrobial peptides and novel proteins under purifying selection.

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
  - [Most Important Output Files](#most-important-output-files)
  - [iPtgxDBs Folder](#iptgxdbs-folder)
  - [Searches Folder](#searches-folder)
    - [Search](#search)
    - [Protein Inference](#protein-inference)
    - [Entrapment](#entrapment)
    - [Post-Processing](#post-processing)
- [Rulegraph](#rulegraph)
- [Citation](#citation)

## Software Dependencies
### Automatically Installed
By default, dependencies are installed automatically using conda when running the pipeline for the first time, except for FragPipe and dependencies used for optional downstream analyses of the identified novel proteins. The path of each manually installed software has to be configured in `config/config.yaml`. To disable automatic dependency management, remove `--use-conda` from `run.sh` and make sure the required tools are installed manually.

If PSORTb is enabled, make sure to [install Docker](https://docs.docker.com/engine/install) and also enable [management of Docker as a non-root user](https://docs.docker.com/engine/install/linux-postinstall).

### Manually Installed
#### Mandatory
FragPipe, MSFragger and Philosopher have to be installed manually. While MSFragger can be installed from within the FragPipe GUI, Philosopher should be installed manually as a newer version was used here to calculate protein-level q-values. The following versions were used for the pipeline:

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
|[BatMass-io](https://batmass.org)|1.33.4|

After installation, the paths to the 6 tools mentioned above have to be adapted in `config/config.yaml`.

#### Optional
These dependencies have to be installed and the path adjusted in `config/config.yaml` or the analysis has to be disabled by setting `enable: False` in the corresponding config section.

Either the software itself has be installed manually:

|Software|Version|Size|Notes|
|--------|------:|---:|-----|
|[FDRBench](https://github.com/Noble-Lab/FDRBench)|0.0.4|74 Mb|Estimation of FDR for a pre-defined subset|
|[AMP-scanner](https://github.com/dan-veltri/amp-scanner-v2)|2|28 Mb|Antimicrobial peptide prediction|
|[InterProScan](https://github.com/ebi-pf-team/interproscan)|5.59-91.0|47 Gb|Identification of functional protein domains|
|[LipoP](https://services.healthtech.dtu.dk/services/LipoP-1.0/)|1.0a|632 Kb|Lipoprotein and signal peptide detection|

Or the software is automatically installed with conda but a dataset has to be downloaded manually:

|Dataset|Size|Software|Notes|
|-------|---:|-------:|-----|
|[core_nt](https://ftp.ncbi.nih.gov/blast/db)|288 Gb|[tblastn](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)|Conservation BLAST
|[eggNOG DB 5.0.2](http://eggnog5.embl.de/download/emapperdb-5.0.2/)|49-90 Gb|[eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)|Functional annotation by orthology
|[UniRef30_2023_02](https://wwwuser.gwdguser.de/~compbiol/uniclust/2023_02/UniRef30_2023_02_hhsuite.tar.gz)|328 Gb|[hhsuite](https://github.com/soedinglab/hh-suite)|Search for homologous proteins

## Input Data and Configuration
### Provided Data to Reproduce Results
All data required to reproduce the analysis of six *M. tuberculosis* clinical reference strains from lineages 1 and 2 ([Heiniger et al., 2026](https://doi.org/10.64898/2026.01.27.701740)) is provided in the `data` subfolder except for the MS/MS raw data which has to be downloaded separately from PRIDE ([Project PXD081163](https://www.ebi.ac.uk/pride/archive/projects/PXD081163)). The downloaded `.d` files have to be moved into the strain specific folders in `data/raw/`. All other input data is provided in the following folders:

|Path|Description|
|------|-----------|
|`data/annotations`|RefSeq annotation and annotations to check for presence of novel proteins|
|`data/genomes`|Genome sequence files used by the pipeline|
|`data/iptgxdbs`|Integrated proteogenomics databases (iPtgxDBs) used for searches|
|`data/contams.fasta`|Contaminant proteins from CRAPome ([Mellacheruvu et al., 2013](https://doi.org/10.1038/nmeth.2557))|

The default strain set is defined in `config/config.yaml` under `strains` and includes the analyzed strains from lineage 1 (`0072`, `0153`, and `0157`) and lineage 2 (`0052`, `0145` and `0155`).

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
If enabled in `config/config.yaml`, the Phyre2 and OperonMapper analyses both need manual intervention, as they can be only used as a website. The pipeline creates input folders in the respective results folders that can be uploaded to the corresponding website. An output folder then has to be created next to the input folder containing the results and the pipeline has to be executed again to include these in the final tables.

|Analysis|Results Folder|Expected Input File(s)|
|--------|--------------|-------------------|
|Phyre2|`results/searches/{search}/postproc/novels/phyre2`|`summaryinfo.txt`|
|OperonMapper|`results/searches/{search}/postproc/operon_mapper`|`list_of_operons`<br>`functional descriptions`|

## Output
The workflow writes all of its results into the `results` directory, which is divided into two subfolders: The `iptgxdbs` subfolder contains the processed integrated proteogenomics databases ([iPtgxDBs](https://iptgxdb.expasy.org/)) against which the mass spectrometry data is searched and the `searches` subfolder contains the search results including protein inference and all post-processing steps.

### Most Important Output Files
The subfolder for each search includes the following files which contain the most important results:

|File|Description|
|----|-----------|
|`protinf/protinf_summary.tsv`|Summary of identified proteins per strain and type of novelty|
|`protinf/split_cats/{strain}/{subset}.tsv`|Lists of identified proteins, separated by the type of novelty|
|`postproc/add_annot/{strain}/{subset}.tsv`|Identified RefSeq proteins and novel starts with annotation information added|
|`postproc/novels/combined/{strain}.tsv`|Novel proteins per strain including all post-processing analyses|
|`postproc/novels/novels_panaroo_postproc.tsv`|Orthogroups of identified novel proteins including all post-processing analyses|

If entrapment is enabled, the FDP curves of each sample can be found in this location:
`entrapment/fdrbench/post_protparam/plot_fdp/{strain}/{sample}/{subset}.png`

### iPtgxDBs Folder
The `results/iptgxdbs` folder contains a subfolder for each iPtgxDB type which differs in the composition of the included annotations (see table below). For each strain, contaminants from the CRAPome dataset ([Mellacheruvu et al., 2013](https://doi.org/10.1038/nmeth.2557)) are added to the database and the physico-chemical parameters of the proteins are calculated.


|Annotation<br>Source|RefSeq<br>iPtgxDB|Standard<br>iPtgxDB|Custom<br>iPtgxDB|
|----------|------|--------|------|
|RefSeq|✔|✔|✔|
|Prodigal| |✔|✔|
|ChemGenome| |✔|✔|
|*in silico* ORFs| |✔| |
|H37Rv Ribo-seq| | |✔|

### Searches Folder
For each search defined in `config/config.yaml`, a corresponding subfolder is generated in `results/searches`. The results are then further categorized into the 4 stages of the analysis:

#### Search
The `search` subfolder contains the results of the proteomics searches with MSFragger, including the search database converted to the required format and peptide rescoring based on MSBooster and Percolator.

|Subfolder|Description|
|---------|-----------|
|`database`|iPtgxDBs converted to format compatible with MSFragger|
|`msfragger`|Results of the MSFragger searches|
|`msbooster`|More accurate prediction of *in silico* spectra with MSBooster|
|`percolator`|Rescoring of peptide identifications with Percolator|

#### Protein Inference
The `protinf` subfolder contains the results of inferring proteins using ProteinProphet, the FDR filtering and consolidation of results from multiple samples with Philosopher and further custom filtering to limit FDR and prevent ambiguous identifications.

|Subfolder|Description|
|---------|-----------|
|`proteinprophet`|Protein inference using ProteinProphet|
|`filter_and_report`|FDR filtering using Philosopher|
|`abacus`|Consolidation of results from multiple samples with Abacus|
|`pepclass_filter`|Filtering of ambiguous protein identifications with PeptideClassifier|
|`protparam_filter`|Ad hoc filtering of proteins based on number of PSMs and peptides|
|`collapse_extensions`|Collapse multiple predictions of novel starts for the same gene|
|`pseudo_overlap`|Reclassify novel proteins as pseudogenes if they have large in-frame overlaps|
|`split_cats`|Split identified proteins into RefSeq and different types of novelty|

#### Entrapment
If entrapment is enabled an `entrapment` folder will be created which also contains `search` and `protinf` subfolders for the respective entrapment searches, matching the descriptions for the normal searches above. However, it additionally contains an `fdrbench` subfolder with the following contents:

|Subfolder|Description|
|---------|-----------|
|`pre_protparam/calc_fdp`|FDRBench based FDP calculation before ad hoc (aka protparam) filter|
|`pre_protparam/plot_fdp`|Plots of FDP vs. FDR curves before ad hoc (aka protparam) filter|
|`post_protparam/calc_fdp`|FDRBench based FDP calculation after ad hoc (aka protparam) filter|
|`post_protparam/plot_fdp`|Plots of FDP vs. FDR curves after ad hoc (aka protparam) filter|

#### Post-Processing
The `postproc` subfolder contains all post-processing that has been performed either on all identified proteins or the subset of novel proteins.

|Subfolder|Description|
|---------|-----------|
|`add_annot`|Add details from annotation to identified proteins|
|`operon_mapper`|Predict operons based on identified RefSeq and novel proteins|
|`panaroo`|Calculate pangenome based on identified RefSeq and novel proteins|
|`novels`|Subfolders of all downstream analyses on novel proteins (see table below)|

The `novels` folder may contain the following subfolders, depending on which downstream analyses were enabled in `config/config.yaml`

|Subfolder|Description|
|---------|-----------|
|`fdp_filter`|Entrapment based classification into high and low confidence novel proteins|
|`annot_overlap`|Overlap with annotations not included in iPtgxDBs|
|`extract_seqs`|Extract protein sequences of novel proteins for further analyses|
|`conservation_blast`|Assess conservation at different taxonomic levels|
|`lipop`|Lipoprotein and signal peptide detection|
|`psortb`|Prediction of subcellular localization|
|`eggnog`|Functional annotation by orthology|
|`interpro`|Identification of functional protein domains|
|`hhsuite`|Search for homologous proteins|
|`phyre2`|Prediction of structure and function|
|`codon_gc_freq`|Signal of selection prediction based on codon GC skew ([Smith et al. 2022](https://doi.org/10.7554/eLife.73980))|
|`ampscanner2`|Antimicrobial peptide prediction|
|`combined`|Integration of all of the analyses above into a summary table

## Rulegraph
This graph shows the dependencies of the defined Snakemake rules. Arrows indicate that the rule from which the arrow originates produces the files that are used as input by the rule the arrow points to.

![Rulegraph](rulegraph.svg)

## Citation
**Stringent proteogenomic discovery of novel small proteins in *Mycobacterium tuberculosis* clinical reference strains**

Benjamin Heiniger, Christian Schori, Mohammad Arefian, Amir Banaei-Esfahani, Martin Schuler, Sonia Borrell, Chloé Loiseau, Daniela Brites, Iñaki Comas, Ruedi Aebersold, Sébastien Gagneux, Ben C. Collins, Christian H. Ahrens

bioRxiv 2026.01.27.701740; doi: https://doi.org/10.64898/2026.01.27.701740