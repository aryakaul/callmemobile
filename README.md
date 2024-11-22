# callmemobile – annotation of putatively mobile bacterial genomic regions
![GitHub Release](https://img.shields.io/github/v/release/baymlab/callmemobile)


Workflow for identification of putatively mobile genetic regions from genome
assemblies and bed files.
callmemobile first runs various prediction tools to classify mobile features
of the input genomes and then compiles these results alongside user options 
and a set of genetic regions to identify which subset of these regions are 
predicted to be mobile (maybe).

If you use this in your own work, please cite the papers whose tools the workflow uses:
* [Phigaro](https://doi.org/10.1093/bioinformatics/btaa250)
* [MOB-suite](https://doi:10.1099/mgen.0.000206)
* [plasmidfinder](https://doi:10.1128/AAC.02412-14)
* [IntegronFinder](https://doi:10.3390/microorganisms10040700)

As well as, the paper for which the workflow was written for: [link preprint] 


<h2>Contents</h2>

<!-- vim-markdown-toc GFM -->

* [1. Introduction](#1-introduction)
* [2. Dependencies](#2-dependencies)
* [3. Installation](#3-installation)
* [4. Usage](#4-usage)
    * [4a. Basic example](#4a-basic-example)
    * [4b. Adjusting configuration](#4b-adjusting-configuration)
    * [4c. List of workflow commands](#4c-list-of-workflow-commands)
    * [4d. Troubleshooting](#4d-troubleshooting)
* [5. Citation](#5-citation)
* [6. Issues](#6-issues)
* [7. Changelog](#7-changelog)
* [8. License](#8-license)
* [9. Contacts](#9-contacts)

<!-- vim-markdown-toc -->


## 1. Introduction

The user provides 2 files in the `input/` directory. The first is a file of 
files containing each genomic assembly to be analyzed and their respective
paths. This file ends with a `.txt`. The second is a file of files containing
the bed files of genomic regions to be queried. This file ends with a `.beds`. 
They shoudl share the same basename (i.e. `nctc3k.txt` and `nctc3k.beds`) and
should have the same number of lines. Order of the inputs is how we match the 
bed file to the genomic assembly, so please ensure they are matched correctly.


## 2. Dependencies

* [Conda](https://docs.conda.io/en/latest/miniconda.html) (unless the use of Conda is switched off in the configuration) and ideally also [Mamba](https://mamba.readthedocs.io/) (>= 0.20.0)
* [GNU Make](https://www.gnu.org/software/make/)
* [Python](https://www.python.org/) (>=3.7)
* [Snakemake](https://snakemake.github.io) (>=8.0.0)

and can be installed by Conda by
```bash
conda install -c conda-forge -c bioconda -c defaults \
  make "python>=3.7" "snakemake-minimal>=8.0.0" "mamba>=0.20.0"
```


The rest of the dependencies are installed automatically by Snakemake
when they are requested. The specifications of individual environments
can be found in [`workflow/envs/`](workflow/envs/),
and they contain:
- bedops
- phigaro
- Biopython 
- blast
- cgecore
- csvtk
- seqkit
- plasmidfinder
- MOB-suite
- mobileelementfinder
- pandas

All dependencies across all rules can also be
installed at once by `make conda`.


## 3. Installation

Clone and enter the repository by

```bash
git clone https://github.com/baymlab/callmemobile
cd callmemobile
```

Alternatively, the repository can also be installed using cURL by
```bash
mkdir callmemobile
cd callmemobile
curl -L https://github.com/baymlab/callmemobile/tarball/main \
    | tar xvf - --strip-components=1
```


## 4. Usage

### 4a. Basic example

* ***Step 1: Provide lists of input files.*** \
  For every input, create a txt list of input assemblies in the `input/`
  directory (i.e., as `input/{batch_name}.txt`. Use either absolute paths (recommended),
  or paths relative to the root of the Github repository (not relative to the txt files).

  Such a list can be generated, for instance, by `find` by
  ```bash
  find ~/dir_with_my_genomes -name '*.fa' > input/my_first_batch.txt
  ```
  The supported input files should be in FASTA format.

  Next, create another txt list of input beds in the `input/` directory.
  ```bash
  find ~/dir_with_my_beds -name '*.bed' > input/my_first_batch.beds
  ```

  These should follow the format of typical bed files:
  ```bash
   ➜ head NCTC10036-assembly_abr-search.abr.bed
  ENA|LR134493|LR134493.1_388     406996  410142  adeF
  ENA|LR134493|LR134493.1_694     718039  718224  rsmA
  ENA|LR134493|LR134493.1_1759    1957864 1958496 CRP
  ```
  
The first column corresponds to the contig, the second column is the position
start coordinate and the third column is the end coordinate. The fourth column 
can be any string description of the region (in this case, the gene name).

* ***Step 2 (optional): Adjust configuration.*** \
  By editing [`config.yaml`](config.yaml) it is possible to specify
  parameters of the run, both in terms of options for
  individual tools and the heuristics employed to deem
  a region as potentially mobile.

* ***Step 3: Run the pipeline.*** \
  Run the pipeline by running `make all`; this is run
  Snakemake with the corresponding parameters.

* ***Step 4: Retrieve the output files.*** \
  All output files will be located in `output/`.


### 4b. Adjusting configuration

The workflow can be configured via the [`config.yaml`](./config.yaml) file, and
all options are documented directly there. The configurable functionality includes:

- prophage_maxdist: base pair distance from a predicted prophage to classify it as potentially mobile
- integron_pctolap: percent overlap of given region with an integron to classify it as potentially mobile
- mobileelement_maxdist:    distance between two mobile elements of the same type to classify it as potentially mobile
- plasmidfinder_mincov:     minimum coverage for plasmidfinder to classify a contig as a possible plasmid
- plasmidfinder_threshold:  minimum threshold for plasmidfinder to classify a contig as a possible plasmid

### 4c. List of workflow commands

callmemobile is executed via [GNU Make](https://www.gnu.org/software/make/), which handles all parameters and passes them to Snakemake.
Here's a list of all implemented commands (to be executed as `make {command}`):


```yaml
######################
## General commands ##
######################
    all                  Run everything (the default subcommand)
    help                 Print help messages
    conda                Create the conda environments
    clean                Clean all output archives and files with statistics
    cleanall             Clean everything but Conda, Snakemake, and input files
    cleanallall          Clean completely everything
###############
## Reporting ##
###############
    viewconf             View configuration without comments
    reports              Create html report
####################
## For developers ##
####################
    test                 Run the workflow on test data (P1)
    bigtest              Run the workflow on test data (P1, P2, P3)
    format               Reformat all source code
    checkformat          Check source code format
```

*Note:* `make format` and `make checkformat` require
[YAPF](https://github.com/google/yapf) and
[Snakefmt](https://github.com/snakemake/snakefmt), which can be installed by
`conda install -c conda-forge -bioconda yapf snakefmt`.


### 4d. Troubleshooting

Tests can be run by `make test`. 


## 5. Citation

[todo]

```bibtex
todo
```


## 6. Issues

Please use [Github issues](https://github.com/karel-brinda/miniphy/issues).



## 7. Changelog

See [Releases](https://github.com/karel-brinda/miniphy/releases).



## 8. License

[MIT](https://github.com/karel-brinda/miniphy/blob/master/LICENSE)



## 9. Contacts

* [Arya Kaul](https://arya.casa) \<arya_kaul@g.harvard.edu\>

