# Pipeline for NGS analysis

This is a simple pipeline that combines several NGS processing pipelines for different formats, including

* RNAseq
* circRNA and circDNA seq
* smRNAseq
* scRNAseq (10X)
* taxprofiling
* differential abundance analysis (to be implemented)

Some of these pipelines are based on `nf-core` modules, and others are developed personally.

## Environment

In order to be able to run the pipelines, some libraries have to be installed first. To do that, we are going to build a `conda` environment that installs all the requirements, as well as some dependencies (like aws cli) and checks of the nf-core pipelines.

To do that, we are going to run the command:

```bash
bash src/install_env.sh
```

Any additional software, like centrifuge, will be installed in `/home/USER/Programs` you can change the directory where
these programs will be installed in the variable `DIR_PROGRAMS` of the `install_env.sh` file.

## Supported pipelines

* RNAseq -> nf-core/RNAseq
* circRNA -> nf-core/RNAseq
* circDNA -> nf-core/circDNA
* taxprofiling -> own pipeline, derived from nf-core/taxprofiler
* smrnaseq -> nf-dore/smrnaseq
* scrnaseq (10X) -> nf-core/scrnaseq

[* diffabundance -> own pipeline derived from nf-core/diffabundance]: # 


## Directory structure and project structuring philosophy

The directory structure is divided in several main folders:

* `src/` contains all the scripts, python files, etc. that are used to run the pipeline. Any script that is specific of a project should NOT be put here. However, scripts that are meant to be used by more than one pipeline should be saved in `src/shared_scripts/`.
* `database/` contains all base files that are used for reference/index building (for instance, genome/transcriptome `.fa`, `.gtf`), as well as built indexes. This folder is divided in subfolders:
  * `genomes/` contains all "basic" files to build indexes.
  * `indexes/` has already built indexes (for instance, STAR, bowtie2, kallisto, etc.).
  * `smRNA/` contains files/indexes used for smRNAseq pipeline.
  * `taxprofiler/` contains indexes related to taxprofiling pipelines. Basic files to build the indices should be stored, if necessary, in `genomes/`.
* `data/` contains all the raw data generated during an experiment (`.fastq`, `.bam`, `.bcl2`, etc.). It is recommended that each folder has a spreadsheet with information regarding the base files (organism, age, sex, condition, QC, etc.).
* `projects/` contains files that will be used by the NGS_pipeline to run. It also contains scripts/notabooks that are used to analyse the results. How projects are organised is explained below.
* `results/` contains the processed files of each project. In fact, each subfolder's name is the name of the project it belongs to.
* `work/` contains intermediate files generated during the run of NGS_pipeline. This directory can be removed after the pipeline has run.

This is an example of the directory structure, with two experiments, database, and two projects.

```bash
data/
    Experiment_1/
        FASTQ/
            Exp1_Sample1_R1.fastq.gz
            Exp1_Sample1_R2.fastq.gz
            ...
        fastq_metadata.csv   # Subject, Age, Sex, QC, etc.
    Experiment_2/
        FASTQ/
            Exp2_Sample1_R1.fastq.gz
            Exp2_Sample1_R2.fastq.gz
            ...
        fastq_metadata.csv
database/
    genomes/
        GRCh38/
            genome.fasta
            genes.gtf
            genes.bed
        GRCm38/
            genome.fasta
            genes.gtf
            genes.bed
    smRNA/
        mirbase/
            mature.fa
            ...
        mirgenedb/
            mature.fa
            ...
    indexes/
        GRCh38/
            STAR
            bowtie2
        GRCm38/
            bowtie2
            bowtie1
    taxprofiler/
        kaiju/
        kraken2/
        krakenuniq/
        centrifuge/
        taxonomy/
projects/
    project_A/
        samplesheet.csv   # Files from Exp1 and Exp2 
        config.yml
        analysis/
            analysis1.ipynb
            analysis2.py
            analysis3.R
        figures/
        tables/
    project_B/
        samplesheet.csv   # Files from Exp1
        config.yml
results/
    project_A/
        RNAseq/
        circRNA/
        taxprofiler/
work/
    project_A/
        ...
src/
    shared_scripts/
    processing_functions/
    notebook_functions/
    config/
        install.sh
        conda_env.yaml
        taxpasta_exclusions.yaml
    run_pipeline.py
README.md
.gitconfig
```

## How to organize each project

First, **projects and data do not need to share the same structure**. While it is common that one project usually uses the files of one experiment, you can use some files of one experiment only, or files from different experiments in the same project. For instance, `project_A` is run with files from `Experiment_1` and `Experiment_2`, while `project_B` only uses files from  `Experiment_1`.

All project folders must contain at least two files: `samplesheet.csv` and  `config.yaml` (they don't need to be named liked that exactly).

The samplesheet file contains the information of the `data/` files that are going to be used in the pipeline. 

The `config.yaml` file contains all the information about parameters that are necesary for each process. Also, you can run more than one process on the same project (for instance, if you want to compare results from different aligners in an rnaseq experiment).

### Preparing the `config.yaml` file
The `config.yaml` file is organised in **processes**. A process is a run of a pipeline with a specific set of parameters. Here you can see one example:

```yaml
EXAMPLE_RNASEQ
  general_config:
    pipeline: rnaseq
    organism: mouse
    max_cpus: 3
  nextflow_config:
    resume: false
    profile: docker
  nfcore_config:
    genome: GRCm38
    aligner: star_salmon
  prescript:
    "bash projects/project_A/rnaseq_preprocessing.sh"
  postscript:
    "bash projects/project_A/rnaseq_postprocessing.sh"

EXAMPLE_TAXPROF
  general_config:
    pipeline: taxprofiler
    organism: human
  profiler_config:
    profilers: kaiju,kraken2
  host_mapping_config:
    1st_map: true
    1st_map_extra_args: "--min_trimmed_reads 250"
    2nd_map: true
  kaiju_config:
    minimum_length: 39
    kaiju_extra_args: "-a mem"
```

In this example two processes are run: `EXAMPLE_RNASEQ` and `EXAMPLE_TAXPROF`.

All processes have a 3 common attributes: `general_config`, `prescript` and `postscript`.
* `general_config` includes `organism` (`human` or `mouse`) and `pipeline` (`rnaseq`, `scrnaseq`, `smrnaseq`, `circrna`, `circdna`, `taxprofiler`).
* `prescript` and `postcript` are text attributes that are included before and after the script of the pipeline in `command.sh`.

There are also specific attributes in nf-core-related pipelines:
* `nextflow_config`: it involves general nextflow attributes, which are followed by one hyphen (`-`) in the nf-core pipelines.
* `nfcore_config`: it involves pipeline-specific attributes, which are followed by two hyphens (`--`) in the nf-core pipelines.

Regarding attribute from `taxprofiler` pipeline, they are the following:
* `profiler_config`: the main argument here is `profilers`, which is the list of profilers that can be run (`kaiju`, `centrifuge`, `kraken2`, `krakenuniq`).
* `host_mapping_config`: it includes attributes related to performing 1 or 2 maps to host genome. 1st map is done by running nf-core/rnaseq with star_salmon aligner; and 2nd map is done by running bowtie2 with CHM13 genome. Currently, only human host mapping is supported.
* `<PROFILER>_config`: specific configurations of each profiler.

To check more specific values, check the `projects/test_human/config.yaml` and `projects/test_mouse/config.yaml`.


### Preparing the `samplesheet.csv` file
Depending on the pipeline, you will need to follow a pattern. The section [Samplesheet columns](#samplesheet-columns) contains all the columns that can be used in each pipeline. You don't need, however, to complete all the columns, but only the basic ones (generally `sample`, `fastq_1` and `fastq_2`). For more info, check on the nf-core pipelines.

Bear in mind that **columns must be comma-separated**.

In the `samplesheet.csv` file you'll need to add another column: `process`. This refers to the process that files are assigned to in the `config.yaml` file. One example of `samplesheet.csv` would be:

| sample | fastq_1                  | fastq_2                  | process            |
|--------|--------------------------|--------------------------|--------------------|
| S1     | data/Exp1/S1_R1.fastq.gz | data/Exp1/S1_R2.fastq.gz | PROC_ALL           |
| S2     | data/Exp1/S2_R1.fastq.gz | data/Exp1/S2_R2.fastq.gz | PROC_ALL;PROC_S2S3 |
| S3     | data/Exp1/S3_R1.fastq.gz | data/Exp1/S3_R2.fastq.gz | PROC_ALL;PROC_S2S3 |

In this case, the process `PROC_ALL` would use all three samples in the pipeline, while `PROC_S2S3` would only use S2 and S3. Processes have to be separated with a **semicolon** (`;`) and no separation between processes. 

**The names of the files must be relative to the root dir.**

## How to run the code
**First move to the source directory** (the one from which you can see `src`, `data`, `database`, etc.).

If you have a prepared conda environment (look at [Environment](#environment)), make sure it is activated.

Then run the command

```
python src/run_pipeline.py --project PROJECT --yaml YAML --samplesheet SAMPLESHEET
```

* `--project` is the name of the folder in `projects/` that you are going to run. This is the only **mandatory** argument.
* `--yaml` is optional, and set to `config.yaml` by default. It refers to the name in the `projects/` folder.
* `--samplesheet` is optional and set to `samplesheet.csv` by default. 

### Intermediate files
Once you run the code, there are some files that are generated in the `work/` directory so that you can check if anything is missed in the execution. 

* First, for each process within a project, the file `work/<PROJECT>/samplesheets/<PROCESS>.csv` stores the information of files that will be used. For nf-core pipelines, this file is the one used as an input.

* `work/<PROJECT>/database_download.sh` is a bash script that downloads and generates all the necessary databases to run all processes within the project. 

* `work/<PROJECT>/command.sh` is the bash script that runs all pipelines.

In case you need to make some quick modifications, you can change these files and then run them as normal.


## Samplesheet columns

* RNAseq -> sample, fastq_1, fastq_2, strandedness
* circRNA -> sample, fastq_1, fastq_2
* circDNA -> sample, fastq_1, fastq_2
* taxprofiling -> sample, fastq_1, fastq_2, strandedness  [derived from nf-core/RNAseq]
* smrnaseq -> sample, fastq_1
* scrnaseq (10X) -> sample, fastq_1, fastq_2, expected_cells

[* diffabundance -> sample, fastq_1, fastq_2, condition, replicate, batch]: #
