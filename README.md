# Pipeline for NGS analysis

This is a simple pipeline that combines several NGS processing pipelines for different formats, including

* RNAseq
* circRNA and circDNA seq
* smRNAseq
* scRNAseq (10X)
* taxprofiling
* differential abundance analysis (to be implemented)

Some of these pipelines are based on `nf-core` modules, and others are developed personally.

## Creating the environment

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

First, **projects and data do not need to share the same structure**. While it is common that one project usually uses the files of one experiment, you can use some files of one experiment, or files from different experiments in the same project. For instance, `project_A` is run with files from `Experiment_1` and `Experiment_2`, while `project_B` only uses files from  `Experiment_1`.

All project folders must contain at least two files: `samplesheet.csv` and  `config.yaml` (they don't need to be named liked that exactly).

The samplesheet file contains the information of the `data/` files that are going to be used in the pipeline.

[SEGUIR POR AQUI]

## Samplesheet columns

* RNAseq -> sample, fastq_1, fastq_2, strandedness
* circRNA -> sample, fastq_1, fastq_2
* circDNA -> sample, fastq_1, fastq_2
* taxprofiling -> sample, fastq_1, fastq_2, strandedness  [derived from nf-core/RNAseq]
* smrnaseq -> sample, fastq_1
* scrnaseq (10X) -> sample, fastq_1, fastq_2, expected_cells

[* diffabundance -> sample, fastq_1, fastq_2, condition, replicate, batch]: #
