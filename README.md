# Pipeline for NGS analysis

This is a simple pipeline that combines several NGS processing pipelines for different formats, including

* RNAseq
* circRNA and circDNA seq
* smRNAseq
* scRNAseq (10X)
* taxprofiling
* differential abundance analysis

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
* diffabundance -> own pipeline derived from nf-core/diffabundance

## Samplesheet columns

* RNAseq -> sample, fastq_1, fastq_2, strandedness
* circRNA -> sample, fastq_1, fastq_2
* circDNA -> sample, fastq_1, fastq_2
* taxprofiling -> sample, fastq_1, fastq_2  [derived from nf-core/RNAseq]
* smrnaseq -> sample, fastq_1
* scrnaseq (10X) -> sample, fastq_1, fastq_2, expected_cells
* diffabundance -> sample, fastq_1, fastq_2, condition, replicate, batch

## Directory structure

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
    processing_functions/
    notebook_functions/
    install/
        install.sh
        database_install.sh
        conda_env.yml
    run_pipeline.py
README.md
.gitconfig
```
