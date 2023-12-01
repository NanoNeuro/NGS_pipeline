Supported pipelines
--------------------
- RNAseq -> nf-core/RNAseq
- circRNA -> nf-core/RNAseq
- circDNA -> nf-core/circDNA
- taxprofiling -> own pipeline, derived from nf-core/taxprofiler
- smrnaseq -> nf-dore/smrnaseq
- scrnaseq (10X) -> nf-core/scrnaseq
- diffabundance -> own pipeline derived from nf-core/diffabundance


Samplesheet columns
-------------------
- RNAseq -> sample, fastq_1, fastq_2
- circRNA -> sample, fastq_1, fastq_2
- circDNA -> sample, fastq_1, fastq_2
- taxprofiling -> sample, fastq_1, fastq_2  [derived from nf-core/RNAseq]
- smrnaseq -> sample, fastq_1
- scrnaseq (10X) -> sample, fastq_1, fastq_2, expected_cells
- diffabundance -> sample, fastq_1, fastq_2, condition, replicate, batch


Directory structure [MIRAR CADA METODO Y ESTRUCTURA EN BASE A LA LÓGICA]
-------------------
data/
    proyect_A/
        FASTQ/
            Sample1_R1.fastq.gz
            Sample1_R2.fastq.gz
            ...
        samplesheet.csv
        project.config
database/
    genomes/
        hg38/
            genome.fasta
            genes.gtf
        mm10/
            genome.fasta
            genes.gtf
    indexes/
        hg38/
            STAR
            bowtie2
        mm10/
            bowtie2
            bowtie1
results/
    proyect_A/
        RNAseq/
        circRNA/
        taxprofiler/
processing/
    proyect_A/
        analysis1.ipynb
        analysis2.ipynb
src/
    processing_functions/
    install/
        install.sh
        conda_env.yml
README.md
.gitconfig



