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


Directory structure
-------------------
├── raw_data/
|   ├── Experiment_1/
|   |   ├── FASTQ/
|   |   |   ├── Exp1_Sample1_R1.fastq.gz
|   |   |   ├── Exp1_Sample1_R2.fastq.gz
|   |   |   ├── ...
|   ├── Experiment_2/
|   |   ├── FASTQ/
|   |   |   ├── Exp2_Sample1_R1.fastq.gz
|   |   |   ├── Exp2_Sample1_R2.fastq.gz
|   |   |   ├── ...
├── database/
|   ├── genomes/
|   |   ├── hg38/
|   |   |   ├── genome.fasta
|   |   |   ├── genes.gtf
|   |   ├── mm10/
|   |   |   ├── genome.fasta
|   |   |   ├── genes.gtf
|   ├── indexes/
|   |   ├── hg38/
|   |   |   ├── STAR/
|   |   |   ├── bowtie2/
|   |   ├── mm10/
|   |   |   ├── bowtie2/
|   |   |   ├── bowtie1/
├── projects/
|   ├── project_A/
|   |   ├── samplesheet.csv (Files from Exp1 and Exp2)
|   |   ├── config.yml
|   |   ├── analysis/
|   |   |   ├── nb1.ipynb
|   ├── project_B/
|   |   ├── samplesheet.csv (Files from Exp1)
|   |   ├── config.yml
├── results/
|   ├── proyect_A/
|   |   ├── RNAseq/
|   |   ├── circRNA/
|   |   ├── taxprofiler/
├── processing/
|   ├── proyect_A/
|   |   ├── figures/
|   |   ├── tables/
|   |   ├── scripts/
|   |   |   ├── analysis1.ipynb
|   |   |   ├── analysis2.py
|   |   |   ├── analysis3.R
├── src/
|   ├── processing_functions/
|   ├── notebook_functions/
|   ├── install/
|   |   ├── install.sh
|   |   ├── conda_env.yml
|   |   ├── run_pipeline.sh
├── README.md
├── .gitconfig



