import psutil


################################################################################################
# GLOBAL FUNCTIONS
################################################################################################

def get_RAM_GB():
    virtual_memory = psutil.virtual_memory()
    # Calculate and return the available RAM in gigabytes
    available_ram_gb = virtual_memory.available / (1024 ** 3)
    return available_ram_gb



def get_available_cpus():
    # Get the number of available CPUs
    available_cpus = psutil.cpu_count(logical=False)
    return available_cpus




################################################################################################
# GLOBAL VARIABLES
################################################################################################
# These pipeline have to be in the correct order of execution!
VALID_PIPELINES = ['rnaseq', 'scrnaseq', 'smrnaseq', 'circrna', 'circdna', 'taxprofiler']

NF_CONFIG_PIPELINES = ['rnaseq', 'scrnaseq', 'smrnaseq', 'circrna', 'circdna']
DEV_PIPELINES = ['circrna']

VALID_ORGANISMS = ['human', 'mouse']
DICT_GENOMES = {'human': 'GRCh38', 'mouse': 'GRCm38'}

DEFAULT_RNASEQ_ALIGNER = 'star_salmon'
DEFAULT_SCRNASEQ_ALIGNER = 'alevin'
DEFAULT_MIRGENEDB = False
DEFAULT_CIRCRNA_MODULE = 'circrna_discovery'
DEFAULT_CIRCRNA_TOOL = 'circexplorer'
DEFAULT_CIRCDNA_TOOL = 'circle_map_realign'

VALID_PROFILERS = ['kaiju', 'krakenuniq', 'kraken2', 'centrifuge']
TAXPASTA_SUMMARISE_AT = 'genus'
TAXPASTA_ADD_NAME = True
TAXPASTA_ADD_LINEAGE = True



MAX_MEMORY_PER = 0.8
MAX_CPU_PER = 0.8
MAX_TIME = 500

HISAT2_LIMIT_MEM = 200
##########


MB = 1000000
GB = 1000000000

MAX_RAM = int(MAX_MEMORY_PER * get_RAM_GB())
MAX_CPU = int(MAX_CPU_PER * get_available_cpus()) 


###########

ENSEMBL_VERSION_HUMAN = 111
ENSEMBL_VERSION_MOUSE = 102

DICT_CONTAINERS = {'ncbi-datasets': 'biocontainers/ncbi-datasets-cli:15.12.0_cv23.1.0-4', 
                   'ncbi-genome-download': 'quay.io/biocontainers/ncbi-genome-download:0.3.3--pyh7cba7a3_0', 
                   'aria2c': 'andrey01/aria2c:latest',
                   'gffread': 'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_3', 
                   'bedops': 'biocontainers/bedops:v2.4.35dfsg-1-deb_cv1', 
                   'aws-cli': 'amazon/aws-cli:2.15.5', 
                   'STAR': 'quay.io/biocontainers/star:2.7.11a--h0033a41_0',
                   'STAR-nfcore': 'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2@sha256:74c2d751f5090bdc543f8922f37ab52cc0fae1f295e73923e00a694a697535d9',
                   'bowtie': 'quay.io/biocontainers/bowtie:1.3.1--py310h7b97f60_6', 
                   'bowtie2': 'quay.io/biocontainers/bowtie2:2.5.3--py310ha0a81b8_0',
                   'bwa': 'quay.io/biocontainers/bwa:0.7.3a--he4a0461_9',
                   'hisat2': 'quay.io/biocontainers/hisat2:2.2.1--hdbdd923_6',
                   'salmon': 'quay.io/biocontainers/salmon:1.10.2--hecfa306_0',
                   'samtools': 'biocontainers/samtools:v1.9-4-deb_cv1', 
                   'kallisto': 'quay.io/biocontainers/kallisto:0.50.1--hc877fd6_0', 
                   'rsem': 'quay.io/biocontainers/rsem:1.3.3--pl526ha52163a_0', 
                   'rsem-star': 'quay.io/biocontainers/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:64aad4a4e144878400649e71f42105311be7ed87-0',
                   'cellranger': 'quay.io/nf-core/cellranger:7.1.0', 
                   'segemehl': 'quay.io/biocontainers/segemehl:0.3.4--hfe57441_9', 
                   'kaiju': 'quay.io/biocontainers/kaiju:1.10.0--h43eeafb_0', 
                   'centrifuge': 'quay.io/biocontainers/centrifuge:1.0.4--hd03093a_0', 
                   'kraken2': 'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0', 
                   'krakenuniq': 'quay.io/biocontainers/krakenuniq:1.0.4--pl5321h6dccd9a_1', 
                   'taxpasta': 'quay.io/biocontainers/taxpasta:0.6.1--pyhdfd78af_0'}
