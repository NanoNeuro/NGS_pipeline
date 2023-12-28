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
    available_cpus = psutil.cpu_count(logical=True)
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


##########


MB = 1000000
GB = 1000000000

MAX_RAM = int(MAX_MEMORY_PER * get_RAM_GB())
MAX_CPU = int(MAX_CPU_PER * get_available_cpus()) 


###########


DICT_CONTAINERS = {'ncbi-datasets': 'biocontainers/ncbi-datasets-cli:15.12.0_cv23.1.0-4', 
                   'gffread': 'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_3', 
                   'bedops': 'biocontainers/bedops:v2.4.35dfsg-1-deb_cv1', 
                   'aws-cli': 'amazon/aws-cli:2.15.5', 
                   'STAR': 'hydragenetics/star:2.7.10a',
                   'bowtie': 'biocontainers/bowtie:v1.2.2dfsg-4-deb_cv1', 
                   'bowtie2': 'biocontainers/bowtie2:v2.4.1_cv1',
                   'samtools': 'biocontainers/samtools:v1.9-4-deb_cv1'}
