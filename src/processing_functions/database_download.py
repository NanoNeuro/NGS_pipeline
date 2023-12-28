import logging
import os

from .versions import VERSION
from .global_vars_and_funcs import MAX_CPU, MAX_RAM, MB, GB, DICT_CONTAINERS



logger = logging.getLogger()




################################################################################################
# PRIMARY FUNCTIONS
################################################################################################
def database_download(list_dbs, project):

    file_text =  f"""
echo "
    _____                                              
   /  /::\        _____                                
  /  /:/\:\      /  /::\                               
 /  /:/  \:\    /  /:/\:\                              
/__/:/ \__\:|  /  /:/~/::\                             
\  \:\ /  /:/ /__/:/ /:/\:|                            
 \  \:\  /:/  \  \:\/:/~/:/                            
  \  \:\/:/    \  \::/ /:/                             
   \  \::/      \  \:\/:/                              
    \__\/        \  \::/                               
                  \__\/                                
    _____          ___           ___           ___     
   /  /::\        /  /\         /__/\         /__/\    
  /  /:/\:\      /  /::\       _\_ \:\        \  \:\   
 /  /:/  \:\    /  /:/\:\     /__/\ \:\        \  \:\  
/__/:/ \__\:|  /  /:/  \:\   _\_ \:\ \:\   _____\__\:\ 
\  \:\ /  /:/ /__/:/ \__\:\ /__/\ \:\ \:\ /__/::::::::\ 
 \  \:\  /:/  \  \:\ /  /:/ \  \:\ \:\/:/ \  \:\~~\~~\/
  \  \:\/:/    \  \:\  /:/   \  \:\ \::/   \  \:\  ~~~ 
   \  \::/      \  \:\/:/     \  \:\/:/     \  \:\     
    \__\/        \  \::/       \  \::/       \  \:\    
                  \__\/         \__\/         \__\/    
                   ___           ___          _____    
                  /  /\         /  /\        /  /::\   
                 /  /::\       /  /::\      /  /:/\:\  
 ___     ___    /  /:/\:\     /  /:/\:\    /  /:/  \:\ 
/__/\   /  /\  /  /:/  \:\   /  /:/~/::\  /__/:/ \__\:|
\  \:\ /  /:/ /__/:/ \__\:\ /__/:/ /:/\:\ \  \:\ /  /:/
 \  \:\  /:/  \  \:\ /  /:/ \  \:\/:/__\/  \  \:\  /:/ 
  \  \:\/:/    \  \:\  /:/   \  \::/        \  \:\/:/  
   \  \::/      \  \:\/:/     \  \:\         \  \::/   
    \__\/        \  \::/       \  \:\         \__\/    
                  \__\/         \__\/                    

------------------------------------------------------------------
| NGS PIPE MERGE v{VERSION} | DATABASE DOWNLOAD
------------------------------------------------------------------

------------------------------------------------------------------
| Databases to check: {', '.join(list_dbs)}
------------------------------------------------------------------
"
"""



    for db in list_dbs:
        path_db, func_db = DICT_DBS[db]
        text_db = func_db(path_db, file_text)
        file_text += text_db



    write_database_file(file_text, project)


# TODO: EN CADA UNA DE LAS BASES DE DATOS, COMPROBAR QUE ESTÁN LOS ARCHIVOS NECESARIOS Y QUE NO ESTÁN VACIOS
# TODO: probar a hacer la descarga desde https://quay.io/repository/biocontainers/rsem?tab=tags y correr el proceso que sea desde docker

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/ # HUMAN
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/ # MOUSE
 # TODO MIRAR CURL
    # https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.40/download?include_annotation_type=FASTA,GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT
    # return ""




################################################################################################
# SECONDARY FUNCTIONS
################################################################################################
def download_fasta_GRCh38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=3 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

        
    path, filename = os.path.split(path_db, file_text)

    title = """echo "### Downloading GRCh38 genome fasta file."\n\n"""

    cmd = f""" docker run -it --rm -v $(pwd)/{path}/:/DB_DOWNLOAD \\\n\
        {DICT_CONTAINERS['ncbi-datasets']} \\\n\
        sh -c "\\\n\
            datasets download genome accession GCF_000001405.40 --include genome;\\\n\
            unzip ncbi_dataset.zip;\\\n\
            cp ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna /DB_DOWNLOAD/{filename}\\\n\
    "\n\n\n"""


    cmd_index_fai = f""" docker run -it --rm -v $(pwd)/{path}/:/DB_DOWNLOAD \\\n\
        {DICT_CONTAINERS['samtools']} \\\n\
        sh -c "\\\n\
            samtools faidx /DB_DOWNLOAD/{filename}\\\n\
    "\n\n\n"""  
    
    
    return title + cmd + cmd_index_fai

   

def download_transcript_fasta_GRCh38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=380 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

        
    path, filename = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)
    
    title = """echo "### Downloading GRCh38 transcript fasta file."\n\n"""

    cmd = f""" docker run -it --rm -v $(pwd)/{path}/:/DB_DOWNLOAD \\\n\
        {DICT_CONTAINERS['ncbi-datasets']} \\\n\
        sh -c "\\\n\
            datasets download genome accession GCF_000001405.40 --include cds;\\\n\
            unzip ncbi_dataset.zip;\\\n\
            cp ncbi_dataset/data/GCF_000001405.40/cds_from_genomic.fna /DB_DOWNLOAD/{filename}\\\n\
    "\n\n\n"""
    
    return title + cmd



def download_gtf_GRCh38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=500 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

        
    path, filename = os.path.split(path_db, file_text)
    os.makedirs(path, exist_ok=True)

    title = """echo "### Downloading GRCh38 gene gtf file."\n\n"""

    cmd = f""" docker run -it --rm -v $(pwd)/{path}/:/DB_DOWNLOAD \\\n\
        {DICT_CONTAINERS['ncbi-datasets']} \\\n\
        sh -c "\\\n\
            datasets download genome accession GCF_000001405.40 --include gff3;\\\n\
            unzip ncbi_dataset.zip;\\\n\
            cp ncbi_dataset/data/GCF_000001405.40/genomic.gff /DB_DOWNLOAD/{filename.replace('gtf', 'gff')}\\\n\
    "\n\n\n"""

    cmd_gff_to_gff = f""" docker run -it --rm -v $(pwd)/{path}/:/DB_DOWNLOAD \\\n\
        {DICT_CONTAINERS['gffread']} \\\n\
        sh -c "\\\n\
            gffread /DB_DOWNLOAD/{filename.replace('gtf', 'gff')} -T -o /DB_DOWNLOAD/{filename}\\\n\
    "\n\n\n"""


    return title + cmd + cmd_gff_to_gff



def download_bed_GRCh38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=1.5 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.
    
    title = """echo "### Downloading GRCh38 gene bed file."\n\n"""
    os.makedirs(path, exist_ok=True)

    # If the gtf is NOT going to be downloaded we are going to perform the download of the gtf first, and then do the transformation to bed

    if 'Downloading GRCh38 gene gtf file' not in file_text:
        download_gtf_GRCh38(path_db.replace('.bed', '.gtf'), file_text)

    path, filename = os.path.split(path_db, file_text)

    cmd_gtf_to_bed = f""" docker run -it --rm -v $(pwd)/{path}/:/DB_DOWNLOAD \\\n\
        {DICT_CONTAINERS['bedops']} \\\n\
        sh -c "\\\n\
            gff2bed < /DB_DOWNLOAD/{filename.replace('bed', 'gff')} > /DB_DOWNLOAD/{filename}\\\n\
    "\n\n\n"""

    return title + cmd_gtf_to_bed



def download_star_index_GRCh38(path_db, file_text):
    subfiles_check = ['chrLength.txt', 'chrNameLength.txt' 'chrName.txt'  'chrStart.txt'  'exonGeTrInfo.tab'  'exonInfo.tab'  'geneInfo.tab'  
                      'Genome'  'genomeParameters.txt'  'Log.out'  'SA'  'SAindex'  'sjdbInfo.txt'  'sjdbList.fromGTF.out.tab'  
                      'sjdbList.out.tab'  'transcriptInfo.tab']
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=20 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    os.makedirs(path_db, exist_ok=True)

    title = """echo "### Generating GRCh38 STAR index file."\n\n"""

    # Since we will need the genome and gtf file to create the index, we will check if those files have been
    # scripted for download

    if 'Downloading GRCh38 transcript fasta file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['fasta_GRCh38'][0], file_text)

    if 'Downloading GRCh38 gene gtf file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['gtf_GRCh38'][0], file_text)


    
    # The number of bases in genomeSAindexNbases 14 is determined by the command 
    # gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${fasta}.fai
    # from https://github.com/nf-core/modules/blob/master/modules/nf-core/star/genomegenerate/main.nf
    # For GRCh38, it gives 14
    star_index_dir = path_db.strip('database/')

    cmd_STAR_build = f""" docker run -it --rm -v $(pwd)/database:/DB_DOWNLOAD \\\n\
        {DICT_CONTAINERS['STAR']} \\\n\
        sh -c " \\\n\
            STAR  \\\n\
            --runMode genomeGenerate \\\n\
            --genomeDir /DB_DOWNLOAD/{star_index_dir} \\\n\
            --genomeFastaFiles {DICT_DBS['fasta_GRCh38']}  \\\n\
            --sjdbGTFfile {DICT_DBS['gtf_GRCh38']} \\\n\
            --runThreadN {MAX_CPU}  \\\n\
            --genomeSAindexNbases 14  \\\n\
            --limitGenomeGenerateRAM {MAX_RAM}  \\\n\
    "\n\n\n"""
   
    return title + cmd_STAR_build



def download_bowtie_index_GRCh38(path_db, file_text):
    subfiles_check = ['genome.1.ebwt', 'genome.2.ebwt', 'genome.3.ebwt', 'genome.4.ebwt', 'genome.rev.1.ebwt', 'genome.rev.2.ebwt',]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=3 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    os.makedirs(path_db, exist_ok=True)
    os.makedirs('aws_downloads', exist_ok=True)

    cmd = f""" docker run -it --rm -v -v ~/.aws:/root/.aws -v $(pwd)/aws_downloads:/aws \\\n\
        {DICT_CONTAINERS['aws-cli']} \\\n\
        --no-sign-request --region eu-west-1 sync \\\n\
            s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/BowtieIndex /aws/BOWTIE
    "\n\n"""

    cmd_copy_rm = f"""mv $(pwd)/aws_downloads/BOWTIE $(pwd)/{path_db}\n\
        rm {path_db}/genome.fa\n\n"""
    
    
    return cmd + cmd_copy_rm



def download_bowtie2_index_GRCh38(path_db, file_text):
    subfiles_check = ['genome.1.ebwt', 'genome.2.ebwt', 'genome.3.ebwt', 'genome.4.ebwt', 'genome.rev.1.ebwt', 'genome.rev.2.ebwt',]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=3 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    os.makedirs(path_db, exist_ok=True)
    os.makedirs('aws_downloads', exist_ok=True)

    cmd = f""" docker run -it --rm -v -v ~/.aws:/root/.aws -v $(pwd)/aws_downloads:/aws \\\n\
        {DICT_CONTAINERS['aws-cli']} \\\n\
        --no-sign-request --region eu-west-1 sync \\\n\
            s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/BowtieIndex /aws/BOWTIE
    "\n\n"""

    cmd_copy_rm = f"""mv $(pwd)/aws_downloads/BOWTIE $(pwd)/{path_db}\n\
        rm {path_db}/genome.fa\n\n"""
    
    
    return cmd + cmd_copy_rm



def download_bwa_index_GRCh38(path_db, file_text):
    return ""

def download_hisat_index_GRCh38(path_db, file_text):
    return ""

def download_hisat2_index_GRCh38(path_db, file_text):
    return ""

def download_salmon_index_GRCh38(path_db, file_text):
    return ""

def download_txp2gene_GRCh38(path_db, file_text):
    return ""

def download_rsem_index_GRCh38(path_db, file_text):
    return ""

def download_kallisto_index_GRCh38(path_db, file_text):
    return ""

def download_kallisto_gene_map_GRCh38(path_db, file_text):
    return ""

def download_cellranger_index_GRCh38(path_db, file_text):
    return ""

def download_universc_index_GRCh38(path_db, file_text):
    return ""

def download_segemehl_index_GRCh38(path_db, file_text):
    return ""

def download_aa_data_repo_index_GRCh38(path_db, file_text):
    return ""




def download_fasta_CHM13(path_db, file_text):
    return ""

def download_bowtie2_index_CHM13(path_db, file_text):
    return ""







def download_fasta_GRCm38(path_db, file_text):
    return ""

def download_transcript_fasta_GRCm38(path_db, file_text):
    return ""

def download_gtf_GRCm38(path_db, file_text):
    return ""

def download_bed_GRCm38(path_db, file_text):
    return ""

def download_star_index_GRCm38(path_db, file_text):
    return ""

def download_bowtie_index_GRCm38(path_db, file_text):
    return ""

def download_bowtie2_index_GRCm38(path_db, file_text):
    return ""

def download_bwa_index_GRCm38(path_db, file_text):
    return ""

def download_hisat_index_GRCm38(path_db, file_text):
    return ""

def download_hisat2_index_GRCm38(path_db, file_text):
    return ""

def download_salmon_index_GRCm38(path_db, file_text):
    return ""

def download_txp2gene_GRCm38(path_db, file_text):
    return ""
    
def download_rsem_index_GRCm38(path_db, file_text):
    return ""

def download_kallisto_index_GRCm38(path_db, file_text):
    return ""

def download_kallisto_gene_map_GRCm38(path_db, file_text):
    return ""

def download_cellranger_index_GRCm38(path_db, file_text):
    return ""

def download_universc_index_GRCm38(path_db, file_text):
    return ""

def download_segemehl_index_GRCm38(path_db, file_text):
    return ""

def download_aa_data_repo_index_GRCm38(path_db, file_text):
    return ""




def download_mirbase_mature(path_db, file_text):
    return ""

def download_mirbase_hairpin(path_db, file_text):
    return ""

def download_mirbase_gtf_GRCh38(path_db, file_text):
    return ""

def download_mirbase_gtf_GRCm38(path_db, file_text):
    return ""


def download_mirgenedb_mature(path_db, file_text):
    return ""

def download_mirgenedb_hairpin(path_db, file_text):
    return ""

def download_mirgenedb_gff_GRCh38(path_db, file_text):
    return ""

def download_mirgenedb_gff_GRCm38(path_db, file_text):
    return ""





def download_kaiju(path_db, file_text):
    return ""

def download_centrifuge(path_db, file_text):
    return ""

def download_kraken2(path_db, file_text):
    return ""

def download_krakenuniq(path_db, file_text):
    return ""

def download_taxpasta(path_db, file_text):
    return ""


################################################################################################
# TERCIARY FUNCTIONS
################################################################################################
def db_check(path_check, subfiles_check=[], min_weight=0):
    """
    This function checks the databases to download. First, it checks if the path. 
        If the path is a file, it checks wheter the file exists. 
        If the path is a directory, it checks whether all files mentioned in subfiles_check exist.
        Lastly, it checks the weight of the file/directory. If it weights less than expected, it will download 
            the database because it is likely that it was not correctly downloaded.     

    The return of the function is True if the database has to be downloaded.
    """

    file_check = (os.path.isfile(path_check)) | (os.path.isdir(path_check))
    
    if not file_check:
        logger.info(f">>> Path {path_check} does not exist. Database will be downloaded.")
        return True
    else:
        if len(subfiles_check) > 0:
            for file in subfiles_check:
                if not os.path.isfile(f"{path_check}/{file}"):
                    logger.info(f">>> File {file} within path {path_check} does not exist. Database will be downloaded.")
                    return True

    # Check dir/file size
    if len(subfiles_check) > 0: # is a directory
        weight = sum(os.path.getsize(f) for f in os.listdir(path_check) if os.path.isfile(f))
    else:
        weight = os.path.getsize(path_check)
    
    if weight < min_weight:
        logger.info(f">>> Database size is lower than expected. Database will be downloaded.")
        return True
    
    return False                


def write_database_file(text, project):
    with open(f'work/{project}/database_download.sh', 'w') as file_out:
        file_out.write(text)


################################################################################################
# VARIABLES
################################################################################################
DICT_DBS = {'fasta_GRCh38': ('database/genomes/GRCh38/genome.fasta', download_fasta_GRCh38), 
            'transcript_fasta_GRCh38': ('database/genomes/GRCh38/transcript.fasta', download_transcript_fasta_GRCh38), 
            'gtf_GRCh38': ('database/genomes/GRCh38/genes.gtf', download_gtf_GRCh38), 
            'bed_GRCh38': ('database/genomes/GRCh38/genes.bed', download_bed_GRCh38), 
            'star_index_GRCh38': ('database/indexes/GRCh38/STAR', download_star_index_GRCh38),
            'bowtie_index_GRCh38': ('database/indexes/GRCh38/BOWTIE', download_bowtie_index_GRCh38),
            'bowtie2_index_GRCh38': ('database/indexes/GRCh38/BOWTIE2', download_bowtie2_index_GRCh38), 
            'bwa_index_GRCh38': ('database/indexes/GRCh38/BWA', download_bwa_index_GRCh38), 
            'hisat_index_GRCh38': ('database/indexes/GRCh38/HISAT', download_hisat_index_GRCh38), 
            'hisat2_index_GRCh38': ('database/indexes/GRCh38/HISAT2', download_hisat2_index_GRCh38), 
            'salmon_index_GRCh38': ('database/indexes/GRCh38/SALMON', download_salmon_index_GRCh38), 
            'txp2gene_GRCh38': ('database/indexes/GRCh38/SALMON/txp2gene', download_txp2gene_GRCh38),
            'rsem_index_GRCh38': ('database/indexes/GRCh38/rsem', download_rsem_index_GRCh38), 
            'kallisto_index_GRCh38': ('database/indexes/GRCh38/KALLISTO', download_kallisto_index_GRCh38), 
            'kallisto_gene_map_GRCh38': ('database/indexes/GRCh38/KALLISTO/gene_map', download_kallisto_gene_map_GRCh38), 
            'cellranger_index_GRCh38': ('database/indexes/GRCh38/cellranger', download_cellranger_index_GRCh38),
            'universc_index_GRCh38': ('database/indexes/GRCh38/UNIVERSC', download_universc_index_GRCh38), 
            'segemehl_index_GRCh38': ('database/indexes/GRCh38/segemehl', download_segemehl_index_GRCh38), 
            'aa_data_repo_index_GRCh38': ('database/indexes/GRCh38/aa_data_repo', download_aa_data_repo_index_GRCh38),
            'fasta_CHM13': ('database/genomes/CHM13/genome.fasta', download_fasta_CHM13),
            'bowtie2_index_CHM13': ('database/indexes/CHM13/BOWTIE2', download_bowtie2_index_CHM13),
            'fasta_GRCm38': ('database/genomes/GRCm38/genome.fasta', download_fasta_GRCm38), 
            'transcript_fasta_GRCm38': ('database/genomes/GRCm38/transcript.fasta', download_transcript_fasta_GRCm38),
            'gtf_GRCm38': ('database/genomes/GRCm38/genes.gtf', download_gtf_GRCm38), 
            'bed_GRCm38': ('database/genomes/GRCm38/genes.bed', download_bed_GRCm38), 
            'star_index_GRCm38': ('database/indexes/GRCm38/STAR', download_star_index_GRCm38),
            'bowtie_index_GRCm38': ('database/indexes/GRCm38/BOWTIE', download_bowtie_index_GRCm38),
            'bowtie2_index_GRCm38': ('database/indexes/GRCm38/BOWTIE2', download_bowtie2_index_GRCm38), 
            'bwa_index_GRCm38': ('database/indexes/GRCm38/BWA', download_bwa_index_GRCm38), 
            'hisat_index_GRCm38': ('database/indexes/GRCm38/HISAT', download_hisat_index_GRCm38), 
            'hisat2_index_GRCm38': ('database/indexes/GRCm38/HISAT2', download_hisat2_index_GRCm38), 
            'salmon_index_GRCm38': ('database/indexes/GRCm38/SALMON', download_salmon_index_GRCm38), 
            'txp2gene_GRCm38': ('database/indexes/GRCm38/SALMON/txp2gene', download_txp2gene_GRCm38),
            'rsem_index_GRCm38': ('database/indexes/GRCm38/rsem', download_rsem_index_GRCm38), 
            'kallisto_index_GRCm38': ('database/indexes/GRCm38/KALLISTO', download_kallisto_index_GRCm38),
            'kallisto_gene_map_GRCm38': ('database/indexes/GRCm38/KALLISTO/gene_map', download_kallisto_gene_map_GRCm38),  
            'cellranger_index_GRCm38': ('database/indexes/GRCm38/cellranger', download_cellranger_index_GRCm38), 
            'universc_index_GRCm38': ('database/indexes/GRCm38/UNIVERSC', download_universc_index_GRCm38), 
            'segemehl_index_GRCm38': ('database/indexes/GRCm38/segemehl', download_segemehl_index_GRCm38), 
            'aa_data_repo_index_GRCm38': ('database/indexes/GRCm38/aa_data_repo', download_aa_data_repo_index_GRCm38), 
            'mirbase_gtf_GRCh38': ('database/smRNA/mirbase/GRCh38.gtf', download_mirbase_gtf_GRCh38),
            'mirbase_gtf_GRCm38': ('database/smRNA/mirbase/GRCm38.gtf', download_mirbase_gtf_GRCm38),
            'mirgenedb_gff_GRCh38': ('database/smRNA/mirgenedb/GRCh38.gff', download_mirgenedb_gff_GRCh38),
            'mirgenedb_gff_GRCm38': ('database/smRNA/mirgenedb/GRCm38.gff', download_mirgenedb_gff_GRCm38),
            'mirbase_mature': ('database/smRNA/mirbase/mature.gff', download_mirbase_mature),
            'mirbase_hairpin': ('database/smRNA/mirbase/hairpin.gff', download_mirbase_hairpin),
            'mirgenedb_mature': ('database/smRNA/mirgenedb/mature.gff', download_mirgenedb_mature),
            'mirgenedb_hairpin': ('database/smRNA/mirgenedb/hairpin.gff', download_mirgenedb_hairpin),
            'mature': ('database/mirbase/mature.fa', download_mirbase_mature),
            'hairpin': ('database/mirbase/hairpin.fa', download_mirbase_hairpin),
            'mirna_gtf_GRCh38': ('database/GRCh38.gtf', download_mirbase_gtf_GRCh38),
            'mirna_gtf_GRCm38': ('database/GRCm38.gtf', download_mirbase_gtf_GRCm38), 
            'kaiju': ('database/taxprofiler/kaiju', download_kaiju), 
            'centrifuge': ('database/taxprofiler/centrifuge', download_centrifuge), 
            'kraken2': ('database/taxprofiler/kraken2', download_kraken2), 
            'krakenuniq': ('database/taxprofiler/krakenuniq', download_krakenuniq), 
            'taxpasta': ('database/taxprofiler/taxpasta', download_taxpasta)}
