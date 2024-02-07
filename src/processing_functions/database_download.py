import logging
import os

from .versions import VERSION
from .global_vars_and_funcs import MAX_CPU, MAX_RAM, MB, GB, DICT_CONTAINERS, ENSEMBL_VERSION_HUMAN, ENSEMBL_VERSION_MOUSE, HISAT2_LIMIT_MEM



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







################################################################################################
# SECONDARY FUNCTIONS
################################################################################################
def download_fasta_GRCh38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=3 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

        
    path, filename = os.path.split(path_db)

    title = """echo "### Downloading GRCh38 genome fasta file."\n\n"""

    cmd = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['aria2c']} \\\n\
        --file-allocation=none -x {MAX_CPU} -d /{path} https://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION_HUMAN}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz\n\
            gzip -d {path}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz\n\
            mv {path}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa {path}/{filename}
    \n\n\n"""


    cmd_index_fai = f""" docker run -it --rm -v $(pwd)/{path}/:/database \\\n\
        {DICT_CONTAINERS['samtools']} \\\n\
        sh -c "\\\n\
            samtools faidx /database/{filename}\\\n\
    "\n\n\n"""  
    
    
    return title + cmd + cmd_index_fai

   

def download_transcript_fasta_GRCh38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=440 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

        
    path, filename = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)
    
    title = """echo "### Downloading GRCh38 transcript fasta file."\n\n"""

    cmd = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['aria2c']} \\\n\
        --file-allocation=none -x {MAX_CPU} -d /{path} https://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION_HUMAN}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz\n\
            gzip -d {path}/Homo_sapiens.GRCh38.cdna.all.fa.gz\n\
            mv {path}/Homo_sapiens.GRCh38.cdna.all.fa {path}/{filename}
    \n\n\n"""
    



    return title + cmd



def download_gtf_GRCh38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=500 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

        
    path, filename = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    title = """echo "### Downloading GRCh38 gene gtf file."\n\n"""

    cmd_gtf = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['aria2c']} \\\n\
        --file-allocation=none -x {MAX_CPU} -d /{path} https://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION_HUMAN}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{ENSEMBL_VERSION_HUMAN}.gtf.gz\n\
            gzip -d {path}/Homo_sapiens.GRCh38.{ENSEMBL_VERSION_HUMAN}.gtf.gz\n\
            mv {path}/Homo_sapiens.GRCh38.{ENSEMBL_VERSION_HUMAN}.gtf {path}/{filename}
    \n\n\n"""


    cmd_gff = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['aria2c']} \\\n\
        --file-allocation=none -x {MAX_CPU} -d /{path} https://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION_HUMAN}/gff3/homo_sapiens/Homo_sapiens.GRCh38.{ENSEMBL_VERSION_HUMAN}.gff3.gz\n\
            gzip -d {path}/Homo_sapiens.GRCh38.{ENSEMBL_VERSION_HUMAN}.gff3.gz\n\
            mv {path}/Homo_sapiens.GRCh38.{ENSEMBL_VERSION_HUMAN}.gff3 {path}/{filename.replace('.gtf', '.gff')}
    \n\n\n"""


    return title + cmd_gtf + cmd_gff



def download_gtf_corrected_GRCh38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=500 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

        
    path, filename = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    if 'Downloading GRCh38 gene gtf file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['gtf_GRCh38'][0], file_text)

    title = """echo "### Applying correction to GRCh38 gtf file."\n\n"""

    cmd_gtf = f""" python src/others/generate_correct_gtf.py \\\n\
        --input {path_db.replace('_corrected', '')}  \\\n\
        --output {path_db}
    \n\n\n"""


    return title + cmd_gtf



def download_bed_GRCh38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=500 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.
    
    title = """echo "### Downloading GRCh38 gene bed file."\n\n"""
    path, filename = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    # If the gtf is NOT going to be downloaded we are going to perform the download of the gtf first, and then do the transformation to bed

    if 'Downloading GRCh38 gene gtf file' not in file_text:
        download_gtf_GRCh38(path_db.replace('.bed', '.gtf'), file_text)

    cmd_curl = f""" docker run -it --rm -v $(pwd)/{path}/:/database \\\n\
        {DICT_CONTAINERS['curl']} \\\n\
        sh -c "\\\n\
           curl https://raw.githubusercontent.com/nf-core/rnaseq/3.14.0/bin/gtf2bed > /database/gtf2bed \n
    "\n\n\n"""

    cmd_gtf_to_bed = f""" docker run -it --rm -v $(pwd)/{path}/:/database \\\n\
        {DICT_CONTAINERS['perl']} \\\n\
        sh -c "\\\n\
           perl /database/gtf2bed /database/{filename.replace('bed', 'gtf')} > /database/{filename}\n\
           rm /database/gtf2bed
    "\n\n\n"""

    return title + cmd_curl + cmd_gtf_to_bed


# BUG: FALLA HACER EL GENOME GENERATE!!!!
def download_star_index_GRCh38(path_db, file_text): 
    subfiles_check = ['chrLength.txt', 'chrNameLength.txt', 'chrName.txt', 'chrStart.txt', 'exonGeTrInfo.tab', 'exonInfo.tab', 'geneInfo.tab', 
                      'Genome', 'genomeParameters.txt', 'SA', 'SAindex', 'sjdbInfo.txt', 'sjdbList.fromGTF.out.tab', 'sjdbList.out.tab', 'transcriptInfo.tab']
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=20 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    os.makedirs(path_db, exist_ok=True)

    title = """echo "### Generating GRCh38 STAR index file."\n\n"""

    # Since we will need the genome and gtf file to create the index, we will check if those files have been
    # scripted for download

    if 'Downloading GRCh38 genome fasta file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['genome_fasta_GRCh38'][0], file_text)

    if 'Downloading GRCh38 gene gtf file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['gtf_GRCh38'][0], file_text)


    
    # The number of bases in genomeSAindexNbases 14 is determined by the command 
    # gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${fasta}.fai
    # from https://github.com/nf-core/modules/blob/master/modules/nf-core/star/genomegenerate/main.nf
    # For GRCh38, it gives 14

    cmd_STAR_build = f""" docker run -it --rm -v $(pwd)/database:/database \\\n\
        --cpu-shares 16384 --memory {MAX_RAM}g \\\n\
        {DICT_CONTAINERS['STAR']} \\\n\
        sh -c " \\\n\
            STAR  \\\n\
            --runMode genomeGenerate \\\n\
            --genomeDir /{path_db} \\\n\
            --genomeFastaFiles /{DICT_DBS['genome_fasta_GRCh38'][0]}  \\\n\
            --sjdbGTFfile /{DICT_DBS['gtf_GRCh38'][0]} \\\n\
            --runThreadN {MAX_CPU}  \\\n\
            --genomeSAindexNbases 14  \\\n\
            --limitGenomeGenerateRAM {MAX_RAM * 1000 ** 3}  \n\
    "\n\n\n"""
   
    return title + cmd_STAR_build




def download_bowtie_index_GRCh38(path_db, file_text):
    subfiles_check = ['genome.1.ebwt', 'genome.2.ebwt', 'genome.3.ebwt', 'genome.4.ebwt', 'genome.rev.1.ebwt', 'genome.rev.2.ebwt',]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=3 * GB): 
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    os.makedirs(path_db, exist_ok=True)

    title = """echo "### Generating GRCh38 BOWTIE index file."\n\n"""

    # Since we will need the genome and gtf file to create the index, we will check if those files have been
    # scripted for download

    if 'Downloading GRCh38 genome fasta file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['genome_fasta_GRCh38'][0], file_text)


    cmd_BOWTIE_build = f""" docker run -it --rm -v $(pwd)/database:/database \\\n\
        {DICT_CONTAINERS['bowtie']} \\\n\
        sh -c " \\\n\
            bowtie-build \\\n\
                --threads {MAX_CPU} \\\n\
                    /{DICT_DBS['genome_fasta_GRCh38'][0]} \\\n\
                    /{DICT_DBS['bowtie_index_GRCh38'][0]}/genome
    "\n\n\n"""

    
    return title + cmd_BOWTIE_build



def download_bowtie2_index_GRCh38(path_db, file_text):
    subfiles_check = ['genome.1.bt2', 'genome.2.bt2', 'genome.3.bt2', 'genome.4.bt2', 'genome.rev.1.bt2', 'genome.rev.2.bt2',]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=4 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    os.makedirs(path_db, exist_ok=True)

    title = """echo "### Generating GRCh38 BOWTIE2 index file."\n\n"""

    # Since we will need the genome and gtf file to create the index, we will check if those files have been
    # scripted for download

    if 'Downloading GRCh38 genome fasta file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['genome_fasta_GRCh38'][0], file_text)


    cmd_BOWTIE2_build = f""" docker run -it --rm -v $(pwd)/database:/database \\\n\
        {DICT_CONTAINERS['bowtie2']} \\\n\
        sh -c " \\\n\
            bowtie2-build \\\n\
                --threads {MAX_CPU} \\\n\
                    /{DICT_DBS['genome_fasta_GRCh38'][0]} \\\n\
                    /{DICT_DBS['bowtie2_index_GRCh38'][0]}/genome
    "\n\n\n"""

    return title + cmd_BOWTIE2_build



def download_bwa_index_GRCh38(path_db, file_text):
    subfiles_check = [f'genome.{termination}' for termination in ['amb', 'ann', 'bwt', 'pac', 'sa']]

    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=5 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    os.makedirs(path_db, exist_ok=True)

    title = """echo "### Generating GRCh38 BWA index file."\n\n"""

    # Since we will need the genome and gtf file to create the index, we will check if those files have been
    # scripted for download

    if 'Downloading GRCh38 genome fasta file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['genome_fasta_GRCh38'][0], file_text)

    cmd_BWA_build = f""" docker run -it --rm -v $(pwd)/database:/database \\\n\
        {DICT_CONTAINERS['bwa']} \\\n\
        sh -c " \\\n\
            bwa \\\n\
                index \\\n\
                -p /{DICT_DBS['bwa_index_GRCh38'][0]}/genome \\\n\
                /{DICT_DBS['genome_fasta_GRCh38'][0]}
    "\n\n\n"""

    return title + cmd_BWA_build


def download_hisat2_index_GRCh38(path_db, file_text):
    subfiles_check = [f"genome.{N}.ht2" for N in range(1,9)]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=4 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCh38 hisat2 index file."\n\n"""

    os.makedirs(path_db, exist_ok=True)

    if 'Downloading GRCh38 genome fasta file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['genome_fasta_GRCh38'][0], file_text)


    if MAX_RAM <= HISAT2_LIMIT_MEM:
        cmd_HISAT2_index_GRCh38 = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['hisat2']} \\\n\
        sh -c "\\\n\
            hisat2-build  \\\n\
            -p {MAX_CPU} \\\n\
            /{DICT_DBS['genome_fasta_GRCh38'][0]}  \\\n\
            /{DICT_DBS['hisat2_index_GRCh38'][0]}/genome
        "\n\n\n"""
    else:

        if 'Downloading GRCh38 gene gtf file' not in file_text:
                download_gtf_GRCh38(path_db.replace('.bed', '.gtf'), file_text)

        cmd_HISAT2_index_GRCh38 = f"""docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['hisat2']} \\\n\
        sh -c "\\\n\
            hisat2_extract_splice_sites.py /database/genomes/GRCh38/genes.gtf > genome.splicesites.txt \n\
            hisat2_extract_exons.py /database/genomes/GRCh38/genes.gtf > genome.exons.txt \n\
            hisat2-build  \\\n\
            -p {MAX_CPU} \\\n\
            --ss genome.splicesites.txt \\\n\
            --exon genome.exons.txt \\\n\
            /{DICT_DBS['genome_fasta_GRCh38'][0]}  \\\n\
            /{DICT_DBS['hisat2_index_GRCh38'][0]}/genome
        "\n\n\n"""


    return title + cmd_HISAT2_index_GRCh38


def download_salmon_index_GRCh38(path_db, file_text):
    subfiles_check = ['complete_ref_lens.bin', 'ctable.bin', 'ctg_offsets.bin', 'duplicate_clusters.tsv', 'info.json', 'mphf.bin', 'pos.bin', 
                      'pre_indexing.log', 'rank.bin', 'refAccumLengths.bin', 'ref_indexing.log', 'reflengths.bin', 'refseq.bin',
                       'seq.bin']
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=700 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCh38 salmon index file."\n\n"""
    path, filename = os.path.split(path_db)

    os.makedirs(path, exist_ok=True)

    if 'Downloading GRCh38 transcript fasta' not in file_text:
        download_gtf_GRCh38(DICT_DBS['transcript_fasta_GRCh38'][0], file_text)

    cmd_salmon_index_GRCh38  = f"""docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['salmon']} \\\n\
        sh -c "\\\n\
            salmon index \\\n\
            -p {MAX_CPU} \\\n\
            -t /{DICT_DBS['transcript_fasta_GRCh38'][0]} \\\n\
            -i /{DICT_DBS['salmon_index_GRCh38'][0]}
        "\n\n\n"""

    return title + cmd_salmon_index_GRCh38



def download_txp2gene_GRCh38(path_db, file_text):
    subfiles_check = []

    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=7 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCh38 salmon txp2gene file."\n\n"""

    if 'Downloading GRCh38 gene gtf file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['gtf_GRCh38'][0], file_text)

    cmd_salmon_txp2gene_GRCh38  = f"""docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['gffread']} \\\n\
        sh -c "\\\n\
            gffread \\\n\
                /{DICT_DBS['gtf_GRCh38'][0]} \\\n\
                --table transcript_id,gene_id \\\n\
                    > /{DICT_DBS['txp2gene_GRCh38'][0]}
        "\n\n\n"""

    return title + cmd_salmon_txp2gene_GRCh38



def download_rsem_index_GRCh38(path_db, file_text):
    subfiles_check = [f'genome.{suffix}' for suffix in ['chrlist', 'grp', 'idx.fa', 'n2g.idx.fa', 'seq', 'ti', 'transcripts.fa']] + \
                     ['chrLength.txt', 'chrNameLength.txt', 'chrName.txt',  'chrStart.txt',  'exonGeTrInfo.tab',  'exonInfo.tab',  'geneInfo.tab',  
                      'Genome',  'genomeParameters.txt',  'Log.out',  'SA',  'SAindex',  'sjdbInfo.txt',  'sjdbList.fromGTF.out.tab',  
                      'sjdbList.out.tab',  'transcriptInfo.tab']
    
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=30 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCh38 rsem index file."\n\n"""
    os.makedirs(path_db, exist_ok=True)

    if 'Downloading GRCh38 genome fasta file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['genome_fasta_GRCh38'][0], file_text)

    if 'Downloading GRCh38 gene gtf file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['gtf_GRCh38'][0], file_text)


    cmd_rsem_index_GRCh38  = f"""docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['rsem-star']} \\\n\
        sh -c "\\\n\
            rsem-prepare-reference \\\n\
            --gtf \{DICT_DBS['gtf_GRCh38'][0]} \\\n\
            --star \\\n\
            --star-path /usr/local/bin \\\n\
            -p {MAX_CPU} \\\n\
            /{DICT_DBS['genome_fasta_GRCh38'][0]} \\\n\
            /{DICT_DBS['rsem_index_GRCh38'][0]}/genome
        "\n\n\n"""

    
    return title + cmd_rsem_index_GRCh38



def download_kallisto_index_GRCh38(path_db, file_text):
    subfiles_check = []
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=300 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCh38 kallisto index file."\n\n"""
    path, filename = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    if 'Downloading GRCh38 gene gtf file' not in file_text:
            download_gtf_GRCh38(path_db.replace('.bed', '.gtf'), file_text)

    cmd_kallisto_index_GRCh38 = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['kallisto']} \\\n\
        sh -c "\\\n\
            kallisto index -i {filename} \\\n\
            /{DICT_DBS['transcript_fasta_GRCh38'][0]} && mv {filename} /{path}/{filename}
    "\n\n\n"""


    return title + cmd_kallisto_index_GRCh38


def download_cellranger_index_GRCh38(path_db, file_text):
    subfiles_check = ['fasta/genome.fa', 'fasta/genome.fa.fai', 'genes/genes.gtf.gz', 'reference.json'] + \
                     [f'star/{file}' for file in ['chrLength.txt', 'chrNameLength.txt', 'chrName.txt', 'chrStart.txt', 'exonGeTrInfo.tab', 
                                                  'exonInfo.tab', 'geneInfo.tab', 'Genome', 'genomeParameters.txt', 'SA', 'SAindex',  'sjdbInfo.txt', 
                                                  'sjdbList.fromGTF.out.tab', 'sjdbList.out.tab', 'transcriptInfo.tab']]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=16 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCh38 cellranger index file."\n\n"""
    path, filename = os.path.split(path_db)
    os.makedirs(path_db, exist_ok=True)

    if 'Downloading GRCh38 genome fasta file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['genome_fasta_GRCh38'][0], file_text)
    
    if 'Downloading GRCh38 gene gtf file' not in file_text:
            download_gtf_GRCh38(path_db.replace('.bed', '.gtf'), file_text)

    cmd_cellranger_index_GRCh38 = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['cellranger']} \\\n\
        sh -c "\\\n\
            cellranger mkref \\\n\
            --genome  {filename} \\\n\
            --fasta /{DICT_DBS['genome_fasta_GRCh38'][0]} \\\n\
            --genes /{DICT_DBS['gtf_GRCh38'][0]} \\\n\
            --nthreads={MAX_CPU}\n\
            mv /cellranger /{path}
    "\n\n\n"""

    return title + cmd_cellranger_index_GRCh38



def download_universc_index_GRCh38(path_db, file_text):
    subfiles_check = ['fasta/genome.fa', 'fasta/genome.fa.fai', 'genes/genes.gtf.gz', 'reference.json'] + \
                     [f'star/{file}' for file in ['chrLength.txt', 'chrNameLength.txt', 'chrName.txt', 'chrStart.txt', 'exonGeTrInfo.tab', 
                                                  'exonInfo.tab', 'geneInfo.tab', 'Genome', 'genomeParameters.txt', 'SA', 'SAindex',  'sjdbInfo.txt', 
                                                  'sjdbList.fromGTF.out.tab', 'sjdbList.out.tab', 'transcriptInfo.tab']]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=16 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCh38 universc index file."\n\n"""
    path, filename = os.path.split(path_db)
    os.makedirs(path_db, exist_ok=True)

    if 'Downloading GRCh38 genome fasta file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['genome_fasta_GRCh38'][0], file_text)
    
    if 'Downloading GRCh38 gene gtf file' not in file_text:
            download_gtf_GRCh38(path_db.replace('.bed', '.gtf'), file_text)

    cmd_universc_index_GRCh38 = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['cellranger']} \\\n\
        sh -c "\\\n\
            cellranger mkref \\\n\
            --genome  {filename} \\\n\
            --fasta /{DICT_DBS['genome_fasta_GRCh38'][0]} \\\n\
            --genes /{DICT_DBS['gtf_GRCh38'][0]} \\\n\
            --nthreads={MAX_CPU}\n\
            mv /universc /{path}
    "\n\n\n"""

    return title + cmd_universc_index_GRCh38


def download_segemehl_index_GRCh38(path_db, file_text):
    subfiles_check = ['genome.idx']
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=45 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCh38 segemehl index file."\n\n"""
    os.makedirs(path_db, exist_ok=True)

    if 'Downloading GRCh38 genome fasta file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['genome_fasta_GRCh38'][0], file_text)
    
    if 'Downloading GRCh38 gene gtf file' not in file_text:
            download_gtf_GRCh38(path_db.replace('.bed', '.gtf'), file_text)

    cmd_segemehl_index_GRCh38 = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['segemehl']} \\\n\
        sh -c "\\\n\
            segemehl.x \\\n\
            -t {MAX_CPU} \\\n\
            -d /{DICT_DBS['genome_fasta_GRCh38'][0]} \\\n\
            -x /{path_db}/genome.idx 
    "\n\n\n"""


    return title + cmd_segemehl_index_GRCh38



def download_aa_data_repo_index_GRCh38(path_db, file_text):
    subfiles_check = ['GRCh38/' + i for i in ['annotations/hg38GenomicSuperDup.tab', 'cancer/oncogene_list_hg38.txt', 'cancer/oncogene_list.txt', 
                      'dummy_ploidy.vcf', 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fa', 
                      'GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai', 'GRCh38_centromere.bed', 'hg38full_k35_noMM.mappability.bedgraph',
                      'exclude.cnvnator_100bp.GRCh38.20170403.bed', 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.amb', 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.pac', 'GRCh38_cnvkit_filtered_ref.cnn', 'last_updated.txt', 
                      'chrom_list.txt', 'file_list.txt', 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.ann', 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.sa', 'GRCh38_merged_centromeres_conserved_sorted.bed', 'refGene.txt', 
                      'conserved_gain5_hg38.bed', 'file_sources.txt', 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.bwt', 'Genes_hg38.gff', 'GRCh38_noAlt.fa.fai']]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=9 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCh38 ampliconarchitect index file."\n\n"""
    os.makedirs(path_db, exist_ok=True)

    cmd_aa_repo_index_GRCh38 = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['aria2c']} \\\n\
         -x {MAX_CPU} \\\n\
         -d /{path_db} \\\n\
         --file-allocation=none \\\n\
         https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/GRCh38_indexed.tar.gz \n\
         tar xzvf {path_db}/GRCh38_indexed.tar.gz -C {path_db}\n\
         rm -rf {path_db}/GRCh38_indexed.tar.gz 
    \n\n\n"""

    return title + cmd_aa_repo_index_GRCh38































def download_fasta_CHM13(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=3 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

        
    path, _ = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    title = """echo "### Downloading CHM13 genome fasta file."\n\n"""


    cmd = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['ncbi-datasets']} \\\n\
        datasets download genome \\\n\
        accession GCA_009914755.3 \\\n\
        --filename {path_db}.zip \n\
        unzip {path_db}.zip -d {path} \n\
        mv {path}/ncbi_dataset/data/GCA_009914755.3/GCA_009914755.3_T2T-CHM13v1.1_genomic.fna {path_db} \n\
        rm -rf {path}/README.md {path_db}.zip {path}/ncbi_dataset
    \n\n\n"""
    
    return title + cmd 


def download_bowtie2_index_CHM13(path_db, file_text):
    subfiles_check = ['genome.1.ebwt', 'genome.2.ebwt', 'genome.3.ebwt', 'genome.4.ebwt', 'genome.rev.1.ebwt', 'genome.rev.2.ebwt',]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=3 * GB): 
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    os.makedirs(path_db, exist_ok=True)

    title = """echo "### Generating CHM13 BOWTIE index file."\n\n"""

    # Since we will need the genome and gtf file to create the index, we will check if those files have been
    # scripted for download

    if 'Downloading CHM13 genome fasta file' not in file_text:
        download_gtf_GRCh38(DICT_DBS['genome_fasta_CHM13'][0], file_text)


    cmd_BOWTIE_build = f""" docker run -it --rm -v $(pwd)/database:/database \\\n\
        {DICT_CONTAINERS['bowtie']} \\\n\
        sh -c " \\\n\
            bowtie-build \\\n\
                --threads {MAX_CPU} \\\n\
                    /{DICT_DBS['genome_fasta_CHM13'][0]} \\\n\
                    /{DICT_DBS['bowtie2_index_CHM13'][0]}/genome
    "\n\n\n"""

    
    return title + cmd_BOWTIE_build





































def download_fasta_GRCm38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=2.7 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

        
    path, filename = os.path.split(path_db)

    title = """echo "### Downloading GRCm38 genome fasta file."\n\n"""

    cmd = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['aria2c']} \\\n\
        --file-allocation=none -x {MAX_CPU} -d /{path} https://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION_MOUSE}/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz\n\
            gzip -d {path}/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz\n\
            mv {path}/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa {path}/{filename}
    \n\n\n"""


    cmd_index_fai = f""" docker run -it --rm -v $(pwd)/{path}/:/database \\\n\
        {DICT_CONTAINERS['samtools']} \\\n\
        sh -c "\\\n\
            samtools faidx /database/{filename}\\\n\
    "\n\n\n"""  
    
    
    return title + cmd + cmd_index_fai

   

def download_transcript_fasta_GRCm38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=250 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

        
    path, filename = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)
    
    title = """echo "### Downloading GRCm38 transcript fasta file."\n\n"""

    cmd = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['aria2c']} \\\n\
        --file-allocation=none -x {MAX_CPU} -d /{path} https://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION_MOUSE}/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz\n\
            gzip -d {path}/Mus_musculus.GRCm38.cdna.all.fa.gz\n\
            mv {path}/Mus_musculus.GRCm38.cdna.all.fa {path}/{filename}
    \n\n\n"""


    return title + cmd



def download_gtf_GRCm38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=500 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

        
    path, filename = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    title = """echo "### Downloading GRCm38 gene gtf file."\n\n"""

    cmd_gtf = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['aria2c']} \\\n\
        --file-allocation=none -x {MAX_CPU} -d /{path} https://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION_MOUSE}/gtf/mus_musculus/Mus_musculus.GRCm38.{ENSEMBL_VERSION_MOUSE}.gtf.gz\n\
            gzip -d {path}/Mus_musculus.GRCm38.{ENSEMBL_VERSION_MOUSE}.gtf.gz\n\
            mv {path}/Mus_musculus.GRCm38.{ENSEMBL_VERSION_MOUSE}.gtf {path}/{filename}
    \n\n\n"""


    cmd_gff = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['aria2c']} \\\n\
        --file-allocation=none -x {MAX_CPU} -d /{path} https://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION_MOUSE}/gff3/mus_musculus/Mus_musculus.GRCm38.{ENSEMBL_VERSION_MOUSE}.gff3.gz\n\
            gzip -d {path}/Mus_musculus.GRCm38.{ENSEMBL_VERSION_MOUSE}.gff3.gz\n\
            mv {path}/Mus_musculus.GRCm38.{ENSEMBL_VERSION_MOUSE}.gff3 {path}/{filename.replace('.gtf', '.gff')}
    \n\n\n"""


    return title + cmd_gtf + cmd_gff



def download_gtf_corrected_GRCm38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=500 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

        
    path, filename = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    if 'Downloading GRCm38 gene gtf file' not in file_text:
        download_gtf_GRCm38(DICT_DBS['gtf_GRCm38'][0], file_text)

    title = """echo "### Applying correction to GRCm38 gtf file."\n\n"""

    cmd_gtf = f""" python src/others/generate_correct_gtf.py \\\n\
        --input {path_db.replace('_corrected', '')}  \\\n\
        --output {path_db}
    \n\n\n"""


    return title + cmd_gtf



def download_bed_GRCm38(path_db, file_text):
    if not db_check(path_db, subfiles_check=[], min_weight=500 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.
    
    title = """echo "### Downloading GRCm38 gene bed file."\n\n"""
    path, filename = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    # If the gtf is NOT going to be downloaded we are going to perform the download of the gtf first, and then do the transformation to bed

    if 'Downloading GRCm38 gene gtf file' not in file_text:
        download_gtf_GRCm38(path_db.replace('.bed', '.gtf'), file_text)


    cmd_curl = f""" docker run -it --rm -v $(pwd)/{path}/:/database \\\n\
        {DICT_CONTAINERS['curl']} \\\n\
        sh -c "\\\n\
           curl https://raw.githubusercontent.com/nf-core/rnaseq/3.14.0/bin/gtf2bed > /database/gtf2bed \n
    "\n\n\n"""

    cmd_gtf_to_bed = f""" docker run -it --rm -v $(pwd)/{path}/:/database \\\n\
        {DICT_CONTAINERS['perl']} \\\n\
        sh -c "\\\n\
           perl /database/gtf2bed /database/{filename.replace('bed', 'gtf')} > /database/{filename}\n\
           rm /database/gtf2bed
    "\n\n\n"""

    return title + cmd_curl + cmd_gtf_to_bed



def download_star_index_GRCm38(path_db, file_text): 
    subfiles_check = ['chrLength.txt', 'chrNameLength.txt', 'chrName.txt', 'chrStart.txt', 'exonGeTrInfo.tab', 'exonInfo.tab', 'geneInfo.tab', 
                      'Genome', 'genomeParameters.txt', 'SA', 'SAindex', 'sjdbInfo.txt', 'sjdbList.fromGTF.out.tab', 'sjdbList.out.tab', 'transcriptInfo.tab']
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=25 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    os.makedirs(path_db, exist_ok=True)

    title = """echo "### Generating GRCm38 STAR index file."\n\n"""

    # Since we will need the genome and gtf file to create the index, we will check if those files have been
    # scripted for download

    if 'Downloading GRCm38 genome fasta file' not in file_text:
        download_gtf_GRCm38(DICT_DBS['genome_fasta_GRCm38'][0], file_text)

    if 'Downloading GRCm38 gene gtf file' not in file_text:
        download_gtf_GRCm38(DICT_DBS['gtf_GRCm38'][0], file_text)


    
    # The number of bases in genomeSAindexNbases 14 is determined by the command 
    # gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${fasta}.fai
    # from https://github.com/nf-core/modules/blob/master/modules/nf-core/star/genomegenerate/main.nf
    # For GRCm38, it gives 14

    cmd_STAR_build = f""" docker run -it --rm -v $(pwd)/database:/database \\\n\
        --cpu-shares 16384 --memory {MAX_RAM}g \\\n\
        {DICT_CONTAINERS['STAR']} \\\n\
        sh -c " \\\n\
            STAR  \\\n\
            --runMode genomeGenerate \\\n\
            --genomeDir /{path_db} \\\n\
            --genomeFastaFiles /{DICT_DBS['genome_fasta_GRCm38'][0]}  \\\n\
            --sjdbGTFfile /{DICT_DBS['gtf_GRCm38'][0]} \\\n\
            --runThreadN {MAX_CPU}  \\\n\
            --genomeSAindexNbases 14  \\\n\
            --limitGenomeGenerateRAM {MAX_RAM * 1000 ** 3}  \n\
    "\n\n\n"""
   
    return title + cmd_STAR_build




def download_bowtie_index_GRCm38(path_db, file_text):
    subfiles_check = ['genome.1.ebwt', 'genome.2.ebwt', 'genome.3.ebwt', 'genome.4.ebwt', 'genome.rev.1.ebwt', 'genome.rev.2.ebwt',]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=2.5 * GB): 
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    os.makedirs(path_db, exist_ok=True)

    title = """echo "### Generating GRCm38 BOWTIE index file."\n\n"""

    # Since we will need the genome and gtf file to create the index, we will check if those files have been
    # scripted for download

    if 'Downloading GRCm38 genome fasta file' not in file_text:
        download_gtf_GRCm38(DICT_DBS['genome_fasta_GRCm38'][0], file_text)


    cmd_BOWTIE_build = f""" docker run -it --rm -v $(pwd)/database:/database \\\n\
        {DICT_CONTAINERS['bowtie']} \\\n\
        sh -c " \\\n\
            bowtie-build \\\n\
                --threads {MAX_CPU} \\\n\
                    /{DICT_DBS['genome_fasta_GRCm38'][0]} \\\n\
                    /{DICT_DBS['bowtie_index_GRCm38'][0]}/genome
    "\n\n\n"""

    
    return title + cmd_BOWTIE_build



def download_bowtie2_index_GRCm38(path_db, file_text):
    subfiles_check = ['genome.1.bt2', 'genome.2.bt2', 'genome.3.bt2', 'genome.4.bt2', 'genome.rev.1.bt2', 'genome.rev.2.bt2',]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=3.5 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    os.makedirs(path_db, exist_ok=True)

    title = """echo "### Generating GRCm38 BOWTIE2 index file."\n\n"""

    # Since we will need the genome and gtf file to create the index, we will check if those files have been
    # scripted for download

    if 'Downloading GRCm38 genome fasta file' not in file_text:
        download_gtf_GRCm38(DICT_DBS['genome_fasta_GRCm38'][0], file_text)


    cmd_BOWTIE2_build = f""" docker run -it --rm -v $(pwd)/database:/database \\\n\
        {DICT_CONTAINERS['bowtie2']} \\\n\
        sh -c " \\\n\
            bowtie2-build \\\n\
                --threads {MAX_CPU} \\\n\
                    /{DICT_DBS['genome_fasta_GRCm38'][0]} \\\n\
                    /{DICT_DBS['bowtie2_index_GRCm38'][0]}/genome
    "\n\n\n"""

    return title + cmd_BOWTIE2_build



def download_bwa_index_GRCm38(path_db, file_text):
    subfiles_check = [f'genome.{termination}' for termination in ['amb', 'ann', 'bwt', 'pac', 'sa']]

    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=4.5 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    os.makedirs(path_db, exist_ok=True)

    title = """echo "### Generating GRCm38 BWA index file."\n\n"""

    # Since we will need the genome and gtf file to create the index, we will check if those files have been
    # scripted for download

    if 'Downloading GRCm38 genome fasta file' not in file_text:
        download_gtf_GRCm38(DICT_DBS['genome_fasta_GRCm38'][0], file_text)

    cmd_BWA_build = f""" docker run -it --rm -v $(pwd)/database:/database \\\n\
        {DICT_CONTAINERS['bwa']} \\\n\
        sh -c " \\\n\
            bwa \\\n\
                index \\\n\
                -p /{DICT_DBS['bwa_index_GRCm38'][0]}/genome \\\n\
                /{DICT_DBS['genome_fasta_GRCm38'][0]}
    "\n\n\n"""

    return title + cmd_BWA_build


def download_hisat2_index_GRCm38(path_db, file_text):
    subfiles_check = [f"genome.{N}.ht2" for N in range(1,9)]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=4 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCm38 hisat2 index file."\n\n"""

    os.makedirs(path_db, exist_ok=True)

    if 'Downloading GRCm38 genome fasta file' not in file_text:
        download_gtf_GRCm38(DICT_DBS['genome_fasta_GRCm38'][0], file_text)


    if MAX_RAM <= HISAT2_LIMIT_MEM:
        cmd_HISAT2_index_GRCm38 = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['hisat2']} \\\n\
        sh -c "\\\n\
            hisat2-build  \\\n\
            -p {MAX_CPU} \\\n\
            /{DICT_DBS['genome_fasta_GRCm38'][0]}  \\\n\
            /{DICT_DBS['hisat2_index_GRCm38'][0]}/genome
        "\n\n\n"""
    else:

        if 'Downloading GRCm38 gene gtf file' not in file_text:
                download_gtf_GRCm38(path_db.replace('.bed', '.gtf'), file_text)

        cmd_HISAT2_index_GRCm38 = f"""docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['hisat2']} \\\n\
        sh -c "\\\n\
            hisat2_extract_splice_sites.py /database/genomes/GRCm38/genes.gtf > genome.splicesites.txt \n\
            hisat2_extract_exons.py /database/genomes/GRCm38/genes.gtf > genome.exons.txt \n\
            hisat2-build  \\\n\
            -p {MAX_CPU} \\\n\
            --ss genome.splicesites.txt \\\n\
            --exon genome.exons.txt \\\n\
            /{DICT_DBS['genome_fasta_GRCm38'][0]}  \\\n\
            /{DICT_DBS['hisat2_index_GRCm38'][0]}/genome
        "\n\n\n"""


    return title + cmd_HISAT2_index_GRCm38


def download_salmon_index_GRCm38(path_db, file_text):
    subfiles_check = ['complete_ref_lens.bin', 'ctable.bin', 'ctg_offsets.bin', 'duplicate_clusters.tsv', 'info.json', 'mphf.bin', 'pos.bin', 
                      'pre_indexing.log', 'rank.bin', 'refAccumLengths.bin', 'ref_indexing.log', 'reflengths.bin', 'refseq.bin',
                       'seq.bin']
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=500 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCm38 salmon index file."\n\n"""
    path, filename = os.path.split(path_db)

    os.makedirs(path, exist_ok=True)

    if 'Downloading GRCm38 transcript fasta' not in file_text:
        download_gtf_GRCm38(DICT_DBS['transcript_fasta_GRCm38'][0], file_text)

    cmd_salmon_index_GRCm38  = f"""docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['salmon']} \\\n\
        sh -c "\\\n\
            salmon index \\\n\
            -p {MAX_CPU} \\\n\
            -t /{DICT_DBS['transcript_fasta_GRCm38'][0]} \\\n\
            -i /{DICT_DBS['salmon_index_GRCm38'][0]}
        "\n\n\n"""

    return title + cmd_salmon_index_GRCm38



def download_txp2gene_GRCm38(path_db, file_text):
    subfiles_check = []

    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=5 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCm38 salmon txp2gene file."\n\n"""

    if 'Downloading GRCm38 gene gtf file' not in file_text:
        download_gtf_GRCm38(DICT_DBS['gtf_GRCm38'][0], file_text)

    cmd_salmon_txp2gene_GRCm38  = f"""docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['gffread']} \\\n\
        sh -c "\\\n\
            gffread \\\n\
                /{DICT_DBS['gtf_GRCm38'][0]} \\\n\
                --table transcript_id,gene_id \\\n\
                    > /{DICT_DBS['txp2gene_GRCm38'][0]}
        "\n\n\n"""

    return title + cmd_salmon_txp2gene_GRCm38



def download_rsem_index_GRCm38(path_db, file_text):
    subfiles_check = [f'genome.{suffix}' for suffix in ['chrlist', 'grp', 'idx.fa', 'n2g.idx.fa', 'seq', 'ti', 'transcripts.fa']] + \
                     ['chrLength.txt', 'chrNameLength.txt', 'chrName.txt',  'chrStart.txt',  'exonGeTrInfo.tab',  'exonInfo.tab',  'geneInfo.tab',  
                      'Genome',  'genomeParameters.txt',  'Log.out',  'SA',  'SAindex',  'sjdbInfo.txt',  'sjdbList.fromGTF.out.tab',  
                      'sjdbList.out.tab',  'transcriptInfo.tab']
    
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=27 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCm38 rsem index file."\n\n"""
    os.makedirs(path_db, exist_ok=True)

    if 'Downloading GRCm38 genome fasta file' not in file_text:
        download_gtf_GRCm38(DICT_DBS['genome_fasta_GRCm38'][0], file_text)

    if 'Downloading GRCm38 gene gtf file' not in file_text:
        download_gtf_GRCm38(DICT_DBS['gtf_GRCm38'][0], file_text)


    cmd_rsem_index_GRCm38  = f"""docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['rsem-star']} \\\n\
        sh -c "\\\n\
            rsem-prepare-reference \\\n\
            --gtf \{DICT_DBS['gtf_GRCm38'][0]} \\\n\
            --star \\\n\
            --star-path /usr/local/bin \\\n\
            -p {MAX_CPU} \\\n\
            /{DICT_DBS['genome_fasta_GRCm38'][0]} \\\n\
            /{DICT_DBS['rsem_index_GRCm38'][0]}/genome
        "\n\n\n"""

    
    return title + cmd_rsem_index_GRCm38



def download_kallisto_index_GRCm38(path_db, file_text):
    subfiles_check = []
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=180 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCm38 kallisto index file."\n\n"""
    path, filename = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    if 'Downloading GRCm38 gene gtf file' not in file_text:
            download_gtf_GRCm38(path_db.replace('.bed', '.gtf'), file_text)

    cmd_kallisto_index_GRCm38 = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['kallisto']} \\\n\
        sh -c "\\\n\
            kallisto index -i {filename} \\\n\
            -t {MAX_CPU} \\\n\
            /{DICT_DBS['transcript_fasta_GRCm38'][0]} && mv {filename} /{path}/{filename}
    "\n\n\n"""


    return title + cmd_kallisto_index_GRCm38


def download_cellranger_index_GRCm38(path_db, file_text):
    subfiles_check = ['fasta/genome.fa', 'fasta/genome.fa.fai', 'genes/genes.gtf.gz', 'reference.json'] + \
                     [f'star/{file}' for file in ['chrLength.txt', 'chrNameLength.txt', 'chrName.txt', 'chrStart.txt', 'exonGeTrInfo.tab', 
                                                  'exonInfo.tab', 'geneInfo.tab', 'Genome', 'genomeParameters.txt', 'SA', 'SAindex',  'sjdbInfo.txt', 
                                                  'sjdbList.fromGTF.out.tab', 'sjdbList.out.tab', 'transcriptInfo.tab']]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=14.5 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCm38 cellranger index file."\n\n"""
    path, filename = os.path.split(path_db)
    os.makedirs(path_db, exist_ok=True)

    if 'Downloading GRCm38 genome fasta file' not in file_text:
        download_gtf_GRCm38(DICT_DBS['genome_fasta_GRCm38'][0], file_text)
    
    if 'Downloading GRCm38 gene gtf file' not in file_text:
            download_gtf_GRCm38(path_db.replace('.bed', '.gtf'), file_text)

    cmd_cellranger_index_GRCm38 = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['cellranger']} \\\n\
        sh -c "\\\n\
            cellranger mkref \\\n\
            --genome  {filename} \\\n\
            --fasta /{DICT_DBS['genome_fasta_GRCm38'][0]} \\\n\
            --genes /{DICT_DBS['gtf_GRCm38'][0]} \\\n\
            --nthreads={MAX_CPU}\n\
            mv /cellranger /{path}
    "\n\n\n"""

    return title + cmd_cellranger_index_GRCm38



def download_universc_index_GRCm38(path_db, file_text):
    subfiles_check = ['fasta/genome.fa', 'fasta/genome.fa.fai', 'genes/genes.gtf.gz', 'reference.json'] + \
                     [f'star/{file}' for file in ['chrLength.txt', 'chrNameLength.txt', 'chrName.txt', 'chrStart.txt', 'exonGeTrInfo.tab', 
                                                  'exonInfo.tab', 'geneInfo.tab', 'Genome', 'genomeParameters.txt', 'SA', 'SAindex',  'sjdbInfo.txt', 
                                                  'sjdbList.fromGTF.out.tab', 'sjdbList.out.tab', 'transcriptInfo.tab']]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=14.5 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCm38 universc index file."\n\n"""
    path, filename = os.path.split(path_db)
    os.makedirs(path_db, exist_ok=True)

    if 'Downloading GRCm38 genome fasta file' not in file_text:
        download_gtf_GRCm38(DICT_DBS['genome_fasta_GRCm38'][0], file_text)
    
    if 'Downloading GRCm38 gene gtf file' not in file_text:
            download_gtf_GRCm38(path_db.replace('.bed', '.gtf'), file_text)

    cmd_universc_index_GRCm38 = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['cellranger']} \\\n\
        sh -c "\\\n\
            cellranger mkref \\\n\
            --genome  {filename} \\\n\
            --fasta /{DICT_DBS['genome_fasta_GRCm38'][0]} \\\n\
            --genes /{DICT_DBS['gtf_GRCm38'][0]} \\\n\
            --nthreads={MAX_CPU}\n\
            mv /universc /{path}
    "\n\n\n"""

    return title + cmd_universc_index_GRCm38


def download_segemehl_index_GRCm38(path_db, file_text):
    subfiles_check = ['genome.idx']
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=40 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCm38 segemehl index file."\n\n"""
    os.makedirs(path_db, exist_ok=True)

    if 'Downloading GRCm38 genome fasta file' not in file_text:
        download_gtf_GRCm38(DICT_DBS['genome_fasta_GRCm38'][0], file_text)
    
    if 'Downloading GRCm38 gene gtf file' not in file_text:
            download_gtf_GRCm38(path_db.replace('.bed', '.gtf'), file_text)

    cmd_segemehl_index_GRCm38 = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['segemehl']} \\\n\
        sh -c "\\\n\
            segemehl.x \\\n\
            -t {MAX_CPU} \\\n\
            -d /{DICT_DBS['genome_fasta_GRCm38'][0]} \\\n\
            -x /{path_db}/genome.idx 
    "\n\n\n"""


    return title + cmd_segemehl_index_GRCm38



def download_aa_data_repo_index_GRCm38(path_db, file_text):
    subfiles_check = ['mm10/' + i for i in ['annotations/mm10GenomicSuperDup.tab', 'cancer/oncogene_list.txt', 
                      'dummy_ploidy.vcf', 'file_list.txt', 'file_sources.txt', 'last_updated.txt', 'mm10-blacklist.v2.bed', 
                      'mm10_cnvkit_filtered_ref.cnn', 'mm10_conserved_gain5_onco_subtract.bed', 'mm10.fa.amb', 'mm10.fa.bwt', 
                      'mm10.fa.pac', 'mm10.Hardison.Excludable.full.bed', 'mm10_merged_centromeres_conserved_sorted.bed', 'onco_bed.bed', 
                      'mm10_centromere.bed', 'mm10_conserved_gain5.bed', 'mm10.fa', 'mm10.fa.ann', 'mm10.fa.fai', 'mm10.fa.sa', 
                      'mm10_k35.mappability.bedgraph', 'mm10_noAlt.fa.fai']]
                      
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=8 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building GRCm38 ampliconarchitect index file."\n\n"""
    os.makedirs(path_db, exist_ok=True)

    cmd_aa_repo_index_GRCm38 = f""" docker run -it --rm -v $(pwd)/database/:/database \\\n\
        {DICT_CONTAINERS['aria2c']} \\\n\
         -x {MAX_CPU} \\\n\
         -d /{path_db} \\\n\
         --file-allocation=none \\\n\
         https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/mm10_indexed.tar.gz \n\
         tar xzvf {path_db}/mm10_indexed.tar.gz -C {path_db}\n\
         rm -rf {path_db}/mm10_indexed.tar.gz 
    \n\n\n"""

    return title + cmd_aa_repo_index_GRCm38



















def download_mirbase_mature(path_db, file_text):
    subfiles_check = []
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=3 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Downloading mirbase mature file."\n\n"""
    path, _ = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    cmd = f"""wget https://mirbase.org/download/mature.fa -O {path_db} \n\n\n"""

    return title + cmd


def download_mirbase_hairpin(path_db, file_text):
    subfiles_check = []
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=6 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Downloading mirbase hairpin file."\n\n"""
    path, _ = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    cmd = f"""wget https://mirbase.org/download/hairpin.fa -O {path_db} \n\n\n"""

    return title + cmd


def download_mirbase_gtf_GRCh38(path_db, file_text):
    subfiles_check = []
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=0.3 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Downloading mirbase hsa.gff file."\n\n"""
    path, _ = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    cmd = f"""wget https://mirbase.org/download/hsa.gff3 -O {path_db} \n\n\n"""

    return title + cmd


def download_mirbase_gtf_GRCm38(path_db, file_text):
    subfiles_check = []
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=0.3 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Downloading mirbase mmu.gff file."\n\n"""
    path, _ = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    cmd = f"""wget https://mirbase.org/download/mmu.gff3 -O {path_db} \n\n\n"""


    return title + cmd





def download_mirgenedb_mature(path_db, file_text):
    subfiles_check = []
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=0.02 * MB): # TOO SMALL
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Downloading mirgenedb mature file."\n\n"""
    path, _ = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    cmd = f"""wget https://mirgenedb.org/fasta/hsa?mat=1 -O {path_db} \n\n\n"""

    return title + cmd


def download_mirgenedb_hairpin(path_db, file_text):
    subfiles_check = []
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=0.04 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Downloading mirgenedb hairpin file."\n\n"""
    path, _ = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    cmd = f"""wget https://mirgenedb.org/static/data/hsa/hsa-pre.fas -O {path_db} \n\n\n"""

    return title + cmd



def download_mirgenedb_gff_GRCh38(path_db, file_text):
    subfiles_check = []
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=0.1 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Downloading mirgenedb hsa.gff file."\n\n"""
    path, _ = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)

    cmd = f"""wget https://mirgenedb.org/gff/hsa?sort=pos -O {path_db} \n\n\n"""


    return title + cmd


def download_mirgenedb_gff_GRCm38(path_db, file_text):
    subfiles_check = []
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=0.1 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Downloading mirgenedb mmu.gff file."\n\n"""
    path, _ = os.path.split(path_db)
    os.makedirs(path, exist_ok=True)


    cmd = f"""wget https://mirgenedb.org/gff/mmu?sort=pos -O {path_db} \n\n\n"""

    return title + cmd





def download_kaiju(path_db, file_text):
    subfiles_check = [f'kaiju_{org}.fmi' for org in ['fungi', 'human', 'plasmids', 'refseq']]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=89 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building kaiju index files."\n\n"""
    os.makedirs(path_db, exist_ok=True)


    cmd_download_refseq = f"""wget -L https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_refseq_2023-05-23.tgz -O {path_db}/kaiju_refseq.tgz \n\
                    tar xzvf {path_db}/kaiju_refseq.tgz -C {path_db} \n\
                    mv {path_db}/kaiju_db_refseq.fmi {path_db}/kaiju_refseq.fmi \n\
    \n"""

    cmd_download_fungi = f"""wget -L https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_fungi_2023-05-26.tgz -O {path_db}/kaiju_fungi.tgz \n\
                    tar xzvf {path_db}/kaiju_fungi.tgz -C {path_db} \n\
                    mv {path_db}/kaiju_db_fungi.fmi {path_db}/kaiju_fungi.fmi \n\
    \n"""

    cmd_download_plasmids = f"""wget -L https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_plasmids_2023-05-26.tgz -O {path_db}/kaiju_plasmids.tgz \n\
                    tar xzvf {path_db}/kaiju_plasmids.tgz -C {path_db} \n\
                    mv {path_db}/kaiju_db_plasmids.fmi {path_db}/kaiju_plasmids.fmi \n\
    \n"""


    cmd_download_human = f"""docker run -it --rm -v $(pwd)/database/:/database \\\n\
            {DICT_CONTAINERS['ncbi-genome-download']} \\\n\
            sh -c " \\\n\
            ncbi-genome-download -F protein-fasta -p {MAX_CPU} -r 10 -P vertebrate_mammalian -t "9606" -R reference \\\n\
            -o {path_db} \n\
            
            zcat {path_db}/refseq/vertebrate_mammalian/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_protein.faa.gz | \
                sed 's/^>.*$/>9606/' > {path_db}/human.faa \n\
    "\n"""

    cmd_make_human = f"""docker run -it --rm -v $(pwd)/database/:/database \\\n\
            {DICT_CONTAINERS['kaiju']} \\\n\
            sh -c " \\\n\
            kaiju-mkbwt -n {MAX_CPU} -a ACDEFGHIKLMNPQRSTVWY -infilename {path_db}/human.faa -o {path_db}/kaiju_human
            kaiju-mkfmi {path_db}/kaiju_human
    "\n"""

    cmd_rm = f"""rm -rf {path_db}/*.tgz {path_db}/refseq {path_db}/*.bwt {path_db}/*.sa  {path_db}/*.faa  \n\n\n"""

    return title + cmd_download_refseq + cmd_download_fungi + cmd_download_plasmids + cmd_download_human + cmd_make_human + cmd_rm



def download_centrifuge(path_db, file_text):
    subfiles_check = [f'f+p.{N}.cf' for N in range (1,5)] + [f'nt.{N}.cf' for N in range (1,5)] + [f'p+h+v.{N}.cf' for N in range (1,4)]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=130 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building centrifuge index files."\n\n"""
    os.makedirs(path_db, exist_ok=True)


    cmd_download_centrifuge_1 = f"""docker run -it --rm -v $(pwd)/database/:/database \
                {DICT_CONTAINERS['aws-cli']} \\\n\
                s3 --no-sign-request --region eu-west-1  \\\n\
                cp s3://genome-idx/centrifuge/nt_2018_3_3.tar.gz \\\n\
                /{path_db}/nt_2018_3_3.tar.gz \n\
                tar xvf {path_db}/nt_2018_3_3.tar.gz -C {path_db}
    \n"""

    cmd_download_centrifuge_2 = f"""docker run -it --rm -v $(pwd)/database/:/database \
                {DICT_CONTAINERS['aws-cli']} \\\n\
                s3 --no-sign-request --region eu-west-1  \\\n\
                cp s3://genome-idx/centrifuge/p+h+v.tar.gz \\\n\
                /{path_db}/p+h+v.tar.gz \n\
                tar xvf {path_db}/p+h+v.tar.gz -C {path_db}
    \n"""

    # I CANT RUN DIRECTLY CENTRIFUGE FROM DOCKER BECUASE IT NEEDS DUSTMASKER 
    # TODO: CREATE A CENTRIFUGE + DUSTMASKER 
    cmd_build_fungi_protozoa = f"""centrifuge-download -o {path_db}/taxonomy taxonomy \n\
            centrifuge-download -o {path_db}/library -P {MAX_CPU} -m -d "fungi,protozoa" -a Any refseq > {path_db}/seqid2taxid.map \n\
            cat {path_db}/library/*/*.fna > {path_db}/input-sequences.fna \n\
            centrifuge-build -p {MAX_CPU} --conversion-table {path_db}/seqid2taxid.map \\\n\
                 --taxonomy-tree {path_db}/taxonomy/nodes.dmp --name-table {path_db}/taxonomy/names.dmp \\\n\
                 {path_db}/input-sequences.fna {path_db}/f+p
    \n"""

    cmd_rm = f"""rm -rf {path_db}/*.tar.gz {path_db}/taxonomy {path_db}/library {path_db}/*.map {path_db}/*.fna\n\n\n"""

    return title + cmd_download_centrifuge_1 + cmd_download_centrifuge_2 + cmd_build_fungi_protozoa + cmd_rm





def download_kraken2(path_db, file_text):
    subfiles_check = ['hash.k2d', 'inspect.txt', 'kraken_2_db_inspect.txt', 'ktaxonomy.tsv', 'opts.k2d', 'seqid2taxid.map', 'taxo.k2d'] + \
                     [f'database{N}mers.kmer_distrib' for N in [50, 75, 100, 150, 200, 250, 300]]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=70 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building kraken2 index files."\n\n"""
    os.makedirs(path_db, exist_ok=True)


    cmd_download_kraken2 = f"""wget -L https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20230314.tar.gz \
                                -O {path_db}/kraken_2_db.tar.gz \n\
                    wget -L https://genome-idx.s3.amazonaws.com/kraken/pluspf_20230314/inspect.txt \
                                -O {path_db}/kraken_2_db_inspect.txt \n\
                    tar xzvf {path_db}/kraken_2_db.tar.gz -C {path_db} \n\
    \n"""


    cmd_rm = f"""rm -rf {path_db}/*.tar.gz \n\n\n"""

    return title + cmd_download_kraken2 + cmd_rm


def download_krakenuniq(path_db, file_text):
    subfiles_check = [f'database.{suff}' for suff in ['idx', 'kdb', 'kdb.counts']] + ['seqid2taxid.map', 'taxDB'] + \
                     [f'database{N}mers.kmer_distrib' for N in [50, 75, 100, 150, 200, 250, 300]]

    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=400 * GB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Building krakenuniq index files."\n\n"""
    os.makedirs(path_db, exist_ok=True)

    cmd_download_krakenuniq_1 = f"""docker run -it --rm -v $(pwd)/database/:/database \
                {DICT_CONTAINERS['aws-cli']} \\\n\
                s3 --no-sign-request --region eu-west-1  \\\n\
                cp s3://genome-idx/kraken/uniq/krakendb-2022-06-16-STANDARD/kuniq_standard_minus_kdb.20220616.tgz \\\n\
                /{path_db}/kuniq_standard_minus_kdb.20220616.tgz \n\
                tar xvf {path_db}/kuniq_standard_minus_kdb.20220616.tgz -C {path_db}
    \n"""

    cmd_download_krakenuniq_2 = f"""docker run -it --rm -v $(pwd)/database/:/database \
                {DICT_CONTAINERS['aws-cli']} \\\n\
                s3 --no-sign-request --region eu-west-1  \\\n\
                cp s3://genome-idx/kraken/uniq/krakendb-2022-06-16-STANDARD/database.kdb \\\n\
                /{path_db}/database.kdb
    \n"""

    cmd_rm = f"""rm -rf {path_db}/*.tgz \n\n\n"""

    return title + cmd_download_krakenuniq_1 + cmd_download_krakenuniq_2 + cmd_rm


def download_taxpasta(path_db, file_text):
    subfiles_check = ['citations.dmp', 'delnodes.dmp', 'division.dmp','gencode.dmp', 'images.dmp','merged.dmp',
                      'names.dmp','nodes.dmp','gc.prt','readme.txt'
]
    
    if not db_check(path_db, subfiles_check=subfiles_check, min_weight=400 * MB):
        logger.info(f'Database {path_db} already exists. It will not be downloaded.')
        return ""  # I write the code as such to avoid encapsulation.

    title = """echo "### Downloading taxpasta files."\n\n"""
    path, _ = os.path.split(path_db)
    os.makedirs(path_db, exist_ok=True)


    cmd = f"""wget ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O {path_db}/taxdump.tar.gz \n\
              tar xvf {path_db}/taxdump.tar.gz -C {path_db}/  \n\
              rm  {path_db}/*.tar.gz \n\n\n"""

    return title + cmd



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
                if not os.path.exists(f"{path_check}/{file}"):
                    logger.info(f">>> File {file} within path {path_check} does not exist. Database will be downloaded.")
                    return True

    # Check dir/file size
    if len(subfiles_check) > 0: # is a directory
        weight = sum([os.path.getsize(f"{path_check}/{f}") for f in subfiles_check if os.path.exists(f"{path_check}/{f}")])
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
DICT_DBS = {'genome_fasta_GRCh38': ('database/genomes/GRCh38/genome.fasta', download_fasta_GRCh38), 
            'transcript_fasta_GRCh38': ('database/genomes/GRCh38/transcript.fasta', download_transcript_fasta_GRCh38), 
            'gtf_GRCh38': ('database/genomes/GRCh38/genes.gtf', download_gtf_GRCh38), 
            'gtf_corrected_GRCh38': ('database/genomes/GRCh38/genes_corrected.gtf', download_gtf_corrected_GRCh38), 
            'bed_GRCh38': ('database/genomes/GRCh38/genes.bed', download_bed_GRCh38), 
            'star_index_GRCh38': ('database/indexes/GRCh38/STAR', download_star_index_GRCh38),
            'bowtie_index_GRCh38': ('database/indexes/GRCh38/BOWTIE', download_bowtie_index_GRCh38),
            'bowtie2_index_GRCh38': ('database/indexes/GRCh38/BOWTIE2', download_bowtie2_index_GRCh38), 
            'bwa_index_GRCh38': ('database/indexes/GRCh38/BWA', download_bwa_index_GRCh38), 
            'hisat2_index_GRCh38': ('database/indexes/GRCh38/HISAT2', download_hisat2_index_GRCh38), 
            'salmon_index_GRCh38': ('database/indexes/GRCh38/SALMON', download_salmon_index_GRCh38), 
            'txp2gene_GRCh38': ('database/indexes/GRCh38/SALMON/txp2gene', download_txp2gene_GRCh38),
            'rsem_index_GRCh38': ('database/indexes/GRCh38/rsem', download_rsem_index_GRCh38), 
            'kallisto_index_GRCh38': ('database/indexes/GRCh38/KALLISTO/index', download_kallisto_index_GRCh38), 
            'cellranger_index_GRCh38': ('database/indexes/GRCh38/cellranger', download_cellranger_index_GRCh38),
            'universc_index_GRCh38': ('database/indexes/GRCh38/universc', download_universc_index_GRCh38), 
            'segemehl_index_GRCh38': ('database/indexes/GRCh38/segemehl', download_segemehl_index_GRCh38), 
            'aa_data_repo_index_GRCh38': ('database/indexes/GRCh38/aa_data_repo', download_aa_data_repo_index_GRCh38),
            'genome_fasta_CHM13': ('database/genomes/CHM13/genome.fasta', download_fasta_CHM13),
            'bowtie2_index_CHM13': ('database/indexes/CHM13/BOWTIE2', download_bowtie2_index_CHM13),
            'genome_fasta_GRCm38': ('database/genomes/GRCm38/genome.fasta', download_fasta_GRCm38), 
            'transcript_fasta_GRCm38': ('database/genomes/GRCm38/transcript.fasta', download_transcript_fasta_GRCm38),
            'gtf_GRCm38': ('database/genomes/GRCm38/genes.gtf', download_gtf_GRCm38), 
            'gtf_corrected_GRCm38': ('database/genomes/GRCm38/genes_corrected.gtf', download_gtf_corrected_GRCm38), 
            'bed_GRCm38': ('database/genomes/GRCm38/genes.bed', download_bed_GRCm38), 
            'star_index_GRCm38': ('database/indexes/GRCm38/STAR', download_star_index_GRCm38),
            'bowtie_index_GRCm38': ('database/indexes/GRCm38/BOWTIE', download_bowtie_index_GRCm38),
            'bowtie2_index_GRCm38': ('database/indexes/GRCm38/BOWTIE2', download_bowtie2_index_GRCm38), 
            'bwa_index_GRCm38': ('database/indexes/GRCm38/BWA', download_bwa_index_GRCm38), 
            'hisat2_index_GRCm38': ('database/indexes/GRCm38/HISAT2', download_hisat2_index_GRCm38), 
            'salmon_index_GRCm38': ('database/indexes/GRCm38/SALMON', download_salmon_index_GRCm38), 
            'txp2gene_GRCm38': ('database/indexes/GRCm38/SALMON/txp2gene', download_txp2gene_GRCm38),
            'rsem_index_GRCm38': ('database/indexes/GRCm38/rsem', download_rsem_index_GRCm38), 
            'kallisto_index_GRCm38': ('database/indexes/GRCm38/KALLISTO/index', download_kallisto_index_GRCm38),
            'cellranger_index_GRCm38': ('database/indexes/GRCm38/cellranger', download_cellranger_index_GRCm38), 
            'universc_index_GRCm38': ('database/indexes/GRCm38/universc', download_universc_index_GRCm38), 
            'segemehl_index_GRCm38': ('database/indexes/GRCm38/segemehl', download_segemehl_index_GRCm38), 
            'aa_data_repo_index_GRCm38': ('database/indexes/GRCm38/aa_data_repo', download_aa_data_repo_index_GRCm38), 
            'mirbase_gff_GRCh38': ('database/smRNA/mirbase/hsa.gff3', download_mirbase_gtf_GRCh38),
            'mirbase_gff_GRCm38': ('database/smRNA/mirbase/mmu.gff3', download_mirbase_gtf_GRCm38),
            'mirgenedb_gff_GRCh38': ('database/smRNA/mirgenedb/hsa.gff', download_mirgenedb_gff_GRCh38),
            'mirgenedb_gff_GRCm38': ('database/smRNA/mirgenedb/mmu.gff', download_mirgenedb_gff_GRCm38),
            'mirbase_mature': ('database/smRNA/mirbase/mature.fa', download_mirbase_mature),
            'mirbase_hairpin': ('database/smRNA/mirbase/hairpin.fa', download_mirbase_hairpin),
            'mirgenedb_mature': ('database/smRNA/mirgenedb/hsa.fas', download_mirgenedb_mature),
            'mirgenedb_hairpin': ('database/smRNA/mirgenedb/hsa-pre.fas', download_mirgenedb_hairpin),
            'kaiju': ('database/taxprofiler/kaiju', download_kaiju), 
            'centrifuge': ('database/taxprofiler/centrifuge', download_centrifuge), 
            'kraken2': ('database/taxprofiler/kraken2', download_kraken2), 
            'krakenuniq': ('database/taxprofiler/krakenuniq', download_krakenuniq), 
            'taxpasta': ('database/taxprofiler/taxpasta', download_taxpasta)}
