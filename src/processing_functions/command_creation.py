import os
import pandas as pd

from .versions import VERSION
from .global_vars_and_funcs import NF_CONFIG_PIPELINES, DICT_GENOMES, MAX_RAM, MAX_CPU, MAX_TIME
from .yaml_config import retrieve_r_tag, check_entry, read_yaml
#from .yaml_config import MAX_MEMORY_PER, MAX_CPU_PER, MAX_TIME, get_available_cpus, get_RAM_GB
from .database_download import DICT_DBS





################################################################################################
# GLOBAL VARIABLES
################################################################################################








################################################################################################
# PRIMARY FUNCTIONS
################################################################################################
def create_command_sh_file(yaml_dict, project):
    """
    Create a command.sh file with all the necessary commands to run each process. This will be run directly afterwards.
    """

    file_text = \
f"""
echo "
     ___           ___           ___                                 
    /\  \         /\__\         /\__\                                
    \:\  \       /:/ _/_       /:/ _/_                               
     \:\  \     /:/ /\  \     /:/ /\  \                              
 _____\:\  \   /:/ /::\  \   /:/ /::\  \                             
/::::::::\__\ /:/__\/\:\__\ /:/_/:/\:\__\                            
\:\~~\~~\/__/ \:\  \ /:/  / \:\/:/ /:/  /                            
 \:\  \        \:\  /:/  /   \::/ /:/  /                             
  \:\  \        \:\/:/  /     \/_/:/  /                              
   \:\__\        \::/  /        /:/  /                               
    \/__/         \/__/         \/__/                                
     ___                     ___         ___                         
    /\  \                   /\  \       /\__\                        
   /::\  \     ___         /::\  \     /:/ _/_                       
  /:/\:\__\   /\__\       /:/\:\__\   /:/ /\__\                      
 /:/ /:/  /  /:/__/      /:/ /:/  /  /:/ /:/ _/_                     
/:/_/:/  /  /::\  \     /:/_/:/  /  /:/_/:/ /\__\                    
\:\/:/  /   \/\:\  \__  \:\/:/  /   \:\/:/ /:/  /                    
 \::/__/       \:\/\__\  \::/__/     \::/_/:/  /                     
  \:\  \        \::/  /   \:\  \      \:\/:/  /                      
   \:\__\       /:/  /     \:\__\      \::/  /                       
    \/__/       \/__/       \/__/       \/__/                        
     ___           ___           ___           ___           ___     
    /\  \         /\__\         /\  \         /\__\         /\__\    
   |::\  \       /:/ _/_       /::\  \       /:/ _/_       /:/ _/_   
   |:|:\  \     /:/ /\__\     /:/\:\__\     /:/ /\  \     /:/ /\__\  
 __|:|\:\  \   /:/ /:/ _/_   /:/ /:/  /    /:/ /::\  \   /:/ /:/ _/_ 
/::::|_\:\__\ /:/_/:/ /\__\ /:/_/:/__/___ /:/__\/\:\__\ /:/_/:/ /\__\  
\:\~~\  \/__/ \:\/:/ /:/  / \:\/::::/__ / \:\  \ /:/  / \:\/:/ /:/  /
 \:\  \        \::/_/:/  /   \::/~~/       \:\  /:/  /   \::/_/:/  / 
  \:\  \        \:\/:/  /     \:\~~\        \:\/:/  /     \:\/:/  /  
   \:\__\        \::/  /       \:\__\        \::/  /       \::/  /   
    \/__/         \/__/         \/__/         \/__/         \/__/    

------------------------------------------------------------------
| NGS PIPE MERGE v{VERSION}
------------------------------------------------------------------

------------------------------------------------------------------
| PIPELINE LIST: {', '.join(list(yaml_dict.keys()))}
------------------------------------------------------------------
"
"""

    # When parsing arguments we will find 3 main options:
    #   1) argument: value  -> we pass the text as it is:  --argument value
    #   2) argument: True   -> we pass it as --argument
    #   3) argument: False  -> we don't pass it
    for process_name in yaml_dict.keys():
        process_dict = yaml_dict[process_name]
        pipeline = process_dict['general_config']['pipeline']

        command_process = f"""echo -e "{'=' * 50}\nPROCESS:  {process_name}\n{'=' * 50}\n" \n\n"""

        #TODO: IMPRIMIR LOS ARGUMENTOS 


        if pipeline in NF_CONFIG_PIPELINES:
            command_process += f"nextflow run nf-core/{pipeline} \\\n\
                                            -r {process_dict['nextflow_config']['r']} \\\n\
                                            -profile {process_dict['nextflow_config']['profile']} \\\n\
                                            {'-resume' if process_dict['nextflow_config']['resume'] else ''} \\\n\
                                            --max_cpus {process_dict['general_config']['max_cpus']} \\\n\
                                            --max_memory {process_dict['general_config']['max_memory']}.GB \\\n\
                                            --max_time {process_dict['general_config']['max_time']}.h \\\n"
            
            for key, val in process_dict['nfcore_config'].items():
                command_process += process_arg_keys(key, val, pipeline)

        elif pipeline == 'taxprofiler':
            # 1st and 2nd mappings
            # Depending on wether 1st and 2nd mappings are activated, we have to root the inputs and outputs of the pipelines.
            # The main problem is that the names of the fastqs and the names of the pools may not be the same, so we have to redirect
            # the inputs and outputs accordingly. We will use the spreadsheets of 1stmap and 2ndmap for that use
            host_config = yaml_dict[process_name]['host_mapping_config']
            map_1, map_2 = host_config['1st_map'], host_config['2nd_map']
            organism = yaml_dict[process_name]['general_config']['organism']
            genome = DICT_GENOMES[organism]
            # Depending on the mapping map_1 and map_2 options we will reroute the inputs and outputs of the processing.
            # The arguments of the commands called may vary depending on whether the reads are paired or unpaired
            # TODO: CREATE COMMAND FOR UNPAIRED READS!!!

            # Make dirs
            if map_1:
                os.makedirs(f'results/{project}/1st_map', exist_ok=True)
            if map_2:
                os.makedirs(f'results/{project}/2nd_map', exist_ok=True)

            # Run the mapping processes
            if map_1 & map_2:
                command_process += run_map_1(genome, project, pipeline, process_name, yaml_dict)
                command_process += run_map_2(map_1, genome, project, pipeline, process_name, yaml_dict)
            elif map_1 & (not map_2):
                command_process += run_map_1(genome, project, pipeline, process_name, yaml_dict)
            elif (not map_1) & map_2:
                command_process += run_map_2(map_1, genome, project, pipeline, process_name, yaml_dict)
            else:
                pass



            # PROFILERS
            # MIRAR COMO SE PASAN LOS ARGUMENTOS EN BASE AL YAML CONFIG COMPLETO

            list_profilers_str = yaml_dict[process_name]['profiler_config']['profilers']
            list_profilers = list_profilers_str.split(',')

            for profiler in list_profilers:
                os.makedirs(f'results/{project}/{pipeline}/{profiler}', exist_ok=True)

                # Depending on whether map_1 or map_2 have been applied, we are going to create different input values.
                # We are going to apply the same pattern, either by reading MAP1, MAP2 or INPUT samplesheet.

                if (map_1 & map_2) | (not map_1) & map_2:
                    samplesheet_path = f"work/{project}/samplesheets/{process_name}_map_2.csv"
                elif map_1 & (not map_2):
                    samplesheet_path = f"work/{project}/samplesheets/{process_name}_map_1.csv"
                else:
                    samplesheet_path = f"work/{project}/samplesheets/{process_name}.csv"
                    
                list_samples, list_fastqs_1, list_fastqs_2 = return_list_fastqs(pd.read_csv(samplesheet_path))


                match profiler:
                    case 'kaiju':
                        command_kaiju = run_kaiju(yaml_dict, process_name, pipeline, list_samples, list_fastqs_1, list_fastqs_2)
                        command_process += command_kaiju
                    case 'kraken2':
                        command_kraken2 = run_kraken2(yaml_dict, process_name, pipeline, list_samples, list_fastqs_1, list_fastqs_2)
                        command_process += command_kraken2
                    case 'krakenuniq':
                        command_krakenuniq = run_krakenuniq(yaml_dict, process_name, pipeline, list_samples, list_fastqs_1, list_fastqs_2)
                        command_process += command_krakenuniq
                    case 'centrifuge':
                        command_centrifuge = run_centrifuge(yaml_dict, process_name, pipeline, list_samples, list_fastqs_1, list_fastqs_2)
                        command_process += command_centrifuge
            
            
            command_taxpasta = run_taxpasta(yaml_dict, list_profilers, process_name, pipeline, list_samples)     
            command_process += command_taxpasta       

        # WRITE COMMAND
        if '\\' in command_process[: -3]: # To remove the last "\" in the command
            command_process =  command_process[: -3]

        file_text += command_process + '\n\n\n'  


    # YASSSSSSSS
    write_file(file_text, project)  
 








################################################################################################
# SECONDARY FUNCTIONS
################################################################################################
def run_map_1(genome, project, pipeline, process_name, yaml_dict):

    samplesheet_in = f"work/{project}/samplesheets/{process_name}.csv"


    command_run_1 = f"""echo -e "{'^' * 50}\nPROCESS:  {process_name} | 1ST MAP \n{'^' * 50}\n"  \n\n"""

    command_run_1 += f"nextflow run nf-core/rnaseq -r {retrieve_r_tag('rnaseq')} -profile docker -resume  \\\n\
    --input {samplesheet_in} \\\n\
    --outdir results/{project}/{pipeline}/1st_map \\\n\
    --aligner star_salmon \\\n\
    --fasta {DICT_DBS[f'genome_fasta_{genome}'][0]}\\\n\
    --gtf {DICT_DBS[f'gtf_{genome}'][0]} \\\n\
    --star_index {DICT_DBS[f'star_index_{genome}'][0]}\\\n\
    --rsem_index {DICT_DBS[f'rsem_index_{genome}'][0]} \\\n\
    --salmon_index {DICT_DBS[f'salmon_index_{genome}'][0]}\\\n\
    --extra_salmon_quant_args '--minAssignedFrags 1' \\\n\
    --skip_bbsplit --skip_qualimap --skip_pseudo_alignment --save_unaligned \\\n"

    if check_entry(yaml_dict[process_name]['host_mapping_config'], '1st_map_extra_args'):
        val = yaml_dict[process_name]['host_mapping_config']['1st_map_extra_args']
        val_strip = val.strip().strip("'").strip('"').strip()
        command_run_1 += f'"{val_strip}"'

    command_run_1 += f"--max_cpus {MAX_CPU} --max_memory {MAX_RAM}.GB \\\n\
                       --max_time {MAX_TIME}.h \n\n"
    
    return command_run_1



def run_map_2(map_1, genome, project, pipeline, process_name, yaml_dict):

    command_run_2 = f"""echo -e "{'^' * 50}\nPROCESS:  {process_name} | 2ND MAP \n{'^' * 50}\n"  \n\n"""

    
    if map_1:
        samplesheet_in = f"work/{project}/samplesheets/{process_name}_map_1.csv"
        list_samples, list_fastqs_1, list_fastqs_2 = return_list_fastqs(pd.read_csv(samplesheet_in))
    else:
        samplesheet_in = f"work/{project}/samplesheets/{process_name}.csv"
        list_samples, list_fastqs_1, list_fastqs_2 = return_list_fastqs(pd.read_csv(samplesheet_in))


    command_run_2 += f"""SAMPLES=({' '.join([f'"{i}"' for i in list_samples])}) \n"""
    command_run_2 += f"""FASTQS1=({' '.join([f'"{i}"' for i in list_fastqs_1])}) \n"""
    command_run_2 += f"""FASTQS2=({' '.join([f'"{i}"' for i in list_fastqs_2])}) \n\n"""

    command_run_2 += f"""for ((i=0; i<${{#SAMPLES[@]}}; i++)) \n""" 
    command_run_2 += f"do \n"

    command_run_2 += f"echo ALIGNING ${{SAMPLES[i]}} WITH BOWTIE2 \\\n\
    bowtie2 -x {DICT_DBS[f'bowtie2_index_{genome}'][0]} \\\n\
    -1 ${{FASTQS1[i]}} \\\n\
    -2 ${{FASTQS2[i]}} \\\n\
    -S results/{project}/{pipeline}/2nd_map/${{SAMPLES[i]}}.aligned.sam \\\n"

    if check_entry(yaml_dict[process_name]['host_mapping_config'], '2nd_map_extra_args'):
        val = yaml_dict[process_name]['host_mapping_config']['2nd_map_extra_args']
        val_strip = val.strip().strip("'").strip('"').strip()
        command_run_2 += f'"{val_strip}"'

    command_run_2 +=  f"    --un-conc-gz results/{project}/{pipeline}/2nd_map/${{SAMPLES[i]}}.unmapped.fastq.gz \\\n\
    --very-sensitive -p {MAX_CPU} \\\n\
    done \n\n"

    return command_run_2


# TODO: outsource profiler commands using quay.io containers
def run_kaiju(yaml_dict, process_name, pipeline, list_samples, list_fastqs_1, list_fastqs_2):
    profiler_type = 'kaiju'
    command_kaiju = f"""echo -e "{'^' * 50}\nPROCESS:  {process_name} | KAIJU \n{'^' * 50}\n"  \n\n"""

    command_kaiju += f"""SAMPLES=({' '.join([f'"{i}"' for i in list_samples])}) \n"""
    command_kaiju += f"""FASTQS1=({' '.join([f'"{i}"' for i in list_fastqs_1])}) \n"""
    command_kaiju += f"""FASTQS2=({' '.join([f'"{i}"' for i in list_fastqs_2])}) \n\n"""
    
    command_kaiju += f"""for ((i=0; i<${{#SAMPLES[@]}}; i++)) \n""" 
    command_kaiju += f"do \n"


    # Check args related to KAIJU PROCESSING first to create the corresponding strings
    e_val_str = f"-E {yaml_dict[process_name]['kaiju_config']['e_value']}" if check_entry(yaml_dict[process_name]['kaiju_config'], 'e_value') else ''
    min_len_str = f"-m {yaml_dict[process_name]['kaiju_config']['minimum_length']}" if check_entry(yaml_dict[process_name]['kaiju_config'], 'minimum_length') else ''
    kaiju_index_str = f"-E {yaml_dict[process_name]['kaiju_config']['kaiju_index']}" if check_entry(yaml_dict[process_name]['kaiju_config'], 'kaiju_index') else f"{DICT_DBS[f'kaiju'][0]}/kaiju_refseq.fmi"   
    kaiju_extra_args_str = f""

    if check_entry(yaml_dict[process_name]['kaiju_config'], 'kaiju_extra_args'):
        val = yaml_dict[process_name]['kaiju_config']['kaiju_extra_args']
        kaiju_extra_args_str = val.strip().strip("'").strip('"').strip()

    

    for kaiju_db_type in ['refseq', 'fungi', 'plasmid', 'human']:
        command_kaiju += f"kaiju  -t {DICT_DBS[f'taxpasta'][0]}/nodes.dmp \\\n -f {kaiju_index_str}/kaiju_{kaiju_db_type}.fmi \\\n\
        -i ${{FASTQS1[i]}} \\\n -j ${{FASTQS2[i]}} \\\n\
        -z {MAX_CPU} {e_val_str} {min_len_str} {kaiju_extra_args_str} \\\n\
        -o results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.{kaiju_db_type}.out \n\n"

    
    command_kaiju += f"echo 'Merging entries...' \n"
    for kaiju_db_type in ['refseq', 'fungi', 'plasmid', 'human']:
        command_kaiju += f"sort -k2,2 results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.{kaiju_db_type}.out \
                             > results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.{kaiju_db_type}.sorted.out"

    kaiju_mergeoutputs_extra_args_str = f""
    if check_entry(yaml_dict[process_name]['kaiju_config'], 'kaiju_mergeoutputs_extra_args'):
        val = yaml_dict[process_name]['kaiju_config']['kaiju_mergeoutputs_extra_args']
        kaiju_mergeoutputs_extra_args_str = val.strip().strip("'").strip('"').strip()

    command_kaiju += f"echo 'Combining merged.out files...' \n"
    command_kaiju += f"kaiju-mergeOutputs \\\n\
    -i results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.refseq.sorted.out \\\n\
    -j results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.fungi.sorted.out \\\n\
    -c 'lca' -t {DICT_DBS[f'taxpasta'][0]}/nodes.dmp {kaiju_mergeoutputs_extra_args_str} \\\n\
    -o results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.refseq-fungi.sorted.out \n\n"

    command_kaiju += f"kaiju-mergeOutputs \\\n\
    -i results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.plasmid.sorted.out \\\n\
    -j results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.human.sorted.out \\\n\
    -c 'lca' -t {DICT_DBS[f'taxpasta'][0]}/nodes.dmp {kaiju_mergeoutputs_extra_args_str} \\\n\
    -o results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.human-plasmid.sorted.out \n\n"

    command_kaiju += f"kaiju-mergeOutputs \\\n\
    -i results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.human-plasmid.sorted.out \\\n\
    -j results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.refseq-fungi.sorted.out \\\n\
    -c 'lca' -t {DICT_DBS[f'taxpasta'][0]}/nodes.dmp {kaiju_mergeoutputs_extra_args_str} \\\n\
    -o results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.merged.out \n\n"


    kaiju2table_extra_args_str = f""
    if check_entry(yaml_dict[process_name]['kaiju_config'], 'kaiju2table_extra_args'):
        val = yaml_dict[process_name]['kaiju_config']['kaiju2table_extra_args']
        kaiju2table_extra_args_str = val.strip().strip("'").strip('"').strip()
    
    min_reads_str = f"-c {yaml_dict[process_name]['kaiju_config']['minimum_reads']}" if check_entry(yaml_dict[process_name]['kaiju_config'], 'minimum_reads') else ''   
    command_kaiju += f"echo 'Creating kaiju table...'"
    command_kaiju += f"kaiju2table -t {DICT_DBS[f'taxpasta'][0]}/nodes.dmp \\\n\
    -n {DICT_DBS[f'taxpasta'][0]}/names.dmp \\\n\
    -r species {min_reads_str} \\\n\
    -e \\\n\
    -p {kaiju2table_extra_args_str}\\\n\
    -o results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.report \\\n\
    results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.merged.out  \\\n"

    

    command_kaiju += f"echo 'Removing .out files...' \n"
    command_kaiju += f"rm results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.*.out \n"

    command_kaiju += f"done \n\n"


    return command_kaiju



def run_kraken2(yaml_dict, process_name, pipeline, list_samples, list_fastqs_1, list_fastqs_2):
    profiler_type = 'kraken2'
    command_kraken2 = f"""echo -e "{'^' * 50}\nPROCESS:  {process_name} | KRAKEN2 \n{'^' * 50}\n"  \n\n"""

    command_kraken2 += f"""SAMPLES=({' '.join([f'"{i}"' for i in list_samples])}) \n"""
    command_kraken2 += f"""FASTQS1=({' '.join([f'"{i}"' for i in list_fastqs_1])}) \n"""
    command_kraken2 += f"""FASTQS2=({' '.join([f'"{i}"' for i in list_fastqs_2])}) \n\n"""
    

    confidence_str = f"--confidence {yaml_dict[process_name]['kraken2_config']['confidence']}" if check_entry(yaml_dict[process_name]['kraken2_config'], 'confidence') else ''
    kraken2_index_str = f"{yaml_dict[process_name]['kraken2_config']['kraken2_index']}" if check_entry(yaml_dict[process_name]['kraken2_config'], 'kraken2_index') else f"{DICT_DBS[f'kraken2'][0]}"   

    kraken2_extra_args_str = f""
    if check_entry(yaml_dict[process_name]['kraken2_config'], 'kraken2_extra_args'):
        val = yaml_dict[process_name]['kraken2_config']['kraken2_extra_args']
        kraken2_extra_args_str = val.strip().strip("'").strip('"').strip()

    command_kraken2 += f"""for ((i=0; i<${{#SAMPLES[@]}}; i++)) \n""" 
    command_kraken2 += f"do \n"


    command_kraken2 += f"kraken2 -db {kraken2_index_str} \\\n\
    --threads {MAX_CPU}\\\n\
    --report results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.report    \\\n\
    --classified-out results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.classified#.fastq   \\\n\
    --unclassified-out results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.unclassified#.fastq   \\\n\
    --output results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.output   \\\n\
    --gzip-compressed  \\\n\
    {confidence_str} {kraken2_extra_args_str}  --paired   \\\n\
    ${{FASTQS1[i]}} \\\n\
    ${{FASTQS2[i]}} \n\n"


    command_kraken2 += f"gzip   results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.classified_1.fastq \\\n\
                                results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.classified_2.fastq \\\n\
                                results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.unclassified_1.fastq \\\n\
                                results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.unclassified_2.fastq \n\n"


    command_kraken2 += f"done \n\n"

    return command_kraken2



def run_krakenuniq(yaml_dict, process_name, pipeline, list_samples, list_fastqs_1, list_fastqs_2):
    profiler_type = 'krakenuniq'
    command_krakenuniq = f"""echo -e "{'^' * 50}\nPROCESS:  {process_name} | KRAKENUNIQ \n{'^' * 50}\n"  \n\n"""

    command_krakenuniq += f"""SAMPLES=({' '.join([f'"{i}"' for i in list_samples])}) \n"""
    command_krakenuniq += f"""FASTQS1=({' '.join([f'"{i}"' for i in list_fastqs_1])}) \n"""
    command_krakenuniq += f"""FASTQS2=({' '.join([f'"{i}"' for i in list_fastqs_2])})\n\n"""
    

    hll_precision_str = f"--hll-precision {yaml_dict[process_name]['krakenuniq_config']['hll_precision']}" if check_entry(yaml_dict[process_name]['krakenuniq_config'], 'hll_precision') else ''
    krakenuniq_index_str = f"-E {yaml_dict[process_name]['krakenuniq_config']['krakenuniq_index']}" if check_entry(yaml_dict[process_name]['krakenuniq_config'], 'krakenuniq_index') else f"{DICT_DBS[f'krakenuniq'][0]}"   

    krakenuniq_extra_args_str = f""
    if check_entry(yaml_dict[process_name]['krakenuniq_config'], 'krakenuniq_extra_args'):
        val = yaml_dict[process_name]['krakenuniq_config']['krakenuniq_extra_args']
        krakenuniq_extra_args_str = val.strip().strip("'").strip('"').strip()

    command_krakenuniq += f"""for ((i=0; i<${{#SAMPLES[@]}}; i++)) \n""" 
    command_krakenuniq += f"do \n"


    command_krakenuniq += f"krakenuniq -db {krakenuniq_index_str} \\\n\
    --threads {MAX_CPU}\\\n\
    --report results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.report    \\\n\
    --classified-out results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.classified#.fastq   \\\n\
    --unclassified-out results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.unclassified#.fastq   \\\n\
    --output results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.output   \\\n\
    --gzip-compressed  {hll_precision_str} {krakenuniq_extra_args_str}  \\\n\
    --paired   \\\n\
    ${{FASTQS1[i]}} \\\n\
    ${{FASTQS2[i]}} \n\n"


    command_krakenuniq += f"gzip   results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.classified.fastq \\\n\
                                results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.unclassified.fastq \n\n"


    command_krakenuniq += f"done \n\n"

    return command_krakenuniq



def run_centrifuge(yaml_dict, process_name, pipeline, list_samples, list_fastqs_1, list_fastqs_2):
    profiler_type = 'centrifuge'

    command_centrifuge = f"""echo -e "{'^' * 50}\nPROCESS:  {process_name} | CENTRIFUGE \n{'^' * 50}\n"  \n\n"""

    command_centrifuge += f"""SAMPLES=({' '.join([f'"{i}"' for i in list_samples])}) \n"""
    command_centrifuge += f"""FASTQS1=({' '.join([f'"{i}"' for i in list_fastqs_1])}) \n"""
    command_centrifuge += f"""FASTQS2=({' '.join([f'"{i}"' for i in list_fastqs_2])}) \n\n"""
    

    minimum_length_str = f"--min-hitlen  {yaml_dict[process_name]['centrifuge_config']['minimum_length']}" if check_entry(yaml_dict[process_name]['centrifuge_config'], 'minimum_length') else ''
    qc_filter_str = f"--qc-filter"if check_entry(yaml_dict[process_name]['centrifuge_config'], 'qc_filter') else ''
    centrifuge_index_str = f"-E {yaml_dict[process_name]['centrifuge_config']['centrifuge_index']}" if check_entry(yaml_dict[process_name]['centrifuge_config'], 'centrifuge_index') else f"{DICT_DBS[f'centrifuge'][0]}"   

    centrigue_extra_args_str = f""
    if check_entry(yaml_dict[process_name]['centrifuge_config'], 'centrifuge_extra_args'):
        val = yaml_dict[process_name]['centrifuge_config']['centrifuge_extra_args']
        centrigue_extra_args_str = val.strip().strip("'").strip('"').strip()

    command_centrifuge += f"""for ((i=0; i<${{#SAMPLES[@]}}; i++)) \n""" 
    command_centrifuge += f"do \n"


    command_centrifuge += f"centrifuge \
    -x {centrifuge_index_str}/p+h+v \\\n\
    -1 ${{FASTQS1[i]}} \\\n\
    -2 ${{FASTQS2[i]}} \\\n\
    -S  results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.classification_HABV \\\n\
    --report-file  results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.report.txt \\\n\
    --threads {MAX_CPU} \\\n\
    {minimum_length_str} {qc_filter_str} {centrigue_extra_args_str}  --mm  \n\n"

    command_centrifuge += f"centrifuge \
    -x {centrifuge_index_str}/f+p \\\n\
    -1 ${{FASTQS1[i]}} \\\n\
    -2 ${{FASTQS2[i]}} \\\n\
    -S  results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.classification_FP \\\n\
    --report-file  results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.report.txt \\\n\
    --threads {MAX_CPU} \\\n\
    {minimum_length_str} {qc_filter_str} {centrigue_extra_args_str}  --mm  \n\n"
    
    
    
    command_centrifuge += f"centrifuge-kreport \\\n\
    -x {centrifuge_index_str}/nt {minimum_length_str} \\\n\
    results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.classification_HABV \\\n\
    results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.classification_FP >> \
        results/{process_name}/{pipeline}/{profiler_type}/${{SAMPLES[i]}}.report"
    

    command_centrifuge += f"done \n\n"

    return command_centrifuge



def run_taxpasta(yaml_dict, list_profilers, process_name, pipeline, list_samples):
    exclusions_yaml_dict = read_yaml('src/config/taxpasta_exclusions.yaml')

    summarise_at_str = f"--summarise_at {yaml_dict[process_name]['taxpasta_config']['summarise_at']}" if check_entry(yaml_dict[process_name]['taxpasta_config'], 'summarise_at') else f"--summarise_at genus"
    output_format_str = f"--output_format {yaml_dict[process_name]['taxpasta_config']['output_format']}" if check_entry(yaml_dict[process_name]['taxpasta_config'], 'output_format') else f"--output_format tsv"
    add_name_str = f"--add_name" if check_entry(yaml_dict[process_name]['taxpasta_config'], 'add_name') else ''
    add_lineage_str = f"--add_lineage" if check_entry(yaml_dict[process_name]['taxpasta_config'], 'add_lineage') else ''
    taxpasta_dir_str = f"--taxonomy {yaml_dict[process_name]['taxpasta_config']['taxpasta_dir']}" if check_entry(yaml_dict[process_name]['taxpasta_config'], 'taxpasta_dir') else f"--taxonomy {DICT_DBS['taxpasta'][0]}"

    taxpasta_extra_args_str = f""
    if check_entry(yaml_dict[process_name]['taxpasta_config'], 'taxpasta_extra_args'):
        val = yaml_dict[process_name]['taxpasta_config']['taxpasta_extra_args']
        taxpasta_extra_args_str = val.strip().strip("'").strip('"').strip()

    command_taxpasta = f"""echo -e "{'^' * 50}\nPROCESS:  {process_name} | TAXPASTA \n{'^' * 50}\n"  \n\n"""
    command_taxpasta += f"""SAMPLES=({' '.join([f'"{i}"' for i in list_samples])}) \n"""

    command_taxpasta += f"""for ((i=0; i<${{#SAMPLES[@]}}; i++)) \n""" 
    command_taxpasta += f"do \n"

    for profiler in list_profilers:
        list_exclusions = exclusions_yaml_dict[profiler]

        if list_exclusions is not None:
            str_exclusions = f'|'.join(map(str, list_exclusions))
            command_taxpasta += f"grep -vE {str_exclusions} results/{process_name}/{pipeline}/{profiler}/${{SAMPLES[i]}}.report > \
                                           results/{process_name}/{pipeline}/{profiler}/${{SAMPLES[i]}}.report.pretaxpasta \n"
        else:
            command_taxpasta += f"cp results/{process_name}/{pipeline}/{profiler}/${{SAMPLES[i]}}.report \
                                     results/{process_name}/{pipeline}/{profiler}/${{SAMPLES[i]}}.report.pretaxpasta    \n"

        command_taxpasta += f"taxpasta standardise -p {profiler} {add_name_str} {add_lineage_str} {summarise_at_str} {taxpasta_dir_str} {output_format_str} \\\n\
                -o results/{process_name}/{pipeline}/{profiler}/${{SAMPLES[i]}}.report.standardised  {taxpasta_extra_args_str} \\\n\
                results/{process_name}/{pipeline}/{profiler}/${{SAMPLES[i]}}.report.pretaxpasta \n\n"

    command_taxpasta += f"done \n\n"

    return command_taxpasta


################################################################################################
# TERCIARY FUNCTIONS
################################################################################################

def write_file(text, project):
    with open(f'work/{project}/command.sh', 'w') as file_out:
        file_out.write(text)


def process_arg_keys(key, val, pipeline, config_arg="--"):
    if val == True:
        return f"{config_arg}{key} \\\n"
    elif val == False:
        return ""
    # FOR TAXPROFILER 1st AND 2nd MAPS!
    elif val in ['1st_map_extra_args', '2nd_map_extra_args']:
        return f"""{val.strip().strip("'").strip('"').strip()} \\\n"""
    elif key == 'genome_fasta':
        return f"{config_arg}fasta {val} \\\n"
    # elif (pipeline == 'circdna') & (key == 'bwa_index'):  # In circna is --bwa and path to dir, and in circdna is --bwa_index and path to .bwt file
    #     return f"{config_arg}{key} {val}/genome.bwt \\\n"
    elif (pipeline == 'circdna') & (key == 'reference_build') & (val == 'GRCm38'):
        return f"{config_arg}{key} mm10 \\\n"
    elif 'gtf_corrected'in key: 
        return f"{config_arg}gtf {val} \\\n"
    elif key in ['aa_data_repo', 'mirbase_gff', 'mirgene_gff', 'mirbase_mature', 'mirbase_hairpin', 
                 'mirgenedb_mature', 'mirgenedb_hairpin']:
                    # The path for aa_data_repo has to be absolute https://github.com/nf-core/circdna/issues/69
        return f"{config_arg}{key} $(pwd)/{val} \\\n"
    elif 'mirbase_' in key:
        if key == 'mirbase_gff':
            return f"{config_arg}mirna_gtf {val} \\\n"
        else:
            return f"{config_arg}{key.replace('mirbase_', '')} {val} \\\n"
    
    # For the rest of arguments
    else:
        if 'extra_' in key:
            val_strip = val.strip().strip("'").strip('"').strip()
            val = f'"{val_strip}"'

        return f"{config_arg}{key} {val} \\\n"
    

def return_list_fastqs(df_samplesheet):
    list_samples = df_samplesheet.drop_duplicates('sample')['sample'].values.tolist()
    list_fastqs_1 = []
    list_fastqs_2 = []

    for sample in list_samples: # BUG: Este código puede dar fallo si no hay fastq_2 en el samplesheet.csv!
        df_sample = df_samplesheet[df_samplesheet['sample'] == sample]
        list_fastqs_1.append(' '.join([f"'{i}'" for i in df_sample['fastq_1'].values.tolist()]))
        list_fastqs_2.append(' '.join([f"'{i}'" for i in df_sample['fastq_2'].values.tolist()]))

    return list_samples, list_fastqs_1, list_fastqs_2