import logging
import os
import psutil
import subprocess
import yaml

from .database_download import dict_dbs as DICT_DBS


logger = logging.getLogger()



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


################################################################################################
# PRIMARY FUNCTIONS
################################################################################################
def parse_yaml_file(yaml_path, project, samplesheet):
    """
    This function parses the yaml file and executes several tasks:
        - Create the necessary samplesheets for each process
        - Create the list of databases to download
        - Create the list of commands that will be executed in order

    The yaml files have to be consistent will be checked to make sure that the structure is correct.

    Args:
        yaml_file (_type_): _description_

    Returns:
        _type_: _description_
    """

    list_samplesheets, list_dbs_to_download = [], []

    yaml_dict = read_yaml(yaml_path)
    with open(yaml_path, 'r') as file:
        try:
            yaml_dict = yaml.safe_load(file)
        except Exception as e: 
            logger.error('YAML configuration is not valid.', exc_info=True)
            raise
            
    try:
        assert yaml_dict is not None
        assert len(yaml_dict) > 0
    except AssertionError:
        logger.error('There are no active processes in the yaml file.', exc_info=True); raise


    logging.info(f'Parsing YAML.')

    # Checks of general arguments
    parse_general_config(yaml_dict, project)


    # Now we are going to order the dictionary of processes. This is because some pipelines may require results from previous 
    # pipelines
    sort_pipelines(yaml_dict)


    # With the ordered dictionary, we are first going  to check and fill the general information
    parse_nfcore_config(yaml_dict)


    # There are specific arguments that should be checked for taxprofiler, so we are going to process them separately
    parse_taxprofiler_config(yaml_dict)


    # Create the list of samplesheets and complete the empty input paths
    parse_input_samplesheet(yaml_dict, project, samplesheet, list_samplesheets)


    # Then, we are going to check the database to download and append that to the list of databases to download
    parse_database_arguments(yaml_dict, list_dbs_to_download)



    return yaml_dict, list_samplesheets, list_dbs_to_download


################################################################################################
# SECONDARY FUNCTIONS
################################################################################################
def parse_general_config(yaml_dict, project):
    for process_name, process_info in yaml_dict.items():
        logging.info(f'Parsing YAML - general config (process {process_name}).')
        # We first are going to check that the pipeline exists and the organism is human or mouse.
        
        if 'general_config' not in process_info:
            e = f'Process {process_name} - general_config must be a element of the process in the yaml.'
            logger.error(e, exc_info=True); raise AssertionError (e)

        try:
            pipeline, organism = process_info['general_config']['pipeline'], process_info['general_config']['organism']
        except:
            e = f'Process {process_name} - general_config must have pipeline and organism elements included.'
            logger.error(e, exc_info=True); raise AssertionError (e)

        if pipeline not in VALID_PIPELINES:
            e = f'Pipeline {pipeline} in process {process_name} is not a valid pipeline {VALID_PIPELINES}.'
            logger.error(e, exc_info=True); raise AssertionError (e)

        if organism not in VALID_ORGANISMS:
            e = f'Organism {organism} in process {process_name} is not a valid organism {VALID_ORGANISMS}.'
            logger.error(e, exc_info=True); raise AssertionError (e)


        # Check for other necessaty variables
        for element in ['outdir', 'max_cpus', 'max_memory', 'max_time']:
            if element not in yaml_dict[process_name]['general_config']:
                yaml_dict[process_name]['general_config'][element] = False

        if yaml_dict[process_name]['general_config']['outdir'] == False:
            # TODO: MKDIR
            yaml_dict[process_name]['general_config']['outdir'] = f"results/{project}/{process_name}"
        else:
            if 'nfcore_config' in yaml_dict:
                if 'outdir' in yaml_dict['nfcore_config']:
                    if yaml_dict['nfcore_config']['outdir'] != False:
                        logger.warning(f"Config mismatch in process {process_name}: outdir in general config ({yaml_dict['general_config']['outdir']}) \
                                     does not match outfir in nfcore_config ({yaml_dict['nfcore_config']['outdir']}). Outdir in nfcore_config will be selected.")

        # TODO: CHECK THAT THE LIMITS ARE CORRECT; AND ELSE ADD A WARNING
        if yaml_dict[process_name]['general_config']['max_cpus'] == False:
            yaml_dict[process_name]['general_config']['max_cpus'] = int(MAX_CPU_PER * get_available_cpus())
        else:
            if yaml_dict[process_name]['general_config']['max_cpus'] > get_available_cpus():
                logger.warning(f"The number of CPUs set ({yaml_dict[process_name]['general_config']['max_cpus']}) \
                               is larger than the maximum allowed number ({get_available_cpus()}).")
            
        if yaml_dict[process_name]['general_config']['max_memory'] == False:
            yaml_dict[process_name]['general_config']['max_memory'] = int(MAX_MEMORY_PER * get_RAM_GB())
        else:
            if yaml_dict[process_name]['general_config']['max_cpus'] > get_available_cpus():
                logger.warning(f"The number of CPUs set ({yaml_dict[process_name]['general_config']['max_cpus']}) \
                               is larger than the maximum allowed number ({get_available_cpus()}).")

        if yaml_dict[process_name]['general_config']['max_time'] == False:
            yaml_dict[process_name]['general_config']['max_time'] = 500

        # Check if other mayor yaml configs are available
        #   - nextflow_config, nfcore_config in nf_core related configs
        #   - profiler_config, host_mapping_config, in taxprofiler
        
        pipeline = yaml_dict[process_name]['general_config']['pipeline']

        if pipeline in NF_CONFIG_PIPELINES:
            configs_to_check = ['nextflow_config', 'nfcore_config']
        elif pipeline == 'taxprofiler':
            configs_to_check = ['profiler_config', 'host_mapping_config']
            # SPECIFIC PROFILER CONFIGS ARE CHECKED in parse_taxprofiler_config() !
        else:
            configs_to_check = []
        
        for config in configs_to_check:
            if config not in yaml_dict[process_name]:
                yaml_dict[process_name][config] = {}
        



def parse_nfcore_config(yaml_dict):
    """
    Check and fill arguments including:
    - r version on nf-core dependent processes
    - profile of nf-core dependent processes
    - output directory
    - genome / species
    - cpu, memory, time
    """
    
    for process_name in yaml_dict.keys():
        pipeline = yaml_dict[process_name]['general_config']['pipeline']

        if pipeline in NF_CONFIG_PIPELINES:
            ### NEXTFLOW CONFIG
            # Check r, resume and pipeline # TODO: REHACER ESTO!!!
            nextflow_config = yaml_dict[process_name]['nextflow_config']
            pipeline = yaml_dict[process_name]['general_config']['pipeline']

            if not check_entry(nextflow_config, 'r'): 
                if pipeline in DEV_PIPELINES:
                    nextflow_config['r'] = 'dev'
                else:
                    nextflow_config['r'] = retrieve_r_tag(pipeline)
                        
            if not check_entry(nextflow_config, 'resume'): nextflow_config['resume'] = True
            if not check_entry(nextflow_config, 'profile'): nextflow_config['profile'] = 'docker'

            ### NF-CORE CONFIG (ALL PROCESSES)
            for element in ['outdir', 'genome',]:
                if element not in yaml_dict[process_name]['nfcore_config']:
                    yaml_dict[process_name]['nfcore_config'][element] = False
            
            if yaml_dict[process_name]['nfcore_config']['outdir'] == False:
                yaml_dict[process_name]['nfcore_config']['outdir'] = f"results/{pipeline}/{process_name}"


            
            if yaml_dict[process_name]['nfcore_config']['genome'] == False:
                organism = yaml_dict[process_name]['general_config']['organism']
                yaml_dict[process_name]['nfcore_config']['genome'] = DICT_GENOMES[organism]



def parse_taxprofiler_config(yaml_dict):
    """
    Processes arguments specific of taxprofiler-derived processes. In this function we check:
    1) The profilers chosen. If any profiler is not within the list of available ones, we will raise an error.
    2) The host mapping configuration. We can accept 4 scenarios:
        - No mapping is performed. This option is not recommended, so it will raise a warning. 
        - 1st mapping (GRCh38) is performed, and not 2nd one.
        - 2nd mapping is performed (CHM13) and not 1st one. This option is not recommended, so it will raise a warning. 
        - Both 1st and 2nd mapping are performed.
    3) Check that the dictionary of profilers' configuration exists.
    4) Taxpasta argument checks. 
    """

    for process_name in yaml_dict.keys():
        pipeline = yaml_dict[process_name]['general_config']['pipeline']

        if pipeline == 'taxprofiler':
            # 1) Profiler checking
            is_profilers = check_entry(yaml_dict[process_name]['profiler_config'], 'profilers')

            if not is_profilers:
                yaml_dict[process_name]['profiler_config']['profilers'] = ','.join(VALID_PROFILERS)
            
            # We will reset the list of profilers by stripping spaces
            list_profilers_srt = yaml_dict[process_name]['profiler_config']['profilers']
            list_profilers = [variable.strip() for variable in list_profilers_srt.split(',')]

            yaml_dict[process_name]['profiler_config']['profilers'] = ','.join(list_profilers)
            
            for profiler in list_profilers:
                if profiler not in VALID_PROFILERS:
                    err = f"Profiler {profiler} not in the list of valid profilers: {', '.join(VALID_PROFILERS)}"
                    logger.error(err); raise AssertionError(err)

            # 2) Host mapping configuration
            is_map1 = check_entry(yaml_dict[process_name]['host_mapping_config'], '1st_map')
            if not is_map1:
                yaml_dict[process_name]['host_mapping_config']['1st_map'] = True

            is_map2 = check_entry(yaml_dict[process_name]['host_mapping_config'], '2nd_map')
            if not is_map2:
                yaml_dict[process_name]['host_mapping_config']['2nd_map'] = True
            
            map1, map2 = yaml_dict[process_name]['host_mapping_config']['1st_map'], yaml_dict[process_name]['host_mapping_config']['2nd_map']

            if ((not map1) & (not map2)) | ((not map1) & (map2)):
                logger.warning('We recommend setting 1st_map to true to remove host reads!')

            # 3) Configuration of profilers
            ## Check that the config dictionary exists
            for profiler in list_profilers:
                if f'{profiler}_config' not in yaml_dict[process_name]:
                    yaml_dict[process_name][f'{profiler}_config'] = {}
            
            # 4) Taxpasta argument check
            if f'taxpasta_config' not in yaml_dict[process_name]:
                    yaml_dict[process_name][f'taxpasta_config'] = {'summarise_at': 'genus', 'add_name': True, 'add_lineage': True}


def parse_input_samplesheet(yaml_dict, project, samplesheet, list_samplesheets):
    """
    Process input samplesheets. All projects must be contained within the master samplesheet, so any input argument
    than contains a samplesheet will raise an error. With this function we will return a list of tuples like that:
    [(process_name, pipeline, path_to_process_csv)]. Also, the input variables will be writen for the necessary 
    """

    for process_name in yaml_dict.keys():
        pipeline = yaml_dict[process_name]['general_config']['pipeline']
        
        path_to_process_csv = f"work/{project}/samplesheets/{process_name}.csv"

        if pipeline in NF_CONFIG_PIPELINES:
            if 'input' not in yaml_dict[process_name]['nfcore_config']:
                yaml_dict[process_name]['nfcore_config']['input'] = False

            yaml_config_input_file = yaml_dict[process_name]['nfcore_config']['input']
            if yaml_config_input_file == False:
                    yaml_dict[process_name]['nfcore_config']['input'] = path_to_process_csv
                    list_samplesheets.append((project, pipeline, path_to_process_csv))
            else:
                if yaml_config_input_file != path_to_process_csv:
                    logger.warning(f"Input file in yaml config ({yaml_config_input_file}) differs from the expected path \
                                   ({path_to_process_csv}). We recommend to add the files in the main spreadsheet \
                                    (/projects/{project}/{samplesheet}), and set this \
                                    configuration as false in the config yaml file.", exc_info=True)
        else:
            list_samplesheets.append((project, pipeline, path_to_process_csv))


def parse_database_arguments(yaml_dict, list_dbs_to_download):
    """
    Parse the YAML to check for undownloaded database files. 

    If the process runs on nf-core:
    - If the database path is present, we don't add the element to the list of databases to download.
        - We check that the path is correct, and raise an error in case it is not.
    - If the database path is absent and the element if downloadable (in the list of downloadable databases), it will be downloaded.
    - If the database path is absent but the element is not downloadable, it will be kept empty (this is for nf-core-related databases).

    For taxprofiler:
    - Check if any of the arguments of index download is false. If it is not, it will check the path of the index.
    """

    for process_name in yaml_dict.keys():
        pipeline = yaml_dict[process_name]['general_config']['pipeline']
        organism = yaml_dict[process_name]['general_config']['organism']

        # TODO: PENSAR COMO HACER EL ARRANGEMENT DE LAS COSAS (COMO DETERMINAR QUÉ BAJAR)
        # SEGUN ALIGNER
        # Hay que mirar si la base de datos es descargable/construible -> si no dejarla en false
        # Y poner un info!
        if pipeline in NF_CONFIG_PIPELINES:
            fillable_args = []
            match pipeline:
                case 'rnaseq':

                    # Check the aligner used to set the databases that need to be downloaded
                    # Pseudoaligners can be used as the only option. In that case, we will use only
                    # the pseudoaligner option

                    is_aligner = check_entry(yaml_dict[process_name]['nfcore_config'], 'aligner')
                    is_pseudoaligner = check_entry(yaml_dict[process_name]['nfcore_config'], 'pseudo_aligner')

                    if (not is_aligner) & (not is_pseudoaligner):
                        yaml_dict[process_name]['nfcore_config']['aligner'] = DEFAULT_RNASEQ_ALIGNER
                        aligner = yaml_dict[process_name]['nfcore_config']['aligner']
                        pseudoaligner = False
                    elif (is_aligner) & (not is_pseudoaligner):
                        aligner = yaml_dict[process_name]['nfcore_config']['aligner']
                        pseudoaligner = False
                    elif (not is_aligner) & (is_pseudoaligner):
                        aligner = False
                        pseudoaligner = yaml_dict[process_name]['nfcore_config']['pseudo_aligner']
                    else :
                        aligner = yaml_dict[process_name]['nfcore_config']['aligner']
                        pseudoaligner = yaml_dict[process_name]['nfcore_config']['pseudo_aligner']

                    fillable_args += ['fasta', 'gtf', 'gene_bed']

                    match aligner:
                        # TODO: checkear si en configs sin salmon hace falta el indice de salmon para el paso de infer_strandness
                        case 'star_salmon':
                            fillable_args += ['star_index', 'salmon_index']
                        case 'star_rsem':
                            fillable_args += ['star_index', 'rsem_index']
                        case 'hisat2':
                            fillable_args += ['hisat2_index']
                        case False:
                            pass
                        case _:
                            err = "rnaseq pipeline yaml currently supports the following aligners: \
                                   star_salmon, star_rsem, hisat2."
                            logger.error(err); raise AssertionError(err)

                    match pseudoaligner:
                        case 'kallisto':
                            fillable_args += ['kallisto_index']
                        case 'salmon': 
                            fillable_args += ['salmon_index']
                        case False:
                            pass
                        case _:
                            err = "rnaseq pipeline yaml currently supports the following pseudoaligners: \
                                   kallisto, salmon."
                            logger.error(err); raise AssertionError(err)


                case 'scrnaseq':
                    is_aligner = check_entry(yaml_dict[process_name]['nfcore_config'], 'aligner')
                    if (not is_aligner):
                        yaml_dict[process_name]['nfcore_config']['aligner'] = DEFAULT_SCRNASEQ_ALIGNER

                    fillable_args += ['fasta', 'transcript_fasta', 'gtf',]
                    match yaml_dict[process_name]['nfcore_config']['aligner']:
                        # TODO: checkear si en configs sin salmon hace falta el indice de salmon para el paso de infer_strandness
                        case 'alevin':
                            fillable_args += ['salmon_index', 'txp2gene']
                        case 'kallisto':
                            # TODO: check how to create kallisto gene map
                            fillable_args += ['kallisto_index', 'kallisto_gene_map']
                        case 'star':
                            fillable_args += ['star_index']
                        case 'cellranger':
                            fillable_args += ['cellranger_index']
                        case 'universc':
                            fillable_args += ['universc_index']
                        case _:
                            err = "scrnaseq pipeline yaml currently supports the following aligners: \
                                   alevin, kallisto, star, cellranger."
                            logger.error(err); raise AssertionError(err)


                case 'smrnaseq':
                    fillable_args += ['fasta', 'bowtie_index']
                    is_mirgenedb = check_entry(yaml_dict[process_name]['nfcore_config'], 'mirgenedb')
                    if (not is_mirgenedb):
                        yaml_dict[process_name]['nfcore_config']['mirgenedb'] = DEFAULT_MIRGENEDB

                    if yaml_dict[process_name]['nfcore_config']['mirgenedb'] == True:
                        fillable_args += ['mirgenedb_gff', 'mirgenedb_mature', 'mirgene_hairpin']
                    else:
                        fillable_args += ['mirna_gtf', 'mature', 'hairpin']


                case 'circrna':
                    # TODO: CHECKEAR QUE INDICES SE EMPLEAN PARA CADA TOOL
                    # TODO: CHECKEAR QUE LOS INDICES EN FALSE NO SE USAN: p.ej si no uso segemehl, que no me cree el index el programa
                    fillable_args += ['fasta', 'gtf', 'star', 'bowtie', 'bowtie2', 'bwa', 'hisat2',]

                    circrna_module = check_entry(yaml_dict[process_name]['nfcore_config'], 'module')
                    if (not circrna_module):
                        yaml_dict[process_name]['nfcore_config']['module'] = DEFAULT_CIRCRNA_MODULE

                    if 'mirna_prediction' in yaml_dict[process_name]['nfcore_config']['module']:
                        fillable_args += ['mature']

                    circrna_tool = check_entry(yaml_dict[process_name]['nfcore_config'], 'tool')
                    if (not circrna_tool):
                        yaml_dict[process_name]['nfcore_config']['tool'] = DEFAULT_CIRCRNA_TOOL

                    if 'segemehl' in yaml_dict[process_name]['nfcore_config']['tool']:
                        fillable_args += ['segemehl']


                case 'circdna':
                    fillable_args += ['fasta', 'bwa']

                    cirdrna_tool = check_entry(yaml_dict[process_name]['nfcore_config'], 'circle_identifier')
                    if (not cirdrna_tool):
                        yaml_dict[process_name]['nfcore_config']['circle_identifier'] = DEFAULT_CIRCDNA_TOOL
                    
                    if 'ampliconarchitect' in yaml_dict[process_name]['nfcore_config']['circle_identifier']:
                        fillable_args += ['aa_data_repo']
                        yaml_dict[process_name]['nfcore_config']['reference_build'] = \
                            yaml_dict[process_name]['nfcore_config']['genome']
                        

                case _:
                    pass
                
            # Check each of the fillable args to see if it needs to be downloaded or not
            for arg in fillable_args:
                is_arg = check_entry(yaml_dict[process_name]['nfcore_config'], arg)

                if is_arg:  # The path exists
                    if not os.path.exists(yaml_dict[process_name]['nfcore_config'][arg]):
                        err = f"Path {yaml_dict[process_name]['nfcore_config'][arg]} does not exist."
                        logger.error(err); raise AssertionError(err)
                
                else:  # There is no path
                    if arg in ['star', 'bowtie', 'bowtie2', 'bwa', 'hisat2', 'segemehl']:
                        arg_db = arg + '_index_' + DICT_GENOMES[organism]
                    elif arg in ['fasta', 'transcript_fasta', 'gtf', 'bed', 'star_index', 'bowtie_index', 'bowtie2_index', \
                                 'bwa_index', 'hisat_index', 'salmon_index', 'rsem_index', 'kallisto_index', 'kallisto_gene_map', \
                                    'cellranger_index', 'universc_index']:
                        arg_db = arg + '_' + DICT_GENOMES[organism]
                    elif arg == 'gene_bed':
                        arg_db = 'bed_' + DICT_GENOMES[organism]
                    elif arg in ['mirna_gtf', 'mature', 'hairpin']:
                        arg_db = 'mirbase_' + arg.replace('mirna_', '')
                        if arg_db == 'mirbase_gtf':
                            arg_db = arg_db + '_' + DICT_GENOMES[organism]
                    elif arg == 'mirgenedb_gff':
                        arg_db = arg + '_' + DICT_GENOMES[organism]
                    elif arg == 'txp2gene':
                        arg_db = arg + '_' + DICT_GENOMES[organism]
                    else:
                        err = f'Argument {arg} not available.'
                        logger.error(err); raise AssertionError(err)
                    

                    # Check that the element is downloadable
                    if arg_db in DICT_DBS.keys():
                        arg_path, _ = DICT_DBS[arg_db]
                        yaml_dict[process_name]['nfcore_config'][arg] = arg_path

                        # Append the arg to the list of downloadable datasets
                        if arg_db not in list_dbs_to_download:
                            list_dbs_to_download.append(arg_db)
                    
                    else:
                        logger.info(f'Argument {arg} from process {process_name} is not downloadable. \
                                      It will be downloadable by nf-core.')

                    
        elif pipeline == 'taxprofiler': # pipeline is taxprofiler
            fillable_args = []
            # databases related to 1st and 2nd maps: 
            if yaml_dict[process_name]['host_mapping_config']['1st_map']:
                for db in ['fasta', 'gtf', 'star_index', 'salmon_index', 'rsem_index']:
                    fillable_args.append(db + '_' + DICT_GENOMES[organism])

            if yaml_dict[process_name]['host_mapping_config']['2nd_map']:
                fillable_args += ['fasta_CHM13', 'bowtie2_index_CHM13']

            list_profilers_str = yaml_dict[process_name]['profiler_config']['profilers']
            list_profilers = list_profilers_str.split(',')

            for profiler in list_profilers:
                if not check_entry(yaml_dict[process_name][f'{profiler}_config'], f'{profiler}_index'):
                    fillable_args +=  [profiler]

            list_dbs_to_download += [i for i in fillable_args if i not in list_dbs_to_download]

    # We sort the values in case there are some elements that have to be downloaded before others.
    # For example, to download txp2gene we need to build the index of salmon first.        
    list_dbs_to_download = sorted(list_dbs_to_download, key=custom_sort_dict_dbs)



################################################################################################
# TERCIARY FUNCTIONS
################################################################################################
def read_yaml(yaml_path):
    with open(yaml_path, 'r') as file:
        try:
            yaml_dict = yaml.safe_load(file)
        except Exception as e: 
            logger.error('YAML configuration is not valid.', exc_info=True)
            raise



def write_yaml_file(yaml_dict, yaml_path):
    with open(yaml_path, 'w') as outfile:
        yaml.dump(yaml_dict, outfile, default_flow_style=False, sort_keys=False)    



def custom_sort_valid_pipelines(item):
    return VALID_PIPELINES.index(item[1])



def custom_sort_dict_dbs(item):
    list_dbs = list(DICT_DBS.keys())
    return list_dbs.index(item)



def sort_pipelines(yaml_dict):
    """
    Sort the entries (processes) of yaml_dict so that secondary processes dependent on primary ones (like differential abundance)
    are run first
    """
    list_pipelines = [yaml_dict[i]['general_config']['pipeline'] for i in yaml_dict.keys()]
    sorted_dict = sorted(dict(zip(yaml_dict.keys(), list_pipelines)).items(), key=custom_sort_valid_pipelines)
    sorted_pipelines = [i[0] for i in sorted_dict]
    
    ordered_yaml_dict = {k: yaml_dict[k] for k in sorted_pipelines}
    return ordered_yaml_dict



def retrieve_r_tag(pipeline):
    # Run the shell command and capture its output
    command = f'curl -s --url https://api.github.com/repos/nf-core/{pipeline}/tags --header "Authorization: Bearer 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" \
                | grep -oP \'"name": "\\K([^\"]*)\' | head -n 1'
    output = subprocess.check_output(command, shell=True, text=True)

    # Store the output in a variable
    version = output.strip()
    return version
    


def get_RAM_GB():
    virtual_memory = psutil.virtual_memory()
    # Calculate and return the available RAM in gigabytes
    available_ram_gb = virtual_memory.available / (1024 ** 3)
    return available_ram_gb



def get_available_cpus():
    # Get the number of available CPUs
    available_cpus = psutil.cpu_count(logical=True)
    return available_cpus



def check_entry(dict_eval, dict_key):
    if dict_key not in dict_eval:
        return False

    if dict_eval[dict_key] == False:
        return False
    else:
        return True