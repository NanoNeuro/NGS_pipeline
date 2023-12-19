import logging
import subprocess
import psutil
import yaml



logger = logging.getLogger()



################################################################################################
# GLOBAL VARIABLES
################################################################################################
# These pipeline have to be in the correct order of execution!
VALID_PIPELINES = ['rnaseq', 'scrnaseq', 'smrnaseq', 'circrna', 'circdna', 'taxprofiler']

NF_CONFIG_PIPELINES = ['rnaseq', 'scrnaseq', 'smrnaseq', 'circrna', 'circdna']
DEV_PIPELINES = ['circrna']

VALID_ORGANISMS = ['human', 'mouse']

DEFAULT_RNASEQ_ALIGNER = 'star_salmon'

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


    # With the ordered dictionary, we are first going to check and fill the general information
    parse_nfcore_config(yaml_dict)


    # Create the list of samplesheets and complete the empty input paths
    parse_input_samplesheet(yaml_dict, project, samplesheet, list_samplesheets)


    # Then, we are going to check the database to download and append that to the list of databases to download
    yaml_dict, list_dbs_to_download = parse_database_arguments(yaml_dict, list_dbs_to_download)



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
            yaml_dict[process_name]['general_config']['outdir'] = f"results/{project}/"
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
                logger.warning(f"The number of CPUs set ({yaml_dict[process_name]['general_config']['max_cpus']}) 
                               is larger than the maximum allowed number ({get_available_cpus()}).")
            
        if yaml_dict[process_name]['general_config']['max_memory'] == False:
            yaml_dict[process_name]['general_config']['max_memory'] = int(MAX_MEMORY_PER * get_RAM_GB())
        else:
            if yaml_dict[process_name]['general_config']['max_cpus'] > get_available_cpus():
                logger.warning(f"The number of CPUs set ({yaml_dict[process_name]['general_config']['max_cpus']}) 
                               is larger than the maximum allowed number ({get_available_cpus()}).")

        if yaml_dict[process_name]['general_config']['max_time'] == False:
            yaml_dict[process_name]['general_config']['max_time'] = 500




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
            # In case nothing is set in the YAML file
            if 'nextflow_config' not in yaml_dict[process_name]:
                yaml_dict[process_name]['nextflow_config'] = {'r': False, 'resume': False, 'profile': False}

            # Check r, resume and pipeline
            nextflow_config = yaml_dict[process_name]['nextflow_config']
            pipeline = yaml_dict[process_name]['general_config']['pipeline']

            r = nextflow_config['r'] if 'r' in nextflow_config else False
            resume = nextflow_config['resume'] if 'resume' in nextflow_config else False
            profile = nextflow_config['profile'] if 'profile' in nextflow_config else False

            if r == False:
                # Check the latest version of the pipeline
                if pipeline in DEV_PIPELINES:
                    r = 'dev'
                else:
                    r = retrieve_r_tag(pipeline)

                nextflow_config['r'] = r

            if profile == False:
                nextflow_config['profile'] = 'docker'
            
            if resume == False:
                nextflow_config['resume'] = True
            
            ### NF-CORE CONFIG (ALL PROCESSES)
            for element in ['outdir', 'genome',]:
                if element not in yaml_dict[process_name]['nfcore_config']:
                    yaml_dict[process_name]['nfcore_config'][element] = False
            
            if yaml_dict[process_name]['nfcore_config']['outdir'] == False:
                yaml_dict[process_name]['nfcore_config']['outdir'] = f"{yaml_dict[process_name]['nfcore_config']['outdir']}/{pipeline}/"


            dict_genomes = {'human': 'GRCh38', 'mouse': 'GRCm38'}
            if yaml_dict[process_name]['nfcore_config']['genome'] == False:
                organism = yaml_dict[process_name]['nfcore_config']['organism']
                yaml_dict[process_name]['nfcore_config']['genome'] = dict_genomes[organism]



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
    - If a database file is absent, the path to the database will be appended. 
    - If the default database file is present but the database vale is absent in the config file, the file name will be added, 
      but the database will not be downloaded.
    - If the value of the database is already filled, we will check that the file/directory exists.
    """

    for process_name in yaml_dict.keys():
        pipeline = yaml_dict[process_name]['general_config']['pipeline']

        # TODO: PENSAR COMO HACER EL ARRANGEMENT DE LAS COSAS (COMO DETERMINAR QUÉ BAJAR)
        # SEGUN ALIGNER
        # Hay que mirar si la base de datos es descargable/construible -> si no dejarla en false
        # Y poner un info!
        if pipeline in NF_CONFIG_PIPELINES:
            match pipeline:
                case 'rnaseq':

                    # Check the aligner used to set the databases that need to be downloaded
                    # Pseudoaligners can be used as the only option. In that case, we will use only
                    # the pseudoaligner option

                    is_aligner = check_entry(yaml_dict[process_name]['general_config'], 'aligner')
                    is_pseudoaligner = check_entry(yaml_dict[process_name]['general_config'], 'pseudo_aligner')

                    if (not is_aligner) & (not is_pseudoaligner):
                        yaml_dict[process_name]['general_config']['aligner'] = DEFAULT_RNASEQ_ALIGNER
                        aligner = yaml_dict[process_name]['general_config']['aligner']
                        pseudoaligner = None
                    elif (is_aligner) & (not is_pseudoaligner):
                        aligner = yaml_dict[process_name]['general_config']['aligner']
                        pseudoaligner = None
                    elif (not is_aligner) & (is_pseudoaligner):
                        aligner = None
                        pseudoaligner = yaml_dict[process_name]['general_config']['pseudo_aligner']
                    else :
                        aligner = yaml_dict[process_name]['general_config']['aligner']
                        pseudoaligner = yaml_dict[process_name]['general_config']['pseudo_aligner']

                    match aligner:
                        case [] 

                    if yaml_dict[nfcore_config]['general_config']
                    list_dbs_to_check = ['star_index', 'salm']
                case 'scrnaseq':
                
                case 'smrnaseq':
                
                case 'circrna':

                case 'circdna':

                case _:
                    pass

    check_entry()
    



def parse_process(yaml_dict_element):
    ...




################################################################################################
# TERCIARY FUNCTIONS
################################################################################################
def custom_sort(item):
    return VALID_PIPELINES.index(item[1])


def sort_pipelines(yaml_dict):
    """
    Sort the entries (processes) of yaml_dict so that secondary processes dependent on primary ones (like differential abundance)
    are run first
    """
    list_pipelines = [yaml_dict[i]['general_config']['pipeline'] for i in yaml_dict.keys()]
    sorted_dict = sorted(dict(zip(yaml_dict.keys(), list_pipelines)).items(), key=custom_sort)
    sorted_pipelines = [i[0] for i in sorted_dict]
    
    ordered_yaml_dict = {k: yaml_dict[k] for k in sorted_pipelines}
    return ordered_yaml_dict


def retrieve_r_tag(pipeline):
    # Run the shell command and capture its output
    command = f'curl https://api.github.com/repos/nf-core/{pipeline}/tags | grep -oP \'"name": "\\K([^\"]*)\' | head -n 1'
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
    available_cpus = psutil.cpu_count(logical=False)  # logical=False returns the number of physical CPUs
    return available_cpus


def check_entry(dict_eval, dict_key):
    if dict_key not in dict_eval:
        return False

    if dict_eval[dict_key] == False:
        return False
    else:
        return True