import logging
import yaml
import subprocess


logger = logging.getLogger()

# These pipeline have to be in the correct order of execution!
VALID_PIPELINES = ['rnaseq', 'scrnaseq', 'smrnaseq', 'circrna', 'circdna', 'taxprofiler']
NF_CONFIG_PIPELINES = ['rnaseq', 'scrnaseq', 'smrnaseq', 'circrna', 'circdna']
DEV_PIPELINES = ['circrna']
VALID_ORGANISMS = ['human', 'mouse']

# TODO: SORTEAR LOS PIPELINES EN UN ORDEN ESPECÍFICO (https://chat.openai.com/c/5253f50d-7810-4f0e-bc16-c7e4a18a02e5)

def parse_yaml_file(yaml_path):
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

    for process_name, process_info in yaml_dict.items():
        logging.info(f'Parsing YAML - general config (process {process_name}).')
        # We first are going to check that the pipeline exists and the organism is human or mouse.
        parse_general_info(process_name, process_info)


    # Now we are going to order the dictionary of processes. This is because some pipelines may require results from previous 
    # pipelines
    sort_pipelines(yaml_dict)


    # With the ordered dictionary, we are first going to check and fill the general information
    parse_general_arguments(yaml_dict, list_samplesheets)

    
    # Then, we are going to check the database to download and append that to the list of databases to download
    yaml_dict, list_dbs_to_download = parse_database_arguments(yaml_dict, list_dbs_to_download)



    return yaml_dict, list_samplesheets, list_dbs_to_download




def parse_general_info(process_name, process_info):
    logging.info(f'Parsing information of process {process_name}.')
        
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


def parse_general_arguments(yaml_dict, list_samplesheets):
    """
    Check and fill general arguments including:
    - r version on nf-core dependent processes
    - profile of nf-core dependent processes
    - input samplesheet and and output directory
    - genome / species
    """
    
    for process_name in yaml_dict.keys():
        pipeline = yaml_dict[process_name]['general_config']['pipeline']

        if pipeline in NF_CONFIG_PIPELINES:
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

                nextflow_config = yaml_dict[process_name]['nextflow_config']
            if profile == False:
                nextflow_config['profile'] = 'docker'


    return yaml_dict



def retrieve_r_tag(pipeline):
    # Run the shell command and capture its output
    command = f'curl https://api.github.com/repos/nf-core/{pipeline}/tags | grep -oP \'"name": "\\K([^\"]*)\' | head -n 1'
    output = subprocess.check_output(command, shell=True, text=True)

    # Store the output in a variable
    version = output.strip()
    return version
    


def parse_database_arguments(yaml_dict, list_dbs_to_download):
    """
    Parse the YAML to check for undownloaded database files. In that case, the path to the database will be appended and 
    the database to be downloaded will be added to the list of to-be-downloaded databases.
    """
    ...


def parse_process(yaml_dict_element):
    ...

