import logging
import yaml


logger = logging.getLogger()


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

    list_samplesheets, list_dbs_to_download, list_pipeline_commands = [], [], []

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
        # We first are going to check that the pipeline exists and the organism is human or mouse.
        parse_general_info(process_name, process_info)



def parse_process(yaml_dict_element):
    ...


def parse_general_info(process_name, process_info):
    print('HOLA')
    logging.info(f'Parsing information of process {process_name}.')
        
    if 'general_config' not in process_info:
        e = 'general_config must be a element of the process in the yaml.'
        logger.error(e, exc_info=True); raise AssertionError (e)

    try:
        pipeline, organism = process_info['general_config']['pipeline'], process_info['general_config']['organism']
    except:
        e = 'general_config must have pipeline and organism elements included.'
        logger.error(e, exc_info=True); raise AssertionError (e)


    valid_pipelines = ['rnaseq', 'circrna', 'cirdna', 'smrnaseq', 'scrnaseq', 'taxprofiler']



    # TODO: SORTEAR LOS PIPELINES EN UN ORDEN ESPECÍFICO (https://chat.openai.com/c/5253f50d-7810-4f0e-bc16-c7e4a18a02e5)



def parse_profiler_config():
    ...