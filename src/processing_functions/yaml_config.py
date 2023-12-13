import yaml


def process_yaml_file(yaml_path):
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
        yaml_file = yaml.safe_load(file)

    print(yaml_file)


    return list_samplesheets, list_dbs_to_download, list_pipeline_commands