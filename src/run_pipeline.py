import argparse
import logging
import os

from processing_functions.logging import setup_logger
from processing_functions.generic_process import create_dirs
from processing_functions.yaml_config import parse_yaml_file, write_yaml_file
from processing_functions.command_creation import create_command_sh_file
from processing_functions.samplesheet_creation import create_samplesheets
from processing_functions.database_download import database_download
from processing_functions.versions import VERSION

# TODO: MIRAR LAS COLUMNAS OPCIONALES Y OBlIGATORIAS DE CADA ELEMENTO


def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description='Python script to run combined NGS pipelines.')
    
    # Add argument for log level (required)
    parser.add_argument('--project', required=True,  help='Name of the project. It has to be the same as the \
                        folder name located in the projects folder.')
    parser.add_argument('--verbose', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='WARNING',
                        help='Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)')
    parser.add_argument('--yaml', default='config.yaml',
                        help='Configuration yaml file.')
    parser.add_argument('--samplesheet', default='samplesheet.csv',
                        help='File with information of fastq files and other parameters for use.')
    


    # Parse command-line arguments
    args = parser.parse_args()

    # Create work and results dir for the project
    create_dirs(args=args)



    # Set up logger based on command-line arguments
    log_file = f"results/{args.project}/{args.project}.log"
    setup_logger(getattr(logging, args.verbose), log_file)
    logging.info(f'Log file created in {log_file}.')

    logger = logging.getLogger()

    # Check that yaml and samplesheet files exist:
    if not os.path.exists(f"projects/{args.project}/{args.yaml}"): 
        err = f'Path to yaml file projects/{args.project}/{args.yaml} does not exist.'
        logger.error(err); raise AssertionError(err)
    if not os.path.exists(f"projects/{args.project}/{args.samplesheet}"):
        err = f'Path to samplesheet file projects/{args.project}/{args.samplesheet} does not exist.'
        logger.error(err); raise AssertionError(err)


    # Parse config.yaml
    # TODO: COMPROBAR OUTPUTS
    yaml_dict, list_samplesheets, list_dbs_to_download = parse_yaml_file(yaml_path=f"projects/{args.project}/{args.yaml}", 
                                                                         project=args.project, 
                                                                         samplesheet=args.samplesheet)
       
    # Write yaml dict into work/PROJECT/config.yaml
    write_yaml_file(yaml_dict=yaml_dict,
                    yaml_path=f"work/{args.project}/config.yaml")


    # create_samplesheets(list_samplesheets)
    create_samplesheets(master_samplesheet_path=f"projects/{args.project}/{args.samplesheet}", 
                        project=args.project,
                        yaml_dict=yaml_dict, 
                        list_samplesheets=list_samplesheets)


    # create command.sh file
    # TODO: Los programas que haya que ejecutar tipo centrifuge, hacerlos mediante DOCKER
    create_command_sh_file(yaml_dict=yaml_dict, 
                           project=args.project)


    # download_DBs(list_dbs_to_download)
    database_download(list_dbs_to_download, project=args.project)
    os.system(f'work/{args.project}/database_download.sh')

    # run_pipeline(list_pipeline_commands)
    os.system(f'work/{args.project}/command.sh')




if __name__ == '__main__':
    main()