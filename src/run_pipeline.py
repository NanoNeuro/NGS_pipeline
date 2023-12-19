import argparse
import logging
import os

from processing_functions.logging import setup_logger
from processing_functions.generic_process import create_dirs
from processing_functions.yaml_config import parse_yaml_file

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

    

    # Parse config.yaml
    # TODO: COMPROBAR OUTPUTS
    yaml_dict, list_samplesheets, list_dbs_to_download = parse_yaml_file(yaml_path=f"projects/{args.project}/{args.yaml}", 
                                                                         project=args.project, 
                                                                         samplesheet=args.samplesheet)
    
    # make directories associated to the config
    # TODO ESCRIBIR config.yaml TEMPORAL EN WORK


    # create_samplesheets(list_samplesheets)


    # download_DBs(list_dbs_to_download)


    # create command.sh file


    # run_pipeline(list_pipeline_commands)





if __name__ == '__main__':
    main()